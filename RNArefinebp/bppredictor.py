#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 22:11:43 2024

@author: tsukiyamashou
"""
import os
import sys
from Bio import PDB
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pathlib import Path
import shutil
import argparse
import warnings
import string
import numpy as np
import pandas as pd
import torch
from pdbparser import PDBLoader
from affine import build_affine3d_from_coordinates
from torch.utils.data import DataLoader
from BP_layers.network import BPDetectNet
from BP_trace_layers.network import BPDetectNet_trace
from bpprediction import predict, predict_trace, predict2, dataloader, dataloader_trace, dataloader2, collate_fn2, parse_seq, get_seq_rep
from RNATransformer import Evoformer
from Net3 import Net
from numpy.lib.stride_tricks import sliding_window_view
import torch.nn.functional as F
import cvxpy as cp
import copy
import pickle

sys.path.append(os.path.join(os.path.dirname(__file__), './bp'))

#np.set_printoptions(threshold=np.inf)
#import npy2pdb

#sys.path.append(os.path.join(os.path.dirname(__file__), './bp/scripts'))
#import Feature

#sys.path.append(os.path.join(os.path.dirname(__file__), './bp/DeepE2EPotential'))
#import predict

rna_codes = {
    "A": "A", "C": "C", "G": "G", "U": "U", "T": "U",
    "DA": "A", "DC": "C", "DG": "G", "DT": "U", "DU": "U"
}

NTS = ["A","C","G","U","N"]
NTS2 = ["A", "U", "G", "C"]

frame1_atoms = ["P", "C4'", "C3'"]
frame2_atoms = ["O5'", "C3'", "O3'"]

device = "cuda" if torch.cuda.is_available() else "cpu"
        
def get_str_seg(dist_map, S):
    segment_idx = list(string.ascii_letters)
    current_idx = 0
    flag = True
    
    segments = [segment_idx[current_idx]]
    for i in range(1, S):
        if dist_map[0][i-1][i] < 8:
            segments.append(segment_idx[current_idx])
        else:
            current_idx += 1
            
            if(current_idx >= len(segment_idx)):
                flag=False
                break
                
            segments.append(segment_idx[current_idx])
    
    seqidx2segidx = dict(zip(np.arange(S).tolist(), segments))
    
    return seqidx2segidx, flag

def save_pickle(filename, data):
    with open(filename, 'wb') as f:
        pickle.dump(data, f)
        
def apply_canonical_condition(_X):
    """_summary_

    Args:
        _X (n,n): Base-pair probabilities

    Returns:
        optimal_X: Binary matrix satisfied secondary structure conditions
    """
    n = _X.shape[0]
    solver_name = "SCIP"
    A = _X
    B = copy.deepcopy(A)
    B[B>0] = 1
    B[B<=0] = 0
    B=B.astype(bool)
    X = cp.Variable((n, n), boolean=True)
    X.value = B

    constraints1 = [cp.sum(X[:, i]) <= 1 for i in range(n)]
    constraints2 = [X == X.T, X >= 0]

    obj = cp.Maximize(cp.sum(cp.multiply(A, X)))

    constraints = constraints1 + constraints2
    prob = cp.Problem(obj, constraints)
    result = prob.solve(solver_name)
    optimal_X = X.value

    return optimal_X
    
        
def predict_bp_from_backbonestr(dist_map, str_info, str_size, str_fe_1, str_fe_2, str_fe_3, pdb_nts, pdb_nt_mask, seqidx2pdbidx, seqidx2nt, modelpath, chain_name, pred_type, outpath):
    max_res_num = 16
    dist_threshold = 20
    masked_neighbor = 2
    
    B, S = str_size
    
    rots1, trans1, coord_mask1 = str_fe_1
    rots2, trans2, coord_mask2 = str_fe_2
    rots3, trans3, coord_mask3 = str_fe_3
    
    dist_map[torch.isnan(dist_map)] = torch.inf
    dist_map = torch.sqrt(dist_map).numpy()
    
    seqidx2segidx, flag = get_str_seg(dist_map, S)
    
    #basepairs = []
    basepairmaps = []
    processflag = False
    if(np.sum(pdb_nt_mask).astype(int) > 1 and flag):
        adj_con_flag = [False] + [True if dist_map[0][j-1][j] < 8 else False for j in range(1, S)]
        
        idx = np.arange(S)
        mask_neighbor = (np.abs(idx[:, None] - idx) <= masked_neighbor)[None, ...] #.expand(B, S, S)
        mask_neighbor = np.broadcast_to(mask_neighbor, (B, S, S))

        dist_map_candiate = dist_map.copy()
        dist_map_candiate[mask_neighbor] = np.inf
        target_pairs = np.array(np.where(dist_map_candiate < dist_threshold)).T
        
        target_pairs = [(target_pairs[j][0], target_pairs[j][1], target_pairs[j][2]) if target_pairs[j][1] < target_pairs[j][2] else (target_pairs[j][0], target_pairs[j][2], target_pairs[j][1]) for j in range(len(target_pairs))]
        target_pairs = list(set(target_pairs))
        
        samples = {"A-A": [], "A-U": [], "A-C": [], "A-G": [], "C-C": [], "G-C": [], "G-G": [], "U-C": [], "U-G": [], "U-U": []}
        basepairs = {"A-A": [], "A-U": [], "A-C": [], "A-G": [], "C-C": [], "G-C": [], "G-G": [], "U-C": [], "U-G": [], "U-U": []}
        basepairmaps = {"A-A": [], "A-U": [], "A-C": [], "A-G": [], "C-C": [], "G-C": [], "G-G": [], "U-C": [], "U-G": [], "U-U": []}
        
        for pairs in target_pairs:
            target_nucleotides = []
            base_i = pairs[1] 
            base_j = pairs[2]
            
            if((not str_info[base_i]["required_flag"]) or (not str_info[base_j]["required_flag"])): continue
            
            base_nt_i = str_info[base_i]["nt_name"]
            base_nt_j = str_info[base_j]["nt_name"]
            
            base_ntidx_i = NTS2.index(base_nt_i)
            base_ntidx_j = NTS2.index(base_nt_j)
            
            if(base_ntidx_i > base_ntidx_j):
                temp = base_i
                base_i = base_j
                base_j = temp
                
                temp = base_nt_i
                base_nt_i = base_nt_j
                base_nt_j = temp
                
                temp = base_ntidx_i
                base_ntidx_i = base_ntidx_j
                base_ntidx_j = temp
                
            samples[f"{NTS2[base_ntidx_i]}-{NTS2[base_ntidx_j]}"].append([base_i, base_j])
            
        for key, val in samples.items():
            if(len(val) == 0): continue
            
            model1 = BPDetectNet().to(device) 
            model1.load_state_dict(torch.load("{}/BP/{}/{}/model.pt".format(modelpath, "_".join(frame1_atoms), key), map_location=torch.device(device)))
            
            data_loader = dataloader(val, rots1, trans1, str_info, seqidx2segidx, seqidx2pdbidx, dist_map, pdb_nts, coord_mask1, adj_con_flag, max_res_num, S, device)
            loader = DataLoader(dataset = data_loader, batch_size = 64, shuffle=False)
        
            pred_cls = predict(model1)
            pred1 = pred_cls.predict(loader)
            pred1 = np.array(pred1)
            
            model2 = BPDetectNet().to(device) 
            model2.load_state_dict(torch.load("{}/BP/{}/{}/model.pt".format(modelpath, "_".join(frame2_atoms), key), map_location=torch.device(device)))
            
            data_loader = dataloader(val, rots2, trans2, str_info, seqidx2segidx, seqidx2pdbidx, dist_map, pdb_nts, coord_mask2, adj_con_flag, max_res_num, S, device)
            loader = DataLoader(dataset = data_loader, batch_size = 64, shuffle=False)
            
            pred_cls = predict(model2)
            pred2 = pred_cls.predict(loader)
            pred2 = np.array(pred2)
            
            
            if(pred_type == "full"):
                model3 = BPDetectNet().to(device) 
                model3.load_state_dict(torch.load("{}/BP/base/{}/model.pt".format(modelpath, key), map_location=torch.device(device)))
                
                data_loader = dataloader(val, rots3, trans3, str_info, seqidx2segidx, seqidx2pdbidx, dist_map, pdb_nts, coord_mask3, adj_con_flag, max_res_num, S, device)
                loader = DataLoader(dataset = data_loader, batch_size = 64, shuffle=False)
                
                pred_cls = predict(model3)
                pred3 = pred_cls.predict(loader)
                pred3 = np.array(pred3)

                pred = ((pred1 + pred2 + pred3)/3).tolist()
                basepairs[key] = [[chain_name, chain_name, val[i][0], val[i][1], seqidx2pdbidx[val[i][0]], seqidx2pdbidx[val[i][1]], seqidx2nt[val[i][0]], seqidx2nt[val[i][1]], pred[i]] for i in range(len(val)) if pred[i] > 0.5]
                basepairmaps[key] = [val, pred1, pred2, pred3]
            else:
                pred = ((pred1 + pred2)/2).tolist()
                basepairs[key] = [[chain_name, chain_name, val[i][0], val[i][1], seqidx2pdbidx[val[i][0]], seqidx2pdbidx[val[i][1]], seqidx2nt[val[i][0]], seqidx2nt[val[i][1]], pred[i]] for i in range(len(val)) if pred[i] > 0.5]
                basepairmaps[key] = [val, pred1, pred2]
            
            
        canonical_pairs = basepairs["A-U"] + basepairs["G-C"] + basepairs["U-G"]
        non_canonical_pairs = basepairs["A-A"] + basepairs["A-C"] + basepairs["A-G"] + basepairs["C-C"] + basepairs["G-G"] + basepairs["U-C"] + basepairs["U-U"]
        
        canonical_pairs = pd.DataFrame(canonical_pairs, columns = ["chain_i", "chain_j", "base_i_seq", "base_j_seq", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
        non_canonical_pairs = pd.DataFrame(non_canonical_pairs, columns = ["chain_i", "chain_j", "base_i_seq", "base_j_seq", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
        
        if(len(canonical_pairs) > 0):
            base_i = canonical_pairs.loc[:, "base_i_seq"].astype(int).tolist()
            base_j = canonical_pairs.loc[:, "base_j_seq"].astype(int).tolist()
            score = canonical_pairs.loc[:, "score"].astype(float).tolist()
            canonical_pair_map = np.zeros((S, S))
            canonical_pair_map[base_i, base_j] = score
            canonical_pair_map[base_j, base_i] = score
        else:
            canonical_pair_map = np.zeros((S, S))
        
        if(len(non_canonical_pairs) > 0):
            base_i = non_canonical_pairs.loc[:, "base_i_seq"].astype(int).tolist()
            base_j = non_canonical_pairs.loc[:, "base_j_seq"].astype(int).tolist()
            score = non_canonical_pairs.loc[:, "score"].astype(float).tolist()
            non_canonical_pair_map = np.zeros((S, S))
            non_canonical_pair_map[base_i, base_j] = score
            non_canonical_pair_map[base_j, base_i] = score
        else:
            non_canonical_pair_map = np.zeros((S, S))
        
        processflag = True
        
    else:
        if(np.sum(pdb_nt_mask).astype(int) > 1):
            print("Warning: No enough backbone atoms to predict base-pairs.")
            #basepairs = pd.DataFrame([], columns = ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
            basepairmaps = []
            canonical_pair_map = np.zeros((S, S))
            non_canonical_pair_map = np.zeros((S, S))
        
        if(not flag):
            print("Warning: Too many sequence segments.")
            #basepairs = pd.DataFrame([], columns = ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
            basepairmaps = []
            canonical_pair_map = np.zeros((S, S))
            non_canonical_pair_map = np.zeros((S, S))
            
    #os.makedirs(outpath, exist_ok = True)
    #basepairs.to_csv(f"{outpath}/RNAref3D_bp_bb.csv", index = False)
    #basepairs.to_csv(f"{outpath}/RNAref3D_bp.csv", index = False)
    
    #save_pickle(f"{outpath}/RNAref3D_bpmap_bb.pkl", basepairmaps)
        
    return basepairmaps, processflag, canonical_pair_map, non_canonical_pair_map
    #return basepairs, basepairmaps, processflag
    
def predict_bp_from_C4trace(dist_map, str_info, str_size, str_fe, pdb_nts, pdb_nt_mask, seqidx2pdbidx, seqidx2nt, modelpath, chain_name, pred_type, outpath):
    max_res_num = 16
    dist_threshold = 20
    masked_neighbor = 2
    
    B, S = str_size
    
    trans, coord_mask = str_fe
    
    dist_map[torch.isnan(dist_map)] = torch.inf
    dist_map = torch.sqrt(dist_map).numpy()
    
    seqidx2segidx, flag = get_str_seg(dist_map, S)
    
    #basepairs = []
    basepairmaps = []
    processflag = False
    if(np.sum(pdb_nt_mask).astype(int) > 1 and flag):
        adj_con_flag = [False] + [True if dist_map[0][j-1][j] < 8 else False for j in range(1, S)]
        
        idx = np.arange(S)
        mask_neighbor = (np.abs(idx[:, None] - idx) <= masked_neighbor)[None, ...] #.expand(B, S, S)
        mask_neighbor = np.broadcast_to(mask_neighbor, (B, S, S))

        dist_map_candiate = dist_map.copy()
        dist_map_candiate[mask_neighbor] = np.inf
        target_pairs = np.array(np.where(dist_map_candiate < dist_threshold)).T
        
        target_pairs = [(target_pairs[j][0], target_pairs[j][1], target_pairs[j][2]) if target_pairs[j][1] < target_pairs[j][2] else (target_pairs[j][0], target_pairs[j][2], target_pairs[j][1]) for j in range(len(target_pairs))]
        target_pairs = list(set(target_pairs))
        
        samples = {"A-A": [], "A-U": [], "A-C": [], "A-G": [], "C-C": [], "G-C": [], "G-G": [], "U-C": [], "U-G": [], "U-U": []}
        basepairs = {"A-A": [], "A-U": [], "A-C": [], "A-G": [], "C-C": [], "G-C": [], "G-G": [], "U-C": [], "U-G": [], "U-U": []}
        basepairmaps = {"A-A": [], "A-U": [], "A-C": [], "A-G": [], "C-C": [], "G-C": [], "G-G": [], "U-C": [], "U-G": [], "U-U": []}
        
        for pairs in target_pairs:
            base_i = pairs[1] 
            base_j = pairs[2]
            
            if((not str_info[base_i]["required_flag"]) or (not str_info[base_j]["required_flag"])): continue
            
            base_nt_i = str_info[base_i]["nt_name"]
            base_nt_j = str_info[base_j]["nt_name"]
            
            base_ntidx_i = NTS2.index(base_nt_i)
            base_ntidx_j = NTS2.index(base_nt_j)
            
            if(base_ntidx_i > base_ntidx_j):
                temp = base_i
                base_i = base_j
                base_j = temp
                
                temp = base_nt_i
                base_nt_i = base_nt_j
                base_nt_j = temp
                
                temp = base_ntidx_i
                base_ntidx_i = base_ntidx_j
                base_ntidx_j = temp
                
            samples[f"{NTS2[base_ntidx_i]}-{NTS2[base_ntidx_j]}"].append([base_i, base_j])
            
        for key, val in samples.items():
            if(len(val) == 0): continue
            
            model = BPDetectNet_trace().to(device) 
            model.load_state_dict(torch.load("{}/BP/C4trace/{}/model.pt".format(modelpath, key), map_location=torch.device(device)))
            
            data_loader = dataloader_trace(val, trans, str_info, seqidx2segidx, seqidx2pdbidx, dist_map, pdb_nts, coord_mask, adj_con_flag, max_res_num, S, device)
            loader = DataLoader(dataset = data_loader, batch_size = 64, shuffle=False)
        
            pred_cls = predict_trace(model)
            pred = pred_cls.predict(loader)
            pred = np.array(pred).tolist()
            
            basepairs[key] = [[chain_name, chain_name, val[i][0], val[i][1], seqidx2pdbidx[val[i][0]], seqidx2pdbidx[val[i][1]], seqidx2nt[val[i][0]], seqidx2nt[val[i][1]], pred[i]] for i in range(len(val)) if pred[i] > 0.5]
            basepairmaps[key] = [val, pred]
            
            
        canonical_pairs = basepairs["A-U"] + basepairs["G-C"] + basepairs["U-G"]
        non_canonical_pairs = basepairs["A-A"] + basepairs["A-C"] + basepairs["A-G"] + basepairs["C-C"] + basepairs["G-G"] + basepairs["U-C"] + basepairs["U-U"]
        
        canonical_pairs = pd.DataFrame(canonical_pairs, columns = ["chain_i", "chain_j", "base_i_seq", "base_j_seq", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
        non_canonical_pairs = pd.DataFrame(non_canonical_pairs, columns = ["chain_i", "chain_j", "base_i_seq", "base_j_seq", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
        
        if(len(canonical_pairs) > 0):
            base_i = canonical_pairs.loc[:, "base_i_seq"].astype(int).tolist()
            base_j = canonical_pairs.loc[:, "base_j_seq"].astype(int).tolist()
            score = canonical_pairs.loc[:, "score"].astype(float).tolist()
            canonical_pair_map = np.zeros((S, S))
            canonical_pair_map[base_i, base_j] = score
            canonical_pair_map[base_j, base_i] = score
        else:
            canonical_pair_map = np.zeros((S, S))
        
        if(len(non_canonical_pairs) > 0):
            base_i = non_canonical_pairs.loc[:, "base_i_seq"].astype(int).tolist()
            base_j = non_canonical_pairs.loc[:, "base_j_seq"].astype(int).tolist()
            score = non_canonical_pairs.loc[:, "score"].astype(float).tolist()
            non_canonical_pair_map = np.zeros((S, S))
            non_canonical_pair_map[base_i, base_j] = score
            non_canonical_pair_map[base_j, base_i] = score
        else:
            non_canonical_pair_map = np.zeros((S, S))
        
        processflag = True
        
    else:
        if(np.sum(pdb_nt_mask).astype(int) > 1):
            print("Warning: No enough backbone atoms to predict base-pairs.")
            #basepairs = pd.DataFrame([], columns = ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
            basepairmaps = []
            canonical_pair_map = np.zeros((S, S))
            non_canonical_pair_map = np.zeros((S, S))
        
        if(not flag):
            print("Warning: Too many sequence segments.")
            #basepairs = pd.DataFrame([], columns = ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
            basepairmaps = []
            canonical_pair_map = np.zeros((S, S))
            non_canonical_pair_map = np.zeros((S, S))
            
    #os.makedirs(outpath, exist_ok = True)
    #basepairs.to_csv(f"{outpath}/RNAref3D_bp_bb.csv", index = False)
    #basepairs.to_csv(f"{outpath}/RNAref3D_bp.csv", index = False)
    
    #save_pickle(f"{outpath}/RNAref3D_bpmap_bb.pkl", basepairmaps)
        
    return basepairmaps, processflag, canonical_pair_map, non_canonical_pair_map
    

def predict_bp_from_initbp(initial_bps, seq, dist_map, modelpath, pred_type, maxlength=1024):
    L = len(seq)
    
    #dist_map_clip = dist_map.detach().clone()
    #dist_map_clip[torch.isinf(dist_map_clip)] = -1
    #dist_map_clip = torch.clamp(dist_map_clip, max=25.0)
    #dist_map_bin = torch.bucketize(dist_map_clip, torch.tensor([0, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]), right=True)
    
    if(pred_type == "full"): pred_map = np.zeros([3, L, L])
    else: pred_map = np.zeros([2, L, L])
    
    for key in ["A-A", "A-U", "A-C", "A-G", "C-C", "G-C", "G-G", "U-C", "U-G", "U-U"]:
    #for l in range(len(initial_bps)):
        for j in range(1, len(initial_bps[key])):
            for k in range(len(initial_bps[key][j])):
                base_i = initial_bps[key][0][k][0]
                base_j = initial_bps[key][0][k][1]
                prob = initial_bps[key][j][k]
                
                #pred_map[j-1][base_i][base_j] = int(prob > 0.5)
                #pred_map[j-1][base_j][base_i] = int(prob > 0.5)
                
                pred_map[j-1][base_i][base_j] += prob
                pred_map[j-1][base_j][base_i] += prob
              
    pred_bp_map = pred_map.mean(0)
    pred_bp_list = np.array(np.where(pred_bp_map > 0.5))
    
    data_loader = dataloader2([[pred_bp_list, seq, dist_map]], maxlength=maxlength, device = device)
    loader = DataLoader(dataset = data_loader, collate_fn=lambda b: collate_fn2(b, maxlength=maxlength, device=device), batch_size = 1, shuffle=False)

    #data_loader = dataloader2([[pred_bp_list, seq, dist_map_bin]], maxlength=maxlength, device = device)
    #loader = DataLoader(dataset = data_loader, batch_size = 1, shuffle=False)

    model = Net(dim_hidden = 96).to(device) 
    #model = Evoformer(m_dim = 128, z_dim = 64).to(device)
    model.load_state_dict(torch.load("{}/BP2/model.pt".format(modelpath), map_location=torch.device(device)))
            
    pred_cls = predict2(model)
    pred, probs = pred_cls.predict(loader, [L])
    
    """
    prob = np.concatenate([pred_map, probs[0][None, :, :]])
    prob = prob.mean(0)
    
    predicts = (prob >= 0.5)
    post_bp_idx_list = np.array(np.where(predicts > 0))
    
    basepairs = []
    if(len(post_bp_idx_list) > 0):
        post_prob_list = prob[post_bp_idx_list[0], post_bp_idx_list[1]]
        post_bp_idx_list = post_bp_idx_list.T
        basepairs = [[seqidx2pdbidx[post_bp_idx_list[i][0]], seqidx2pdbidx[post_bp_idx_list[i][1]], seqidx2nt[post_bp_idx_list[i][0]], seqidx2nt[post_bp_idx_list[i][1]], post_prob_list[i]] for i in range(len(post_bp_idx_list)) if (seqidx2pdbidx[post_bp_idx_list[i][0]] < seqidx2pdbidx[post_bp_idx_list[i][1]])]
    """
    #print(probs)
    
    aa_type = parse_seq(seq)
    #A:0, G:1, C:2, U:3
    onehot_=torch.from_numpy(aa_type)
    onehot=F.one_hot(onehot_.long(),6).float()[:, 1:5]
    nn_pairs = torch.einsum('ia,jb->abij', onehot, onehot)
    nn_pairs = nn_pairs.reshape(16, L, L)
    #AA:0, AG:1, AC:2, AU(T):3, GA:4, GG:5, GC:6, GU(T):7, CA:8, CG:9, CC:10, CU(T):11, UA:12, UG:13, UC:14, UU(T):15
    
    nt_dict = {"A":0, "G":1, "C":2, "U":3}
    canonical_pairs = [["A", "U"], ["G", "C"], ["U", "G"]]
    canonical_map = torch.zeros(L,L)
    for nt1, nt2 in canonical_pairs:
        canonical_map += nn_pairs[nt_dict[nt1]+4*nt_dict[nt2]]
        canonical_map += nn_pairs[nt_dict[nt2]+4*nt_dict[nt1]]
    
    rem_onehot=F.one_hot(onehot_.long(),6).float()[:, [0,5]] #padding and non-standard nn
    rem_pairs = torch.einsum('ia,jb->abij', rem_onehot, rem_onehot)
    rem_pairs = rem_pairs.reshape(4, L, L).sum(0)
    
    non_canonical_map = torch.ones(L,L) - canonical_map - rem_pairs
    
    canonical_map = canonical_map.detach().clone().numpy()
    non_canonical_map = non_canonical_map.detach().clone().numpy()
    probs = np.array(probs[0])
    
    probs[probs <= 0.5] = 0.0
    
    return probs*canonical_map, probs*non_canonical_map
    
    """
    print(len(probs))
    print(probs[0])
    print(canonical_map.size())
    print(non_canonical_map.size())
    print((canonical_map + non_canonical_map).size())
    
    sys.exit()
    
    #canonical_pairs = basepairs["A-U"] + basepairs["G-C"] + basepairs["U-G"]
    #non_canonical_pairs = basepairs["A-A"] + basepairs["A-C"] + basepairs["A-G"] + basepairs["C-C"] + basepairs["G-G"] + basepairs["U-C"] + basepairs["U-U"]
    
    
    
    print(onehot)
    print(nn_pairs[0])
    print(nn_pairs[1])
    print(nn_pairs[2])
    print(nn_pairs[3])
    print(nn_pairs.sum(0))
    print(onehot.size())
    print(nn_pairs.size())
    sys.exit()
    
    basepairmaps = {"A-A": [], "A-U": [], "A-C": [], "A-G": [], "C-C": [], "G-C": [], "G-G": [], "U-C": [], "U-G": [], "U-U": []}
    
    pred = pred[0]
    post_bp_idx_list = np.array(np.where(pred > 0))
    basepairs = []
    if(len(post_bp_idx_list) > 0):
        probs = np.array(probs[0])
        
        post_prob_list = probs[post_bp_idx_list[0], post_bp_idx_list[1]]
        post_bp_idx_list = post_bp_idx_list.T
        basepairs = [[chain_name, chain_name, seqidx2pdbidx[post_bp_idx_list[i][0]], seqidx2pdbidx[post_bp_idx_list[i][1]], seqidx2nt[post_bp_idx_list[i][0]], seqidx2nt[post_bp_idx_list[i][1]], post_prob_list[i]] for i in range(len(post_bp_idx_list)) if (seqidx2pdbidx[post_bp_idx_list[i][0]] < seqidx2pdbidx[post_bp_idx_list[i][1]])]
    
    basepairs = pd.DataFrame(basepairs, columns = ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
    
    basepairs.loc[:, "chain_i"] = basepairs.loc[:, "chain_i"].astype(str)
    basepairs.loc[:, "chain_j"] = basepairs.loc[:, "chain_j"].astype(str)
    basepairs.loc[:, "base_i"] = basepairs.loc[:, "base_i"].astype(str)
    basepairs.loc[:, "base_j"] = basepairs.loc[:, "base_j"].astype(str)
    basepairs.loc[:, "score"] = basepairs.loc[:, "score"].astype(str)
    interact = (basepairs.loc[:, ["chain_i", "chain_j", "base_i", "base_j", "score"]]).values.tolist()

    interact = [",".join(interact[i]) for i in range(len(interact))]

    with open(f"{outpath}/RNAref3D_bp.txt", "w") as f:
        f.write("\n".join(interact))
    
    #os.makedirs(outpath, exist_ok = True)
    #basepairs.to_csv(f"{outpath}/RNAref3D_bp.csv", index = False)
    
    return basepairs
    """
    
def transform_bp(initial_bps, frame_type_num, S):
    pred_map = np.zeros([frame_type_num, S, S])
    
    for key in ["A-A", "A-U", "A-C", "A-G", "C-C", "G-C", "G-G", "U-C", "U-G", "U-U"]:
        for j in range(1, len(initial_bps[key])):
            for k in range(len(initial_bps[key][j])):
                base_i = initial_bps[key][0][k][0]
                base_j = initial_bps[key][0][k][1]
                prob = initial_bps[key][j][k]
                
                pred_map[j-1][base_i][base_j] = prob
                pred_map[j-1][base_j][base_i] = prob
                
    return pred_map
    
def predict_bp_ml(initial_bps_, dist_map, str_info, str_size, seqidx2pdbidx, seqidx2nt, modelpath, chain_name, pred_type, c_trace_flag):
    basepairs = {"A-A": [], "A-C": [], "A-G": [], "A-U": [], "C-C": [], "C-G": [], "C-U": [], "G-G": [], "G-U": [], "U-U": []}
    basepairmaps = {"A-A": [], "A-C": [], "A-G": [], "A-U": [], "C-C": [], "C-G": [], "C-U": [], "G-G": [], "G-U": [], "U-U": []}
    NTS_ML = ["A","C","G","U"]
        
    B, S = str_size
    
    dist_threshold = 20
    masked_neighbor = 2
    max_dist = 30
    window_size = 4
    patch_len = window_size * 2 + 1
        
    idx = np.arange(S)
    mask = (np.abs(idx[:, None] - idx) <= masked_neighbor)[None, ...] #.expand(B, S, S)
    mask = np.broadcast_to(mask, (B, S, S))
    
    dist_nonneigh = np.sqrt(dist_map).numpy()
    dist_nonneigh[mask] = np.inf
    
    pdb_idx = [str_info[j]["seq_idx"] for j in range(len(str_info))]
    seq = np.array([NTS_ML.index(str_info[j]["nt_name"]) if str_info[j]["nt_name"] in NTS_ML else len(NTS_ML)-1 for j in range(len(str_info))])
    seq_pair = seq[:, None] * 4 + seq[None, :]
            
    #pdbidx2seqidx = dict(zip(pdb_idx, np.arange(S).tolist()))
    #seqidx2pdbidx = dict(zip(np.arange(S).tolist(), pdb_idx))
            
    if(c_trace_flag):
        initial_bps = transform_bp(initial_bps_, 1, S)
        model_name = "C4trace"
    elif(pred_type == "full"):
        initial_bps = transform_bp(initial_bps_, 3, S)
        model_name = "bb_base"
    else:
        initial_bps = transform_bp(initial_bps_, 2, S)
        model_name = "bb"
        
    initial_bps = initial_bps.mean(0)
    
    initial_bps[initial_bps > 0.5] = 1.0
    initial_bps[initial_bps <= 0.5] = 0.0
    
    dist_map_fe = dist_nonneigh[0]
    dist_map_fe = np.nan_to_num(dist_map_fe, nan=-1, posinf=-1, neginf=-1)
    dist_map_fe[dist_map_fe > max_dist] = max_dist
        
    dist_map_windows = sliding_window_view(
        np.pad(dist_map_fe, window_size, mode='constant', constant_values=-1), 
        window_shape=(patch_len, patch_len)
        )
    
    seq_pair_windows = sliding_window_view(
        np.pad(seq_pair, window_size, mode='constant', constant_values=-1),
        window_shape=(patch_len, patch_len)
        )
    
    initial_bps_windows = sliding_window_view(
        np.pad(initial_bps, window_size, mode='constant', constant_values=-1), 
        window_shape=(patch_len, patch_len)
        )        
    
    lower_tri_mask = np.tril(np.ones((S,S), dtype=int), k=0)
    
    is_dist_ok  = (dist_nonneigh[0] < dist_threshold)
    
    target_pair = np.array(np.where(lower_tri_mask & is_dist_ok)).T

    if(len(target_pair) > 0):
        dist_map_fe = dist_map_windows[target_pair[:, 0], target_pair[:, 1]].reshape(len(target_pair), -1)
        seq_pair_fe = seq_pair_windows[target_pair[:, 0], target_pair[:, 1]].reshape(len(target_pair), -1)
        initial_bps = initial_bps_windows[target_pair[:, 0], target_pair[:, 1]].reshape(len(target_pair), -1)
                   
        sample_fe = np.concatenate([dist_map_fe, seq_pair_fe, initial_bps], axis = 1)     
        
        seq_pair_type = seq_pair[target_pair[:, 0], target_pair[:, 1]]
        
        for k in range(len(NTS_ML)):
            for l in range(k, len(NTS_ML)):
                target_flag = ((seq_pair_type ==4*k+l) | (seq_pair_type ==4*l+k))
                            
                #if(len(sample_fe[target_flag]) > 0):
                    #print(sample_fe[target_flag].shape)
                    
                #ml_model = pickle.load(
                    #"{}/BPML2/{}/{}-{}/model.pt".format(modelpath, model_name, NTS_ML[k], NTS_ML[l])
                #)
                
                if np.sum(target_flag) <= 0: continue
                
                with open("{}/BPML2/LGBM/all/{}-{}/model.pt".format(modelpath, NTS_ML[k], NTS_ML[l]), 'rb') as f:
                    ml_model = pickle.load(f)
                
                pred = ml_model.predict(sample_fe[target_flag])
                
                #pred_proba = ml_model.predict_proba(sample_fe[target_flag])
                #pred = pred_proba[:, 1]
                                
                val = target_pair[target_flag, :].tolist()
                
                basepairs[f"{NTS_ML[k]}-{NTS_ML[l]}"] = [[chain_name, chain_name, val[i][0], val[i][1], seqidx2pdbidx[val[i][0]], seqidx2pdbidx[val[i][1]], seqidx2nt[val[i][0]], seqidx2nt[val[i][1]], pred[i]] for i in range(len(val)) if pred[i] > 0.5]
                basepairmaps[f"{NTS_ML[k]}-{NTS_ML[l]}"] = [val, pred.tolist()]
        
        canonical_pairs = basepairs["A-U"] + basepairs["C-G"] + basepairs["G-U"]
        non_canonical_pairs = basepairs["A-A"] + basepairs["A-C"] + basepairs["A-G"] + basepairs["C-C"] + basepairs["G-G"] + basepairs["C-U"] + basepairs["U-U"]
        
        canonical_pairs = pd.DataFrame(canonical_pairs, columns = ["chain_i", "chain_j", "base_i_seq", "base_j_seq", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
        non_canonical_pairs = pd.DataFrame(non_canonical_pairs, columns = ["chain_i", "chain_j", "base_i_seq", "base_j_seq", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
        
        if(len(canonical_pairs) > 0):
            base_i = canonical_pairs.loc[:, "base_i_seq"].astype(int).tolist()
            base_j = canonical_pairs.loc[:, "base_j_seq"].astype(int).tolist()
            score = canonical_pairs.loc[:, "score"].astype(float).tolist()
            canonical_pair_map = np.zeros((S, S))
            canonical_pair_map[base_i, base_j] = score
            canonical_pair_map[base_j, base_i] = score
        else:
            canonical_pair_map = np.zeros((S, S))
            
        if(len(non_canonical_pairs) > 0):
            base_i = non_canonical_pairs.loc[:, "base_i_seq"].astype(int).tolist()
            base_j = non_canonical_pairs.loc[:, "base_j_seq"].astype(int).tolist()
            score = non_canonical_pairs.loc[:, "score"].astype(float).tolist()
            non_canonical_pair_map = np.zeros((S, S))
            non_canonical_pair_map[base_i, base_j] = score
            non_canonical_pair_map[base_j, base_i] = score
        else:
            non_canonical_pair_map = np.zeros((S, S))
        
        processflag = True
        
    else:
        print("Warning: No nucleotide pairs to predict base-pairs.")
        #basepairs = pd.DataFrame([], columns = ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
        basepairmaps = []
        canonical_pair_map = np.zeros((S, S))
        non_canonical_pair_map = np.zeros((S, S))
    
    return canonical_pair_map, non_canonical_pair_map
    
    
def apply_bp_canonical_constraint(canonical_pair_map):
    post_proc_flag = len(np.where(canonical_pair_map > 0)[0]) > 0
    
    if(post_proc_flag):
        canonical_pair_map_bi = apply_canonical_condition(canonical_pair_map)
    
        return canonical_pair_map*canonical_pair_map_bi

    return canonical_pair_map

def get_basepair_df(chain_name, canonical_pair_map, non_canonical_pair_map, seqidx2pdbidx, seqidx2nt):
    canonical_pair_list = np.array(np.where(canonical_pair_map > 0))
    
    if(len(canonical_pair_list) > 0):
        canonical_pair_prob_list = canonical_pair_map[canonical_pair_list[0], canonical_pair_list[1]]
        
    canonical_pair_list = canonical_pair_list.T
    canonical_basepairs = [[chain_name, chain_name, seqidx2pdbidx[canonical_pair_list[i][0]], seqidx2pdbidx[canonical_pair_list[i][1]], seqidx2nt[canonical_pair_list[i][0]], seqidx2nt[canonical_pair_list[i][1]], canonical_pair_prob_list[i]] for i in range(len(canonical_pair_list)) if (seqidx2pdbidx[canonical_pair_list[i][0]] < seqidx2pdbidx[canonical_pair_list[i][1]])]
    canonical_basepairs = pd.DataFrame(canonical_basepairs, columns = ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
    
    non_canonical_pair_list = np.array(np.where(non_canonical_pair_map > 0))
    if(len(non_canonical_pair_list) > 0):
        non_canonical_pair_prob_list = non_canonical_pair_map[non_canonical_pair_list[0], non_canonical_pair_list[1]]
        
    non_canonical_pair_list = non_canonical_pair_list.T
    non_canonical_basepairs = [[chain_name, chain_name, seqidx2pdbidx[non_canonical_pair_list[i][0]], seqidx2pdbidx[non_canonical_pair_list[i][1]], seqidx2nt[non_canonical_pair_list[i][0]], seqidx2nt[non_canonical_pair_list[i][1]], non_canonical_pair_prob_list[i]] for i in range(len(non_canonical_pair_list)) if (seqidx2pdbidx[non_canonical_pair_list[i][0]] < seqidx2pdbidx[non_canonical_pair_list[i][1]])]
    non_canonical_basepairs = pd.DataFrame(non_canonical_basepairs, columns = ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
    
    basepairs = pd.concat([canonical_basepairs, non_canonical_basepairs])
    basepairs = basepairs.reset_index(drop = True)
    
    return basepairs

def save_bpinfo(outpath, basepairs):
    os.makedirs(outpath, exist_ok = True)
    basepairs.to_csv(f"{outpath}/RNAref3D_bp.csv", index = False)
    
