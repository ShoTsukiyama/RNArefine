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
from bpprediction import predict, predict2, dataloader, dataloader2, collate_fn2
from RNATransformer import Evoformer
from Net3 import Net
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
    
        
def predict_bp_from_backbonestr(dist_map, str_info, str_size, str_fe_1, str_fe_2, pdb_nts, pdb_nt_mask, seqidx2pdbidx, seqidx2nt, modelpath, chain_name, outpath):
    max_res_num = 16
    dist_threshold = 20
    masked_neighbor = 2
    
    B, S = str_size
    
    rots1, trans1, coord_mask1 = str_fe_1
    rots2, trans2, coord_mask2 = str_fe_2
    
    dist_map[torch.isnan(dist_map)] = torch.inf
    dist_map = torch.sqrt(dist_map).numpy()
    
    seqidx2segidx, flag = get_str_seg(dist_map, S)
    
    basepairs = []
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
            
            canonical_pair_map_bi = np.zeros((S, S))
            canonical_pair_map_bi[base_i, base_j] = 1
            canonical_pair_map_bi[base_j, base_i] = 1
            canonical_pair_vec_bi = canonical_pair_map_bi.sum(axis = 0)
            post_proc_flag = len(np.where(canonical_pair_vec_bi > 1)[0]) > 0
            
            if(post_proc_flag):
                canonical_pair_map = np.zeros((S, S))
                canonical_pair_map[base_i, base_j] = score
                canonical_pair_map[base_j, base_i] = score
                
                canonical_pair_map = apply_canonical_condition(canonical_pair_map)
                pair_idx = np.where(canonical_pair_map > 0)
                
                pair_idx = np.array(pair_idx).T
                keys = pd.DataFrame(pair_idx, columns=["base_i_seq", "base_j_seq"])

                canonical_pair_map_post_proc = canonical_pairs.merge(keys, on=["base_i_seq", "base_j_seq"], how="inner")
                
                canonical_pairs = canonical_pair_map_post_proc.loc[:, ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"]]
        
            else:
                canonical_pairs = canonical_pairs.loc[:, ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"]]
            
        if(len(non_canonical_pairs) > 0):
            non_canonical_pairs = non_canonical_pairs.loc[:, ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"]]
        
        basepairs = pd.concat([canonical_pairs, non_canonical_pairs])
        basepairs = basepairs.reset_index(drop = True)
        
        #pd.DataFrame(basepairs, columns = ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
        processflag = True
        
    else:
        if(np.sum(pdb_nt_mask).astype(int) > 1):
            print("Warning: No enough backbone atoms to predict base-pairs.")
            basepairs = pd.DataFrame([], columns = ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
            basepairmaps = []
        
        if(not flag):
            print("Warning: Too many sequence segments.")
            basepairs = pd.DataFrame([], columns = ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
            basepairmaps = []
            
    #os.makedirs(outpath, exist_ok = True)
    #basepairs.to_csv(f"{outpath}/RNAref3D_bp_bb.csv", index = False)
    #basepairs.to_csv(f"{outpath}/RNAref3D_bp.csv", index = False)
    
    #save_pickle(f"{outpath}/RNAref3D_bpmap_bb.pkl", basepairmaps)
        
    return basepairs, basepairmaps, processflag

def predict_bp_from_initbp(initial_bps, seq, dist_map, seqidx2pdbidx, seqidx2nt, modelpath, chain_name, outpath, maxlength=1024):
    L = len(seq)
    
    #dist_map_clip = dist_map.detach().clone()
    #dist_map_clip[torch.isinf(dist_map_clip)] = -1
    #dist_map_clip = torch.clamp(dist_map_clip, max=25.0)
    #dist_map_bin = torch.bucketize(dist_map_clip, torch.tensor([0, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]), right=True)
    
    pred_map = np.zeros([2, L, L])
    for l in range(len(initial_bps)):
        for j in range(1, len(initial_bps[l])):
            for k in range(len(initial_bps[l][j])):
                base_i = initial_bps[l][0][k][0]
                base_j = initial_bps[l][0][k][1]
                prob = initial_bps[l][j][k]
                
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
    model.load_state_dict(torch.load("{}/BP/post/model.pt".format(modelpath), map_location=torch.device(device)))
            
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
    













