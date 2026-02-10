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
from torch.utils.data import DataLoader
from Stack_layers.network import BPDetectNet
from Stack_trace_layers.network import BPDetectNet_trace
from stackprediction import predict, predict_trace, predict2, adjdataloader, nonadjdataloader, adjdataloader_trace, nonadjdataloader_trace, dataloader2, collate_fn2
from Net3 import Net
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

def get_bpmap(basepairs, pdbidx2seqidx, L):
    bpmap = torch.zeros((L+1, L+1))
    
    nt2pnt = {}
    if(len(basepairs) > 0):
        base_i = np.array(basepairs["base_i"].astype(int).map(pdbidx2seqidx).astype(int))
        base_j = np.array(basepairs["base_j"].astype(int).map(pdbidx2seqidx).astype(int))
        bpmap[base_i, base_j] = 1
        bpmap[base_j, base_i] = 1
        
        nt2pnt = dict(zip(base_i.tolist(), base_j.tolist()))
        nt2pnt_inv = dict(zip(base_j.tolist(), base_i.tolist()))
        
        nt2pnt = {**nt2pnt, **nt2pnt_inv}
    
    return bpmap, nt2pnt
    
        
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
        
def predict_adjstack_from_str(bp_map, nt2pnt, dist_map, str_info, str_size, str_fe_1, str_fe_2, str_fe_3, pdb_nts, pdb_nt_mask, seqidx2pdbidx, seqidx2nt, modelpath, chain_name, pred_type, outpath):
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
    
    stackpairs = []
    stackmaps = []
    if(np.sum(pdb_nt_mask).astype(int) > 1 and flag):
        adj_con_flag = [False] + [True if dist_map[0][j-1][j] < 8 else False for j in range(1, S)]
        
        target_pairs = [(0, j, j+1) for j in range(S-1)]
        
        samples = {"A-A": [], "A-U": [], "A-C": [], "A-G": [], "C-A": [], "C-U": [], "C-C": [], "C-G": [], "G-A": [], "G-U": [], "G-C": [], "G-G": [], "U-A": [], "U-U": [], "U-C": [], "U-G": []}
        
        for pairs in target_pairs:
            base_i = pairs[1] 
            base_j = pairs[2]
            
            if((not str_info[base_i]["required_flag"]) or (not str_info[base_j]["required_flag"])): continue
            
            base_nt_i = str_info[base_i]["nt_name"]
            base_nt_j = str_info[base_j]["nt_name"]
            
            base_ntidx_i = NTS2.index(base_nt_i)
            base_ntidx_j = NTS2.index(base_nt_j)
                
            samples[f"{NTS2[base_ntidx_i]}-{NTS2[base_ntidx_j]}"].append([base_i, base_j])
            
        for key, val in samples.items():
            if(len(val) == 0): continue
            
            model1 = BPDetectNet().to(device) 
            model1.load_state_dict(torch.load("{}/STACK/adj/{}/{}/model.pt".format(modelpath, "_".join(frame1_atoms), key), map_location=torch.device(device)))
            
            data_loader = adjdataloader(val, rots1, trans1, str_info, bp_map, nt2pnt, seqidx2segidx, seqidx2pdbidx, dist_map, pdb_nts, coord_mask1, adj_con_flag, max_res_num, S, device)
            loader = DataLoader(dataset = data_loader, batch_size = 64, shuffle=False)
        
            pred_cls = predict(model1)
            pred1 = pred_cls.predict(loader)
            pred1 = np.array(pred1)
            
            model2 = BPDetectNet().to(device) 
            model2.load_state_dict(torch.load("{}/STACK/adj/{}/{}/model.pt".format(modelpath, "_".join(frame2_atoms), key), map_location=torch.device(device)))
            
            data_loader = adjdataloader(val, rots2, trans2, str_info, bp_map, nt2pnt, seqidx2segidx, seqidx2pdbidx, dist_map, pdb_nts, coord_mask2, adj_con_flag, max_res_num, S, device)
            loader = DataLoader(dataset = data_loader, batch_size = 64, shuffle=False)
            
            pred_cls = predict(model2)
            pred2 = pred_cls.predict(loader)
            pred2 = np.array(pred2)
            
            if(pred_type == "full"):
                model3 = BPDetectNet().to(device) 
                model3.load_state_dict(torch.load("{}/STACK/adj/base/{}/model.pt".format(modelpath, key), map_location=torch.device(device)))

                data_loader = adjdataloader(val, rots3, trans3, str_info, bp_map, nt2pnt, seqidx2segidx, seqidx2pdbidx, dist_map, pdb_nts, coord_mask3, adj_con_flag, max_res_num, S, device)
                loader = DataLoader(dataset = data_loader, batch_size = 64, shuffle=False)
            
                pred_cls = predict(model3)
                pred3 = pred_cls.predict(loader)
                pred3 = np.array(pred3)
                
                pred = ((pred1 + pred2 + pred3)/3).tolist()
                stackpairs.extend([[chain_name, chain_name, seqidx2pdbidx[val[i][0]], seqidx2pdbidx[val[i][1]], seqidx2nt[val[i][0]], seqidx2nt[val[i][1]], pred[i]] for i in range(len(val)) if pred[i] > 0.5])
                stackmaps.append([val, pred1, pred2, pred3])
                
            else:
                pred = ((pred1 + pred2)/2).tolist()
                stackpairs.extend([[chain_name, chain_name, seqidx2pdbidx[val[i][0]], seqidx2pdbidx[val[i][1]], seqidx2nt[val[i][0]], seqidx2nt[val[i][1]], pred[i]] for i in range(len(val)) if pred[i] > 0.5])
                stackmaps.append([val, pred1, pred2])
            
        stackpairs = pd.DataFrame(stackpairs, columns = ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
            
        #os.makedirs(outpath, exist_ok = True)
        #stackpairs.to_csv(f"{outpath}/RNAref3D_stack_bb.csv", index = False)
        
        #save_pickle(f"{outpath}/RNAref3D_stackmap_bb.pkl", stackmaps)
        
    else:
        if(np.sum(pdb_nt_mask).astype(int) > 1):
            print("Warning: No enough backbone atoms to predict base-pairs.")
        
        if(not flag):
            print("Warning: Too many sequence segments.")
        
    return stackpairs, stackmaps

def predict_adjstack_from_C4trace(bp_map, nt2pnt, dist_map, str_info, str_size, str_fe, pdb_nts, pdb_nt_mask, seqidx2pdbidx, seqidx2nt, modelpath, chain_name):
    max_res_num = 16
    
    B, S = str_size
    
    trans, coord_mask = str_fe
    
    dist_map[torch.isnan(dist_map)] = torch.inf
    dist_map = torch.sqrt(dist_map).numpy()
    
    seqidx2segidx, flag = get_str_seg(dist_map, S)
    
    stackpairs = []
    stackmaps = []
    if(np.sum(pdb_nt_mask).astype(int) > 1 and flag):
        adj_con_flag = [False] + [True if dist_map[0][j-1][j] < 8 else False for j in range(1, S)]
        
        target_pairs = [(0, j, j+1) for j in range(S-1)]
        
        samples = {"A-A": [], "A-U": [], "A-C": [], "A-G": [], "C-A": [], "C-U": [], "C-C": [], "C-G": [], "G-A": [], "G-U": [], "G-C": [], "G-G": [], "U-A": [], "U-U": [], "U-C": [], "U-G": []}
        
        for pairs in target_pairs:
            base_i = pairs[1] 
            base_j = pairs[2]
            
            if((not str_info[base_i]["required_flag"]) or (not str_info[base_j]["required_flag"])): continue
            
            base_nt_i = str_info[base_i]["nt_name"]
            base_nt_j = str_info[base_j]["nt_name"]
            
            base_ntidx_i = NTS2.index(base_nt_i)
            base_ntidx_j = NTS2.index(base_nt_j)
                
            samples[f"{NTS2[base_ntidx_i]}-{NTS2[base_ntidx_j]}"].append([base_i, base_j])
            
        for key, val in samples.items():
            if(len(val) == 0): continue
            
            model = BPDetectNet_trace().to(device) 
            model.load_state_dict(torch.load("{}/STACK/adj/C4trace/{}/model.pt".format(modelpath, key), map_location=torch.device(device)))
            
            data_loader = adjdataloader_trace(val, trans, str_info, bp_map, nt2pnt, seqidx2segidx, seqidx2pdbidx, dist_map, pdb_nts, coord_mask, adj_con_flag, max_res_num, S, device)
            loader = DataLoader(dataset = data_loader, batch_size = 64, shuffle=False)
        
            pred_cls = predict_trace(model)
            pred = pred_cls.predict(loader)
            pred = np.array(pred).tolist()
            
            stackpairs.extend([[chain_name, chain_name, seqidx2pdbidx[val[i][0]], seqidx2pdbidx[val[i][1]], seqidx2nt[val[i][0]], seqidx2nt[val[i][1]], pred[i]] for i in range(len(val)) if pred[i] > 0.5])
            stackmaps.append([val, pred])
            
        stackpairs = pd.DataFrame(stackpairs, columns = ["chain_i", "chain_j", "base_i", "base_j", "Base_i_nt", "base_j_nt", "score"])
            
        #os.makedirs(outpath, exist_ok = True)
        #stackpairs.to_csv(f"{outpath}/RNAref3D_stack_bb.csv", index = False)
        
        #save_pickle(f"{outpath}/RNAref3D_stackmap_bb.pkl", stackmaps)
        
    else:
        if(np.sum(pdb_nt_mask).astype(int) > 1):
            print("Warning: No enough backbone atoms to predict base-pairs.")
        
        if(not flag):
            print("Warning: Too many sequence segments.")
        
    return stackpairs, stackmaps

def predict_nonadjstack_from_str(bp_map, nt2pnt, dist_map, str_info, str_size, str_fe_1, str_fe_2, str_fe_3, pdb_nts, pdb_nt_mask, seqidx2pdbidx, seqidx2nt, modelpath, chain_name, pred_type, outpath):
    max_res_num = 16
    dist_threshold = 15
    masked_neighbor = 2
    
    B, S = str_size
    
    rots1, trans1, coord_mask1 = str_fe_1
    rots2, trans2, coord_mask2 = str_fe_2
    rots3, trans3, coord_mask3 = str_fe_3
    
    dist_map[torch.isnan(dist_map)] = torch.inf
    dist_map = torch.sqrt(dist_map).numpy()
    
    seqidx2segidx, flag = get_str_seg(dist_map, S)
    
    stackpairs = []
    stackmaps = []
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
        
        target_pairs = [(0, j, j+1) for j in range(S-1)]
        
        idx = np.arange(S)
        adj_mask = (np.abs(idx[:, None] - idx) <= 1)[None, ...] #.expand(B, S, S)
        adj_mask = np.broadcast_to(adj_mask, (B, S, S))
        dist_nonadj = dist_map
        dist_nonadj[adj_mask] = np.inf
        target_pairs = np.array(np.where(dist_nonadj < dist_threshold)).T.tolist()
        target_pairs = [(target_pairs[j][0], target_pairs[j][1], target_pairs[j][2]) for j in range(len(target_pairs)) if target_pairs[j][1] < target_pairs[j][2]]
        
        samples = {"A-A": [], "A-U": [], "A-C": [], "A-G": [], "C-C": [], "G-C": [], "G-G": [], "U-C": [], "U-G": [], "U-U": []}
        
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
            
            model1 = BPDetectNet().to(device) 
            model1.load_state_dict(torch.load("{}/STACK/nonadj/{}/{}/model.pt".format(modelpath, "_".join(frame1_atoms), key), map_location=torch.device(device)))
            
            data_loader = nonadjdataloader(val, rots1, trans1, str_info, bp_map, nt2pnt, seqidx2segidx, seqidx2pdbidx, dist_map, pdb_nts, coord_mask1, adj_con_flag, max_res_num, S, device)
            loader = DataLoader(dataset = data_loader, batch_size = 64, shuffle=False)
        
            pred_cls = predict(model1)
            pred1 = pred_cls.predict(loader)
            pred1 = np.array(pred1)
            
            model2 = BPDetectNet().to(device) 
            model2.load_state_dict(torch.load("{}/STACK/nonadj/{}/{}/model.pt".format(modelpath, "_".join(frame2_atoms), key), map_location=torch.device(device)))
            
            data_loader = nonadjdataloader(val, rots2, trans2, str_info, bp_map, nt2pnt, seqidx2segidx, seqidx2pdbidx, dist_map, pdb_nts, coord_mask2, adj_con_flag, max_res_num, S, device)
            loader = DataLoader(dataset = data_loader, batch_size = 64, shuffle=False)
            
            pred_cls = predict(model2)
            pred2 = pred_cls.predict(loader)
            pred2 = np.array(pred2)
            
            if(pred_type == "full"):
                model3 = BPDetectNet().to(device) 
                model3.load_state_dict(torch.load("{}/STACK/nonadj/base/{}/model.pt".format(modelpath, key), map_location=torch.device(device)))

                data_loader = nonadjdataloader(val, rots3, trans3, str_info, bp_map, nt2pnt, seqidx2segidx, seqidx2pdbidx, dist_map, pdb_nts, coord_mask3, adj_con_flag, max_res_num, S, device)
                loader = DataLoader(dataset = data_loader, batch_size = 64, shuffle=False)
            
                pred_cls = predict(model3)
                pred3 = pred_cls.predict(loader)
                pred3 = np.array(pred3)
                
                pred = ((pred1 + pred2 + pred3)/3).tolist()
                stackpairs.extend([[chain_name, chain_name, seqidx2pdbidx[val[i][0]], seqidx2pdbidx[val[i][1]], seqidx2nt[val[i][0]], seqidx2nt[val[i][1]], pred[i]] for i in range(len(val)) if pred[i] > 0.5])
                stackmaps.append([val, pred1, pred2, pred3])
                
            else:
                pred = ((pred1 + pred2)/2).tolist()
                stackpairs.extend([[chain_name, chain_name, seqidx2pdbidx[val[i][0]], seqidx2pdbidx[val[i][1]], seqidx2nt[val[i][0]], seqidx2nt[val[i][1]], pred[i]] for i in range(len(val)) if pred[i] > 0.5])
                stackmaps.append([val, pred1, pred2])
            
        stackpairs = pd.DataFrame(stackpairs, columns = ["chain_i", "chain_j", "base_i", "base_j", "base_i_nt", "base_j_nt", "score"])
            
        #os.makedirs(outpath, exist_ok = True)
        #stackpairs.to_csv(f"{outpath}/RNAref3D_stack_bb.csv", index = False)
        
        #save_pickle(f"{outpath}/RNAref3D_stackmap_bb.pkl", stackmaps)
        
    else:
        if(np.sum(pdb_nt_mask).astype(int) > 1):
            print("Warning: No enough backbone atoms to predict base-pairs.")
        
        if(not flag):
            print("Warning: Too many sequence segments.")
        
    return stackpairs, stackmaps

def predict_nonadjstack_from_C4trace(bp_map, nt2pnt, dist_map, str_info, str_size, str_fe, pdb_nts, pdb_nt_mask, seqidx2pdbidx, seqidx2nt, modelpath, chain_name):
    max_res_num = 16
    dist_threshold = 15
    masked_neighbor = 2
    
    B, S = str_size
    
    trans, coord_mask = str_fe
    
    dist_map[torch.isnan(dist_map)] = torch.inf
    dist_map = torch.sqrt(dist_map).numpy()
    
    seqidx2segidx, flag = get_str_seg(dist_map, S)
    
    stackpairs = []
    stackmaps = []
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
        
        target_pairs = [(0, j, j+1) for j in range(S-1)]
        
        idx = np.arange(S)
        adj_mask = (np.abs(idx[:, None] - idx) <= 1)[None, ...] #.expand(B, S, S)
        adj_mask = np.broadcast_to(adj_mask, (B, S, S))
        dist_nonadj = dist_map
        dist_nonadj[adj_mask] = np.inf
        target_pairs = np.array(np.where(dist_nonadj < dist_threshold)).T.tolist()
        target_pairs = [(target_pairs[j][0], target_pairs[j][1], target_pairs[j][2]) for j in range(len(target_pairs)) if target_pairs[j][1] < target_pairs[j][2]]
        
        samples = {"A-A": [], "A-U": [], "A-C": [], "A-G": [], "C-C": [], "G-C": [], "G-G": [], "U-C": [], "U-G": [], "U-U": []}
        
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
            model.load_state_dict(torch.load("{}/STACK/nonadj/C4trace/{}/model.pt".format(modelpath, key), map_location=torch.device(device)))
            
            data_loader = nonadjdataloader_trace(val, trans, str_info, bp_map, nt2pnt, seqidx2segidx, seqidx2pdbidx, dist_map, pdb_nts, coord_mask, adj_con_flag, max_res_num, S, device)
            loader = DataLoader(dataset = data_loader, batch_size = 64, shuffle=False)
        
            pred_cls = predict_trace(model)
            pred = pred_cls.predict(loader)
            pred = np.array(pred).tolist()
            
            stackpairs.extend([[chain_name, chain_name, seqidx2pdbidx[val[i][0]], seqidx2pdbidx[val[i][1]], seqidx2nt[val[i][0]], seqidx2nt[val[i][1]], pred[i]] for i in range(len(val)) if pred[i] > 0.5])
            stackmaps.append([val, pred])
            
        stackpairs = pd.DataFrame(stackpairs, columns = ["chain_i", "chain_j", "base_i", "base_j", "base_i_nt", "base_j_nt", "score"])
            
        #os.makedirs(outpath, exist_ok = True)
        #stackpairs.to_csv(f"{outpath}/RNAref3D_stack_bb.csv", index = False)
        
        #save_pickle(f"{outpath}/RNAref3D_stackmap_bb.pkl", stackmaps)
        
    else:
        if(np.sum(pdb_nt_mask).astype(int) > 1):
            print("Warning: No enough backbone atoms to predict base-pairs.")
        
        if(not flag):
            print("Warning: Too many sequence segments.")
        
    return stackpairs, stackmaps





def predict_stack_from_initstack(initial_bps, seq, seqidx2pdbidx, seqidx2nt, modelpath, outpath, maxlength=1024):
    L = len(seq)
    
    pred_map = np.zeros([2, L, L])
    for l in range(len(initial_bps)):
        for j in range(1, len(initial_bps[l])):
            for k in range(len(initial_bps[l][j])):
                base_i = initial_bps[l][0][k][0]
                base_j = initial_bps[l][0][k][1]
                prob = initial_bps[l][j][k]
                
                pred_map[j-1][base_i][base_j] = int(prob > 0.5)
                pred_map[j-1][base_j][base_i] = int(prob > 0.5)
                
    pred_bp_list = np.array(np.where(pred_map > 0))
    
    data_loader = dataloader2([[pred_bp_list, seq]], maxlength=maxlength, device = device)
    loader = DataLoader(dataset = data_loader, collate_fn=lambda b: collate_fn2(b, maxlength=maxlength, device=device), batch_size = 16, shuffle=False)

    model = Net(dim_seq_fe=4, dim_hidden = 196).to(device) 
    model.load_state_dict(torch.load("{}/STACK/post/model.pt".format(modelpath), map_location=torch.device(device)))
            
    pred_cls = predict2(model)
    pred, probs = pred_cls.predict(loader)
    pred = pred[0]
    
    post_bp_idx_list = np.array(np.where(pred > 0))
    
    basepairs = []
    if(len(post_bp_idx_list) > 0):
        probs = np.array(probs[0][0])
        post_prob_list = probs[post_bp_idx_list[0], post_bp_idx_list[1]]
        post_bp_idx_list = post_bp_idx_list.T
        basepairs = [[seqidx2pdbidx[post_bp_idx_list[i][0]], seqidx2pdbidx[post_bp_idx_list[i][1]], seqidx2nt[post_bp_idx_list[i][0]], seqidx2nt[post_bp_idx_list[i][1]], post_prob_list[i]] for i in range(len(post_bp_idx_list)) if (seqidx2pdbidx[post_bp_idx_list[i][0]] < seqidx2pdbidx[post_bp_idx_list[i][1]])]

    basepairs = pd.DataFrame(basepairs, columns = ["base_i", "base_j", "base_i_nt", "base_j_nt", "score"])
    
    os.makedirs(outpath, exist_ok = True)
    basepairs.to_csv(f"{outpath}/RNAref3D_bp.csv", index = False)
    













