#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 17:45:09 2023

@author: tsukiyamashou
"""
import os
import sys
import pandas as pd
import numpy as np
import torch
import torch.utils.data as data
from sklearn.model_selection import train_test_split
import torch
from torch import nn
import torch.nn.functional as F
#from RCLM.Model import RNA2nd
from featurizer import get_paring_feature, get_paring_outerprot
import warnings
warnings.filterwarnings("ignore")

pd.set_option('display.max_rows', None)        # 行数制限を解除
pd.set_option('display.max_columns', None)     # 列数制限を解除
pd.set_option('display.width', None)           # 自動折り返し無効
pd.set_option('display.max_colwidth', None)    # 列内の文字列の最大幅を解除

NTS = ["A","C","G","U","N"]
NTS2 = ["A", "U", "G", "C"]

class dataloader(data.Dataset):
    def __init__(self, pairs, rots, trans, str_info, seqidx2segidx, seqidx2pdbidx, dist_map, pdb_nts, mask, adj_con_flag, max_res_num, seq_len, device):
                 #structure_info, max_res_num=200, device = "cpu"):
        super().__init__()
        self.pairs=pairs
        self.rots = rots
        self.trans = trans
        self.str_info=str_info
        self.seqidx2segidx = seqidx2segidx
        self.seqidx2pdbidx = seqidx2pdbidx
        self.dist_map = dist_map
        self.adj_con_flag = adj_con_flag
        self.max_res_num = max_res_num
        self.pdb_nts = pdb_nts
        self.mask = mask
        self.S = seq_len
        self.device = device
        
    def __len__(self):
        return len(self.pairs)
    
    def __getitem__(self, idx):
        base_i = self.pairs[idx][0] 
        base_j = self.pairs[idx][1]
        
        base_nt_i = self.str_info[base_i]["nt_name"]
        base_nt_j = self.str_info[base_j]["nt_name"]
        
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
            
        base_seg_i = self.seqidx2segidx[base_i]
        base_seg_j = self.seqidx2segidx[base_j]
        
        base_pdb_i = self.seqidx2pdbidx[base_i]
        base_pdb_j = self.seqidx2pdbidx[base_j]
        
        target_nucleotides = []
        target_nucleotides.append(base_i)
        target_nucleotides.append(base_j)
        
        dist_ij = np.concatenate([self.dist_map[0, base_i, :].reshape([1,-1]), self.dist_map[0, base_j, :].reshape([1,-1])], axis = 0).T
        dist_ij = np.min(dist_ij, axis = 1) 
        dist_ij[base_i] = np.inf
        dist_ij[base_j] = np.inf 
        
        if(base_i-1 >= 0): 
            if(self.adj_con_flag[base_i]):
                base_pdb = base_i-1
                if(self.str_info[base_i-1]["required_flag"] and base_pdb not in target_nucleotides):
                    target_nucleotides.append(base_pdb)
                    dist_ij[base_pdb] = np.inf
                
                if(base_i-2 >= 0):
                    if(self.adj_con_flag[base_i-1]):
                        base_pdb = base_i-2
                        if(self.str_info[base_i-2]["required_flag"] and base_pdb not in target_nucleotides):
                            target_nucleotides.append(base_pdb)
                            dist_ij[base_i-2] = np.inf
                            
        if(base_i+1 < self.S):
            if(self.adj_con_flag[base_i+1]):
                base_pdb = base_i+1
                if(self.str_info[base_i+1]["required_flag"] and base_pdb not in target_nucleotides):
                    target_nucleotides.append(base_pdb)
                    dist_ij[base_i+1] = np.inf
                
                if(base_i+2 < self.S):
                    if(self.adj_con_flag[base_i+2]):
                        base_pdb = base_i+2
                        if(self.str_info[base_i+2]["required_flag"] and base_pdb not in target_nucleotides):
                            target_nucleotides.append(base_pdb)
                            dist_ij[base_i+2] = np.inf
                        
        if(base_j-1 >= 0): 
            if(self.adj_con_flag[base_j]):
                base_pdb = base_j-1
                if(self.str_info[base_j-1]["required_flag"] and base_pdb not in target_nucleotides):
                    target_nucleotides.append(base_pdb)
                    dist_ij[base_j-1] = np.inf
                
                if(base_j-2 >= 0):
                    if(self.adj_con_flag[base_j-1]):
                        base_pdb = base_j-2
                        if(self.str_info[base_j-2]["required_flag"] and base_pdb not in target_nucleotides):
                            target_nucleotides.append(base_pdb)
                            dist_ij[base_j-2] = np.inf
            
        if(base_j+1 < self.S):
            if(self.adj_con_flag[base_j+1]):
                base_pdb = base_j+1
                if(self.str_info[base_j+1]["required_flag"] and base_pdb not in target_nucleotides):
                    target_nucleotides.append(base_pdb)
                    dist_ij[base_j+1] = np.inf
                
                if(base_j+2 < self.S):
                    if(self.adj_con_flag[base_j+2]):
                        base_pdb = base_j+2
                        if(self.str_info[base_j+2]["required_flag"] and base_pdb not in target_nucleotides):
                            target_nucleotides.append(base_pdb)
                            dist_ij[base_j+2] = np.inf
                                  
        dist_ij_idx = np.argsort(dist_ij)
        dist_ij = dist_ij[dist_ij_idx]
        dist_ij_idx = dist_ij_idx[~np.isinf(dist_ij)]  
        
        spatial_neighbor_num = max(self.max_res_num - len(target_nucleotides), 0)   
        spatial_neighbors = []         
        if(spatial_neighbor_num > 0): spatial_neighbors = dist_ij_idx[:spatial_neighbor_num]
        
        if(len(spatial_neighbors) > 0): target_nucleotides.extend(spatial_neighbors)   
        
        #target_nucleotides_seg = np.array([seqidx2segidx[pdbidx2seqidx[j]] for j in target_nucleotides]) 
        
        pad = [np.nan for j in range(self.max_res_num - len(target_nucleotides))]
        relpos_i = np.array([self.seqidx2pdbidx[j] - base_pdb_i if(self.seqidx2segidx[j] == base_seg_i) else 1000000 for j in target_nucleotides] + pad)  
        relpos_j = np.array([self.seqidx2pdbidx[j] - base_pdb_j if(self.seqidx2segidx[j] == base_seg_j) else 1000000 for j in target_nucleotides] + pad)
        #relpos_i[np.isnan(relpos_i)] = -1000000
        #relpos_j[np.isnan(relpos_j)] = -1000000
        
        relpos_idx_rel_i = torch.tensor(np.where(np.isnan(relpos_i), 0, np.where(relpos_i > 999999, 2, 1)))
        relpos_idx_rel_j = torch.tensor(np.where(np.isnan(relpos_j), 0, np.where(relpos_j > 999999, 2, 1)))
        
        relpos_i = torch.tensor(np.nan_to_num(relpos_i, nan=1000000)).float()
        relpos_j = torch.tensor(np.nan_to_num(relpos_j, nan=1000000)).float()
        
        seq_indice_pad = [-1 for j in range(self.max_res_num - len(target_nucleotides))]
        
        seq_indices = torch.tensor([NTS.index(self.pdb_nts[j]) if self.pdb_nts[j] in NTS else 4 for j in target_nucleotides] + seq_indice_pad) + 1
        
        mask = np.concatenate([np.array(self.mask)[0][target_nucleotides], pad], axis = 0)
        mask = np.nan_to_num(mask, False)
        mask = torch.tensor(mask)
        
        target_nucleotides = np.nan_to_num(target_nucleotides + pad, nan=self.S).astype(int).tolist()
        
        rots = self.rots[:, target_nucleotides]
        trans = self.trans[:, target_nucleotides]
            
        relpos_idx_rel_i[mask == 0] = 0
        relpos_idx_rel_j[mask == 0] = 0
        seq_indices[mask == 0] = 0
        
        if(len(rots.size()) == 4):
            rots = rots.squeeze(0)

        if(len(trans.size()) == 3):
            trans = trans.squeeze(0)
        
        return seq_indices.long().to(self.device), rots.float().to(self.device), trans.float().to(self.device), relpos_i.float().to(self.device), relpos_j.float().to(self.device), mask.bool().to(self.device), relpos_idx_rel_i.long().to(self.device), relpos_idx_rel_j.long().to(self.device)

class dataloader_trace(data.Dataset):
    def __init__(self, pairs, trans, str_info, seqidx2segidx, seqidx2pdbidx, dist_map, pdb_nts, mask, adj_con_flag, max_res_num, seq_len, device):
                 #structure_info, max_res_num=200, device = "cpu"):
        super().__init__()
        self.pairs=pairs
        self.trans = trans
        self.str_info=str_info
        self.seqidx2segidx = seqidx2segidx
        self.seqidx2pdbidx = seqidx2pdbidx
        self.dist_map = dist_map
        self.adj_con_flag = adj_con_flag
        self.max_res_num = max_res_num
        self.pdb_nts = pdb_nts
        self.mask = mask
        self.S = seq_len
        self.device = device
        
    def __len__(self):
        return len(self.pairs)
    
    def __getitem__(self, idx):
        base_i = self.pairs[idx][0] 
        base_j = self.pairs[idx][1]
        
        base_nt_i = self.str_info[base_i]["nt_name"]
        base_nt_j = self.str_info[base_j]["nt_name"]
        
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
            
        base_seg_i = self.seqidx2segidx[base_i]
        base_seg_j = self.seqidx2segidx[base_j]
        
        base_pdb_i = self.seqidx2pdbidx[base_i]
        base_pdb_j = self.seqidx2pdbidx[base_j]
        
        target_nucleotides = []
        target_nucleotides.append(base_i)
        target_nucleotides.append(base_j)
        
        dist_ij = np.concatenate([self.dist_map[0, base_i, :].reshape([1,-1]), self.dist_map[0, base_j, :].reshape([1,-1])], axis = 0).T
        dist_ij = np.min(dist_ij, axis = 1) 
        dist_ij[base_i] = np.inf
        dist_ij[base_j] = np.inf 
        
        if(base_i-1 >= 0): 
            if(self.adj_con_flag[base_i]):
                base_pdb = base_i-1
                if(self.str_info[base_i-1]["required_flag"] and base_pdb not in target_nucleotides):
                    target_nucleotides.append(base_pdb)
                    dist_ij[base_pdb] = np.inf
                
                if(base_i-2 >= 0):
                    if(self.adj_con_flag[base_i-1]):
                        base_pdb = base_i-2
                        if(self.str_info[base_i-2]["required_flag"] and base_pdb not in target_nucleotides):
                            target_nucleotides.append(base_pdb)
                            dist_ij[base_i-2] = np.inf
                            
        if(base_i+1 < self.S):
            if(self.adj_con_flag[base_i+1]):
                base_pdb = base_i+1
                if(self.str_info[base_i+1]["required_flag"] and base_pdb not in target_nucleotides):
                    target_nucleotides.append(base_pdb)
                    dist_ij[base_i+1] = np.inf
                
                if(base_i+2 < self.S):
                    if(self.adj_con_flag[base_i+2]):
                        base_pdb = base_i+2
                        if(self.str_info[base_i+2]["required_flag"] and base_pdb not in target_nucleotides):
                            target_nucleotides.append(base_pdb)
                            dist_ij[base_i+2] = np.inf
                        
        if(base_j-1 >= 0): 
            if(self.adj_con_flag[base_j]):
                base_pdb = base_j-1
                if(self.str_info[base_j-1]["required_flag"] and base_pdb not in target_nucleotides):
                    target_nucleotides.append(base_pdb)
                    dist_ij[base_j-1] = np.inf
                
                if(base_j-2 >= 0):
                    if(self.adj_con_flag[base_j-1]):
                        base_pdb = base_j-2
                        if(self.str_info[base_j-2]["required_flag"] and base_pdb not in target_nucleotides):
                            target_nucleotides.append(base_pdb)
                            dist_ij[base_j-2] = np.inf
            
        if(base_j+1 < self.S):
            if(self.adj_con_flag[base_j+1]):
                base_pdb = base_j+1
                if(self.str_info[base_j+1]["required_flag"] and base_pdb not in target_nucleotides):
                    target_nucleotides.append(base_pdb)
                    dist_ij[base_j+1] = np.inf
                
                if(base_j+2 < self.S):
                    if(self.adj_con_flag[base_j+2]):
                        base_pdb = base_j+2
                        if(self.str_info[base_j+2]["required_flag"] and base_pdb not in target_nucleotides):
                            target_nucleotides.append(base_pdb)
                            dist_ij[base_j+2] = np.inf
                                  
        dist_ij_idx = np.argsort(dist_ij)
        dist_ij = dist_ij[dist_ij_idx]
        dist_ij_idx = dist_ij_idx[~np.isinf(dist_ij)]  
        
        spatial_neighbor_num = max(self.max_res_num - len(target_nucleotides), 0)   
        spatial_neighbors = []         
        if(spatial_neighbor_num > 0): spatial_neighbors = dist_ij_idx[:spatial_neighbor_num]
        
        if(len(spatial_neighbors) > 0): target_nucleotides.extend(spatial_neighbors)   
        
        #target_nucleotides_seg = np.array([seqidx2segidx[pdbidx2seqidx[j]] for j in target_nucleotides]) 
        
        pad = [np.nan for j in range(self.max_res_num - len(target_nucleotides))]
        relpos_i = np.array([self.seqidx2pdbidx[j] - base_pdb_i if(self.seqidx2segidx[j] == base_seg_i) else 1000000 for j in target_nucleotides] + pad)  
        relpos_j = np.array([self.seqidx2pdbidx[j] - base_pdb_j if(self.seqidx2segidx[j] == base_seg_j) else 1000000 for j in target_nucleotides] + pad)
        #relpos_i[np.isnan(relpos_i)] = -1000000
        #relpos_j[np.isnan(relpos_j)] = -1000000
        
        relpos_idx_rel_i = torch.tensor(np.where(np.isnan(relpos_i), 0, np.where(relpos_i > 999999, 2, 1)))
        relpos_idx_rel_j = torch.tensor(np.where(np.isnan(relpos_j), 0, np.where(relpos_j > 999999, 2, 1)))
        
        relpos_i = torch.tensor(np.nan_to_num(relpos_i, nan=1000000)).float()
        relpos_j = torch.tensor(np.nan_to_num(relpos_j, nan=1000000)).float()
        
        seq_indice_pad = [-1 for j in range(self.max_res_num - len(target_nucleotides))]
        
        seq_indices = torch.tensor([NTS.index(self.pdb_nts[j]) if self.pdb_nts[j] in NTS else 4 for j in target_nucleotides] + seq_indice_pad) + 1
        
        mask = np.concatenate([np.array(self.mask)[0][target_nucleotides], pad], axis = 0)
        mask = np.nan_to_num(mask, False)
        mask = torch.tensor(mask)
        
        target_nucleotides = np.nan_to_num(target_nucleotides + pad, nan=self.S).astype(int).tolist()
        
        trans = self.trans[:, target_nucleotides]
            
        relpos_idx_rel_i[mask == 0] = 0
        relpos_idx_rel_j[mask == 0] = 0
        seq_indices[mask == 0] = 0

        if(len(trans.size()) == 3):
            trans = trans.squeeze(0)
        
        return seq_indices.long().to(self.device), trans.float().to(self.device), relpos_i.float().to(self.device), relpos_j.float().to(self.device), mask.bool().to(self.device), relpos_idx_rel_i.long().to(self.device), relpos_idx_rel_j.long().to(self.device)

class predict():
    def __init__(self, model):
        self.model = model
        
    def predict(self, loader):
        self.model.eval()
        pred = []
        for i, (seq_indices, rots, trans, seqpos_rel_i, seqpos_rel_j, mask, seqpos_idx_rel_i, seqpos_idx_rel_j) in enumerate(loader):
            
            y_hat = self.model(seq_indices, rots, trans, seqpos_rel_i, seqpos_rel_j, mask, seqpos_idx_rel_i, seqpos_idx_rel_j)
            y_hat = nn.Sigmoid()(y_hat).cpu().detach().clone().squeeze(1).numpy().tolist()
            pred.extend(y_hat)
            
        return pred
    
class predict_trace():
    def __init__(self, model):
        self.model = model
        
    def predict(self, loader):
        self.model.eval()
        pred = []
        for i, (seq_indices, trans, seqpos_rel_i, seqpos_rel_j, mask, seqpos_idx_rel_i, seqpos_idx_rel_j) in enumerate(loader):
            
            y_hat = self.model(seq_indices, trans, seqpos_rel_i, seqpos_rel_j, mask, seqpos_idx_rel_i, seqpos_idx_rel_j)
            y_hat = nn.Sigmoid()(y_hat).cpu().detach().clone().squeeze(1).numpy().tolist()
            pred.extend(y_hat)
            
        return pred

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

#lmcfg={}
#lmcfg['s_in_dim']=5
#lmcfg['z_in_dim']=2
#lmcfg['s_dim']= 512
#lmcfg['z_dim']= 128
#lmcfg['N_elayers']=18
#RNAlm = RNA2nd(lmcfg)

#saved_model =  os.path.join( os.path.dirname(os.path.abspath(__file__)), 'RCLM', 'epoch_67000')
#RNAlm.load_state_dict(torch.load(saved_model,map_location=torch.device('cpu')),strict=False)
#RNAlm.to(device)
#RNAlm.eval()
lmaadic = {'A':0,'G':1,'C':2,'U':3,'a':0,'g':1,'c':2,'u':3,'T':3,'t':3,'-':4}

def parse_seq(inseq):
    seqnpy=np.zeros(len(inseq))
    seq1=np.array(list(inseq))
    seqnpy[seq1=='A']=1
    seqnpy[seq1=='G']=2
    seqnpy[seq1=='C']=3
    seqnpy[seq1=='U']=4
    seqnpy[seq1=='T']=4
    return seqnpy

def get_seq_rep(seq, device):
    aa_type = parse_seq(seq)
    seq_idx = np.arange(len(seq)) + 1

    onehot=torch.from_numpy(aa_type).to(device)
    onehot=F.one_hot(onehot.long(),6).float()

    seq_idx = torch.from_numpy(seq_idx).long().to(device)
    
    return onehot,seq_idx

    """
    alphas = np.array(list(seq))

    seqnpy=np.zeros(len(alphas),dtype=int) + lmaadic['-']
    seq1=np.array(list(alphas))
    keys = list(lmaadic.keys())
    for akey in keys:
        seqnpy[seq1==akey] = lmaadic[akey]
    seqnpy = np.eye(lmaadic['-']+1)[seqnpy]
    fea  = {'aa':seqnpy,'idx':seq_idx,'mask':np.zeros(len(alphas))}
    with torch.no_grad():
        lms,lmz = RNAlm.embedding({
            'aa':torch.FloatTensor(fea['aa']).to(device ),
            'idx':seq_idx,
            'mask':torch.FloatTensor(fea['mask']).to(device ),
            })

    return lms,lmz,onehot,seq_idx
    """
    
class RBFExpand(nn.Module):
    def __init__(self, dmax: float = 30.0, K: int = 64, sigma_scale: float = 0.5, centers_mode: str = "midpoint"):
        super().__init__()
        self.dmax = dmax
        self.K = K

        if centers_mode == "midpoint":
            bin_width = dmax / K
            centers = (torch.arange(K, dtype=torch.float32) + 0.5) * bin_width
            sigma = bin_width * sigma_scale
        elif centers_mode == "edges":
            centers = torch.linspace(0.0, dmax, steps=K, dtype=torch.float32)
            bin_width = dmax / max(1, K - 1)
            sigma = bin_width * sigma_scale
        else:
            raise ValueError("centers_mode must be 'midpoint' or 'edges'")

        self.register_buffer("centers", centers)
        self.register_buffer("sigma", torch.tensor(float(sigma)))

    def forward(self, d: torch.Tensor):
        d = d.clamp(0.0, self.dmax)
        x = d.unsqueeze(-1) - self.centers
        rbf = torch.exp(-0.5 * (x / self.sigma)**2)
        return rbf
    
class dataloader2(data.Dataset):
    def __init__(self, data, maxlength=1200, device = "cpu"):
        super().__init__()
        self.data=data
        self.maxlength=maxlength
        self.rbf_obj = RBFExpand(dmax=30.0, K=32, sigma_scale=1.0)
        self.device = device

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        pred_bp_list, seq, dist_map = self.data[idx]
        
        #dist_map = dist_map.squeeze(0)
        #dist_map[dist_map > 30] = 30.0
        #dist_map = (dist_map.squeeze(0)/30).unsqueeze(-1)
        
        dist_map[dist_map > 30] = 30.0
        dist_map = self.rbf_obj(dist_map)
        dist_map = dist_map.squeeze(0)

        n = len(seq)

        init_bp_map = torch.zeros([n,n])
        if(len(pred_bp_list) > 0):
            init_bp_map[pred_bp_list[0], pred_bp_list[1]] = 1
            init_bp_map[pred_bp_list[1], pred_bp_list[0]] = 1
        init_bp_map = init_bp_map.unsqueeze(-1)

        onehot,seq_idx = get_seq_rep(seq, self.device)
        #lms = lms.clone().detach().cpu()
        #lmz = lmz.clone().detach().cpu()
        onehot = onehot.clone().detach().cpu()

        pfe = get_paring_feature(onehot).unsqueeze(-1)
        outerprot = get_paring_outerprot(onehot)
        
        if(n > self.maxlength):
            init_bp_map = init_bp_map[:self.maxlength, :self.maxlength, :]
            dist_map = dist_map[:self.maxlength, :self.maxlength, :]
            #lmz = lmz[:self.maxlength, :self.maxlength, :]
            #lms = lms[:self.maxlength, :]
            seq_idx = seq_idx[:self.maxlength]
            pfe = pfe[:self.maxlength, :self.maxlength, :]
            outerprot = outerprot[:self.maxlength, :self.maxlength, :]
            n = self.maxlength

        init_bp_map = init_bp_map.to(torch.float32).to(self.device)
        dist_map = dist_map.long().to(self.device)
        #lmz = lmz.to(torch.float32).to(self.device)
        #lms = lms.to(torch.float32).to(self.device)
        seq_idx = seq_idx.to(torch.float32).to(self.device)
        pfe = pfe.to(torch.float32).to(self.device)
        outerprot = outerprot.to(torch.float32).to(self.device)

        return init_bp_map, dist_map, pfe, outerprot, n
    
def round_to_pow2(x, max_len):
    pow2_len = 1 << (x - 1).bit_length()
    return min(pow2_len, max_len)
    
def collate_fn2(batch, maxlength=1024, device="cpu"):
    # batch: list of (init_bp_map, native_bp_map, lmz, pfe, outerprot, weight_bp, seq_len)
    seq_lens = [x[-1] for x in batch]
    max_len = min(max(seq_lens), maxlength)
    max_len = round_to_pow2(max_len, maxlength)

    outs = []

    for init_bp_map, dist_map, pfe, outerprot, n in batch:
        pad_len = max_len - n

        init_bp_map = F.pad(init_bp_map, (0, 0, 0, pad_len, 0, pad_len))
        dist_map = F.pad(dist_map, (0, 0, 0, pad_len, 0, pad_len))
        #lmz = F.pad(lmz, (0, 0, 0, pad_len, 0, pad_len))
        pfe = F.pad(pfe, (0, 0, 0, pad_len, 0, pad_len))
        outerprot = F.pad(outerprot, (0, 0, 0, pad_len, 0, pad_len))

        outs.append((init_bp_map, dist_map, pfe, outerprot))

    init_bp_map, dist_map, pfe, outerprot = zip(*outs)
    return (
        torch.stack(init_bp_map).to(device),
        torch.stack(dist_map).to(device),
        #torch.stack(lmz).to(device),
        torch.stack(pfe).to(device),
        torch.stack(outerprot).to(device),
    )
    
class predict2():
    def __init__(self, model):
        self.model = model
        
    def predict(self, loader, L):
        self.model.eval()
        pred, probs = [], []
        for i, (x, dist_map, pfe, outerprot) in enumerate(loader):
            #y_hat = self.model(lms, lmz.squeeze(0), x.squeeze(0), seq_idx.squeeze(0))
            y_hat = self.model(x, dist_map, pfe, outerprot)
            y_hat = y_hat.squeeze(0)[:L[i], :L[i]]
            
            prob = nn.Sigmoid()(y_hat)
            predicts = (prob >= 0.5).to(torch.long).cpu().detach().clone()
            pred.append(predicts)
            probs.append(prob.cpu().detach().clone().numpy())
            
        return pred, probs
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    