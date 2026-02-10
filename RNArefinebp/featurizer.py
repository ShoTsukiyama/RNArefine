#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 12:46:23 2024

@author: tsukiyamashou
"""
import os
import sys
import random
import numpy as np
from sklearn.model_selection import train_test_split
import torch
from torch import nn
from torch.nn import functional as F
import torch.utils.data as data
from torch.utils.data import DataLoader
from torch.autograd import Variable
from torch.utils.data import BatchSampler
#from torchsummary import summary
import joblib
import pickle
import time
import copy
#from ResNetBP import Net 
#from ResNetMAP import Net 

DISTANCE_EPS = 0.001

def discritization(tens, num_bins, min_val_, max_val_):
    min_val = min_val_+(max_val_-min_val_)/num_bins
    max_val = max_val_-(max_val_-min_val_)/num_bins

    bin_edges = torch.linspace(min_val, max_val, num_bins-1)
    indices = torch.bucketize(tens,bin_edges, right=False)
    one_hot = F.one_hot(indices, num_classes=num_bins)
    return one_hot

def get_paring_outerprot(seq_onehot):
    seq_map = []
    for i in range(4):
        for j in range(i, 4):
            outer = (torch.outer(seq_onehot[:, i], seq_onehot[:, j]))[:,:,None] + (torch.outer(seq_onehot[:, j], seq_onehot[:, i]))[:,:,None]
            
            if(i == j): outer = outer/2
            
            seq_map.append(outer)
    seq_map = torch.cat(seq_map, dim = -1)

    return seq_map

def get_paring_feature(seq_onehot):
    gc_map = seq_onehot[None, :, 1] * seq_onehot[:, None, 2]
    au_map = seq_onehot[None, :, 0] * seq_onehot[:, None, 3]
    gu_map = seq_onehot[None, :, 1] * seq_onehot[:, None, 3]
    cg_map = seq_onehot[None, :, 2] * seq_onehot[:, None, 1]
    ua_map = seq_onehot[None, :, 3] * seq_onehot[:, None, 0]
    ug_map = seq_onehot[None, :, 3] * seq_onehot[:, None, 1]
    init_pair_map = (3*gc_map + 2*au_map + 0.8*gu_map + 3*cg_map + 2*ua_map + 0.8*ug_map).float()
    init_pair_map = init_pair_map * torch.tril(torch.ones(init_pair_map.size(0), init_pair_map.size(1)), diagonal=-1)

    L = init_pair_map.size(0)
    pair_map = copy.deepcopy(init_pair_map)
    pair_map_shift = copy.deepcopy(init_pair_map)
    pair_idx = torch.zeros_like(init_pair_map)
    pair_idx[pair_map > 0] = 1
    
    i=1
    while True:
        pair_map_shift = pair_map_shift[1:,:L-1]
        pair_map_shift=F.pad(pair_map_shift, (1, 0, 0, 1))
        pair_idx = pair_idx*(pair_map_shift > 0)
        
        pair_map = pair_map + (np.exp(-(i**2))) * pair_idx * pair_map_shift
        if(pair_idx.sum() == 0):
            break
        
        i+=1
        
    pair_map_shift = copy.deepcopy(init_pair_map)
    pair_idx = torch.zeros_like(init_pair_map)
    pair_idx[init_pair_map > 0] = 1
       
    i=1
    while True:
        pair_map_shift = pair_map_shift[:L-1,1:]
        pair_map_shift = F.pad(pair_map_shift, (0, 1, 1, 0))
        pair_idx = pair_idx*(pair_map_shift > 0)
        pair_map = pair_map + (np.exp(-(i**2))) * pair_idx * pair_map_shift
        
        if(pair_idx.sum() == 0):
            break
        
        i+=1
    
    pair_map = pair_map + pair_map.T
    return pair_map

    
    
