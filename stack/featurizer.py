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
from pytorch3d import transforms
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
        for j in range(4):
            outer = (torch.outer(seq_onehot[:, i], seq_onehot[:, j]))[:,:,None]
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


def normed_vec(V: torch.Tensor, distance_eps: float = DISTANCE_EPS) -> torch.Tensor:
    """Normalized vectors with distance smoothing.

    This normalization is computed as `U = V / sqrt(|V|^2 + eps)` to avoid cusps
    and gradient discontinuities.

    Args:
        V (Tensor): Batch of vectors with shape `(..., num_dims)`.
        distance_eps (float, optional): Distance smoothing parameter for
            for computing distances as `sqrt(sum_sq) -> sqrt(sum_sq + eps)`.
            Default: 1E-3.

    Returns:
        U (Tensor): Batch of normalized vectors with shape `(..., num_dims)`.
    """
    # Unit vector from i to j
    mag_sq = (V ** 2).sum(dim=-1, keepdim=True)
    mag = torch.sqrt(mag_sq + distance_eps)
    U = V / mag
    return U

def normed_cross(
    V1: torch.Tensor, V2: torch.Tensor, distance_eps: float = DISTANCE_EPS
) -> torch.Tensor:
    """Normalized cross product between vectors.

    This normalization is computed as `U = V / sqrt(|V|^2 + eps)` to avoid cusps
    and gradient discontinuities.

    Args:
        V1 (Tensor): Batch of vectors with shape `(..., 3)`.
        V2 (Tensor): Batch of vectors with shape `(..., 3)`.
        distance_eps (float, optional): Distance smoothing parameter for
            for computing distances as `sqrt(sum_sq) -> sqrt(sum_sq + eps)`.
            Default: 1E-3.

    Returns:
        C (Tensor): Batch of cross products `v_1 x v_2` with shape `(..., 3)`.
    """
    C = normed_vec(torch.cross(V1, V2, dim=-1), distance_eps=distance_eps)
    return C

def lengths(
    atom_i: torch.Tensor, atom_j: torch.Tensor, distance_eps: float = DISTANCE_EPS
) -> torch.Tensor:
    """Batched bond lengths given batches of atom i and j.

    Args:
        atom_i (Tensor): Atom `i` coordinates with shape `(..., 3)`.
        atom_j (Tensor): Atom `j` coordinates with shape `(..., 3)`.
        distance_eps (float, optional): Distance smoothing parameter for
            for computing distances as `sqrt(sum_sq) -> sqrt(sum_sq + eps)`.
            Default: 1E-3.

    Returns:
        L (Tensor): Elementwise bond lengths `||x_i - x_j||` with shape `(...)`.
    """
    # Bond length of i-j
    dX = atom_j - atom_i
    L = torch.sqrt((dX ** 2).sum(dim=-1) + distance_eps)
    return L

def angles(
    atom_i: torch.Tensor,
    atom_j: torch.Tensor,
    atom_k: torch.Tensor,
    distance_eps: float = DISTANCE_EPS,
    degrees: bool = False,
) -> torch.Tensor:
    """Batched bond angles given atoms `i-j-k`.

    Args:
        atom_i (Tensor): Atom `i` coordinates with shape `(..., 3)`.
        atom_j (Tensor): Atom `j` coordinates with shape `(..., 3)`.
        atom_k (Tensor): Atom `k` coordinates with shape `(..., 3)`.
        distance_eps (float, optional): Distance smoothing parameter for
            for computing distances as `sqrt(sum_sq) -> sqrt(sum_sq + eps)`.
            Default: 1E-3.
        degrees (bool, optional): If True, convert to degrees. Default: False.

    Returns:
        A (Tensor): Elementwise bond angles with shape `(...)`.
    """
    # Bond angle of i-j-k
    U_ji = normed_vec(atom_i - atom_j, distance_eps=distance_eps)
    U_jk = normed_vec(atom_k - atom_j, distance_eps=distance_eps)
    inner_prod = torch.einsum("ix,ix->i", U_ji, U_jk)
    inner_prod = torch.clamp(inner_prod, -1, 1)
    A = torch.acos(inner_prod)
    if degrees:
        A = A * 180.0 / np.pi
    return A

def dihedrals(
    atom_i: torch.Tensor,
    atom_j: torch.Tensor,
    atom_k: torch.Tensor,
    atom_l: torch.Tensor,
    distance_eps: float = DISTANCE_EPS,
    degrees: bool = False,
) -> torch.Tensor:
    """Batched bond dihedrals given atoms `i-j-k-l`.

    Args:
        atom_i (Tensor): Atom `i` coordinates with shape `(..., 3)`.
        atom_j (Tensor): Atom `j` coordinates with shape `(..., 3)`.
        atom_k (Tensor): Atom `k` coordinates with shape `(..., 3)`.
        atom_l (Tensor): Atom `l` coordinates with shape `(..., 3)`.
        distance_eps (float, optional): Distance smoothing parameter for
            for computing distances as `sqrt(sum_sq) -> sqrt(sum_sq + eps)`.
            Default: 1E-3.
        degrees (bool, optional): If True, convert to degrees. Default: False.

    Returns:
        D (Tensor): Elementwise bond dihedrals with shape `(...)`.
    """
    U_ij = normed_vec(atom_j - atom_i, distance_eps=distance_eps)
    U_jk = normed_vec(atom_k - atom_j, distance_eps=distance_eps)
    U_kl = normed_vec(atom_l - atom_k, distance_eps=distance_eps)
    normal_ijk = normed_cross(U_ij, U_jk, distance_eps=distance_eps)
    normal_jkl = normed_cross(U_jk, U_kl, distance_eps=distance_eps)
    # _inner_product = lambda a, b: torch.einsum("bix,bix->bi", a, b)
    _inner_product = lambda a, b: (a * b).sum(-1)
    cos_dihedrals = _inner_product(normal_ijk, normal_jkl)
    angle_sign = _inner_product(U_ij, normal_jkl)
    cos_dihedrals = torch.clamp(cos_dihedrals, -1, 1)
    D = torch.sign(angle_sign) * torch.acos(cos_dihedrals)
    if degrees:
        D = D * 180.0 / np.pi
    return D

def get_local_system(atom_cent, atom_r, atom_l):
    vec1 = atom_r - atom_cent
    vec2 = atom_l - atom_cent
    vec_z = vec1 + vec2
    
    vec_z = normed_vec(vec_z)
    vec_y = normed_cross(vec1, vec2)
    vec_x = normed_cross(vec_y, vec_z)
    
    return torch.concat([vec_x.unsqueeze(-2), vec_y.unsqueeze(-2), vec_z.unsqueeze(-2)], dim = -2)
   
def get_Rmatrix(lsystem):
    lsystem_obj=torch.permute(lsystem,  (0, 2, 1))
    invlsystemA = torch.linalg.inv(lsystem_obj)
    lsystemB = lsystem_obj
    R_matrices = torch.einsum('mik,nkj->nmji', lsystemB, invlsystemA)
    return R_matrices

def get_Tvector(base_i, base_j, R_matrices):
    T_vector = base_j[None,:,:].repeat(base_j.size(0),1, 1) - torch.einsum('mnij,mni->mnj', R_matrices, base_i[:,None,:].repeat(1, base_i.size(0), 1))
    return T_vector



def get_quaternion(lsystem):
    R_matrices = get_Rmatrix(lsystem)
    #lsystem_obj=torch.permute(lsystem,  (0, 2, 1))
    #invlsystemA = torch.linalg.inv(lsystem_obj)
    #lsystemB = lsystem_obj
    #R_matrices = torch.einsum('mik,nkj->nmji', lsystemB, invlsystemA)
    
    quaternions = transforms.matrix_to_quaternion(R_matrices)
    quaternions = torch.permute(quaternions,  (1,0,2))
    
    return quaternions

def get_geopair_fe(coords):
    lsystem = get_local_system(coords[:,1,:], coords[:,0,:], coords[:,2,:])
    R_matrices = get_Rmatrix(lsystem)
    T_vector = get_Tvector(lsystem[:, 0, :], lsystem[:, 0, :], R_matrices)

    #R_matrices = R_matrices.unsqueeze(2).expand(-1, -1, times, -1, -1)
    #T_vector = T_vector.unsqueeze(2).expand(-1, -1, times, -1)
    
    #print(R_matrices.size())
    #print(T_vector.size())

    return R_matrices, T_vector
    #print(lsystem[:,None,2,:].repeat(1,lsystem.size(0), 1).size())
    #lsystem_tmp = lsystem[:,None,2,:].repeat(1,lsystem.size(0), 1)
    #data = torch.einsum('ijmn,ijm->ijn',R_matrices,lsystem_tmp+T_vector)

    sys.exit()
    quaternions = get_quaternion(lsystem)
    
    #plengthsP = lengths(coords[:, None, 0, :].repeat(1, coords.size(0), 1), coords[None, :, 0, :].repeat(coords.size(0), 1, 1))
    plengthsC4 = lengths(coords[:, None, 1, :].repeat(1, coords.size(0), 1), coords[None, :, 1, :].repeat(coords.size(0), 1, 1))
    #plengthsN = lengths(coords[:, None, 2, :].repeat(1, coords.size(0), 1), coords[None, :, 2, :].repeat(coords.size(0), 1, 1))
    
    center_coords = coords.mean(1)
    plengthscent = lengths(center_coords[:, None, :].repeat(1, coords.size(0), 1), center_coords[None, :, :].repeat(coords.size(0), 1, 1))
    
    #plengthsC4[plengthsC4 > 20] = 20
    #plengthscent[plengthscent > 20] = 20
    #plengthsC4 = plengthsC4/20
    #plengthscent = plengthscent/20
    
    plengths = torch.concat([plengthsC4.unsqueeze(-1), plengthscent.unsqueeze(-1)], dim = -1)
    #plengths[plengths > 20] = 20
    #print(torch.max(plengths)) 
    #plengths=torch.log(plengths + 1e-5)
  
    quaternions_norm = torch.norm(quaternions, p=2, dim=-1, keepdim=True)
    quaternions/=quaternions_norm

    quaternions_tens = discritization(quaternions, 50, -1, 1)
    quaternions_tens = quaternions_tens.reshape(quaternions_tens.size(0), quaternions_tens.size(1), quaternions_tens.size(2)*quaternions_tens.size(3))

    plengths_tens = discritization(plengths, 50, 0, 20)
    plengths_tens = plengths_tens.reshape(plengths_tens.size(0), plengths_tens.size(1), plengths_tens.size(2)*plengths_tens.size(3))

    return quaternions, torch.permute(quaternions,  (1,0,2)), plengths_tens

def _pad_pack(D, A, L, N_D, N_A, N_L, num_residues):
    # Pad and pack together the components
    D = F.pad(D, (1, 2))
    A = F.pad(A, (1, 1))
    L = F.pad(L, (0, 1))
    N_D = F.pad(N_D, (0, 1))
    D, A, L = [x.reshape(num_residues, 2) for x in [D, A, L]]
    _pack = lambda a, b: torch.cat([a, b.unsqueeze(-1)], dim=-1)
    L = _pack(L, N_L)
    A = _pack(A, N_A)
    D = _pack(D, N_D)
    return D, A, L

def get_geoseq_fe(coords):
    _lengths = lambda Xi, Xj: lengths(Xi, Xj, distance_eps=1e-5)
    _angles = lambda Xi, Xj, Xk: angles(Xi, Xj, Xk, distance_eps=1e-5)
    _dihedrals = lambda Xi, Xj, Xk, Xl: dihedrals(Xi, Xj, Xk, Xl, distance_eps=1e-5)
    
    bbcoords = coords[:, :2, :]
    num_residues, _, _ = bbcoords.shape
    bbcoords = bbcoords.reshape(2 * num_residues, 3)
    
    PC4p_L = _lengths(bbcoords[1:, :], bbcoords[:-1, :])
    PC4p_A = _angles(bbcoords[:-2, :], bbcoords[1:-1, :], bbcoords[2:, :])
    PC4p_D = _dihedrals(bbcoords[:-3, :], bbcoords[1:-2, :], bbcoords[2:-1, :], bbcoords[3:, :])
    
    X_P, X_C4p, X_N = coords.unbind(dim=1)
    X_P_next = coords[1:, 0, :]
    
    N_L = _lengths(X_C4p, X_N)
    N_A = _angles(X_P, X_C4p, X_N)
    N_D = _dihedrals(X_P_next, X_N[:-1, :], X_C4p[:-1, :], X_P[:-1, :])
    
    D, A, L = _pad_pack(PC4p_D, PC4p_A, PC4p_L, N_D, N_A, N_L, num_residues)
    
    angle_stack = torch.cat([D, A], dim=-1)
    L = torch.log(L + 1e-5)
    internal_coords_feat = torch.cat([torch.cos(angle_stack), torch.sin(angle_stack), L], dim=-1)
    
    return internal_coords_feat

def parse_seq(inseq):
    seqnpy=np.zeros(len(inseq))
    seq1=np.array(list(inseq))
    seqnpy[seq1=='A']=0
    seqnpy[seq1=='G']=1
    seqnpy[seq1=='C']=2
    seqnpy[seq1=='U']=3
    seqnpy[seq1=='T']=3
    return seqnpy

def get_seq_fe(seq):
    seq_onehot = torch.tensor(parse_seq(seq))
    seq_onehot = F.one_hot(seq_onehot.long(),4).float()
    
    return seq_onehot

    
    
