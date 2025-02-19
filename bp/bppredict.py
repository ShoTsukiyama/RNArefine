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
from sklearn.model_selection import train_test_split
#from pytorch3d import transforms
#from torchsummary import summary
import joblib
import pickle
import time
import cpdb
import pandas as pd
import json
from featurizer import get_geopair_fe, get_seq_fe, get_paring_feature, get_paring_outerprot, get_geoseq_fe
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, matthews_corrcoef
from get_ss import main_get_seq_ss, main_get_struct_ss, main_get_init_struct_ss, bp2ss 
from typing import Any, List, Literal, Optional
from Net3 import Net
#from ResNetBP import Net 
#from ResNetMAP import Net 
import subprocess
import cvxpy as cp
import copy
#CSSRPATH = "/RNArefinement/program/src/single/bp_stack_contact/CSSR-Linux64"

RNA_ATOMS = [
    'P',"C4'",'N1','N9'
]
FILL_VALUE = 0

RESIDUE_TO_ATOMS = {
    'A': ['P', 'C4\'', 'N9'],
    'G': ['P', 'C4\'', 'N9'],
    'C': ['P', 'C4\'', 'N1'],
    'U': ['P', 'C4\'', 'N1']
}
ATOM_TO_IDX={"P":0,"C4'":1,"N9":2, "N1":2}
KEEP_ATOM_NUM = 3

def write_txt(filepath, data):
    with open(filepath, "w") as f:
        f.write(data)
        
def load_pickle(filename):
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    return data

def save_pickle(filename, data):
    with open(filename, 'wb') as f:
        pickle.dump(data, f)


def get_ss_info(brk, seqnum):
    pairs = []
    brk1, brk2, brk3, brk4 = 0, 0, 0, 0
    brk1_p, brk2_p, brk3_p, brk4_p = [], [], [], []
    pairs_mat = np.zeros([seqnum, seqnum])
    
    for i in range(len(brk)):
        if(brk[i] == "("):
            brk1 += 1
            brk1_p.append(i)
        
        elif(brk[i] == "["):
            brk2 += 1
            brk2_p.append(i)
        elif(brk[i] == "<"):
            brk3 += 1
            brk3_p.append(i)
        elif(brk[i] == "{"):
            brk4 += 1
            brk4_p.append(i)
        if(brk[i] == ")"):
            brk1 -= 1
            pairs.append((brk1_p[-1], i))
            pairs.append((i, brk1_p[-1]))
            pairs_mat[brk1_p[-1], i] = 1
            pairs_mat[i, brk1_p[-1]] = 1
            brk1_p.pop()
        elif(brk[i] == "]"):
            brk2 -= 1
            pairs.append((brk2_p[-1], i))
            pairs.append((i, brk2_p[-1]))
            pairs_mat[brk2_p[-1], i] = 1
            pairs_mat[i, brk2_p[-1]] = 1
            brk2_p.pop()
        elif(brk[i] == ">"):
            brk3 -= 1
            pairs.append((brk3_p[-1], i))
            pairs.append((i, brk3_p[-1]))
            pairs_mat[brk3_p[-1], i] = 1
            pairs_mat[i, brk3_p[-1]] = 1
            brk3_p.pop()
        elif(brk[i] == "}"):
            brk4 -= 1
            pairs.append((brk4_p[-1], i))
            pairs.append((i, brk4_p[-1]))
            pairs_mat[brk4_p[-1], i] = 1
            pairs_mat[i, brk4_p[-1]] = 1
            brk4_p.pop()  
        
    return pairs, pairs_mat

def filter_dataframe(
    dataframe: pd.DataFrame,
    by_column: str,
    list_of_values: List[Any],
    boolean: bool,
) -> pd.DataFrame:
    """
    Filter function for DataFrame.

    Filters the DataFrame such that the ``by_column`` values have to be
    in the ``list_of_values`` list if ``boolean == True``, or not in the list
    if ``boolean == False``.

    Credit: Arian Jamasb, graphein (https://github.com/a-r-j/graphein)

    :param dataframe: pd.DataFrame to filter.
    :type dataframe: pd.DataFrame
    :param by_column: str denoting column of DataFrame to filter.
    :type by_column: str
    :param list_of_values: List of values to filter with.
    :type list_of_values: List[Any]
    :param boolean: indicates whether to keep or exclude matching
        ``list_of_values``. ``True`` -> in list, ``False`` -> not in list.
    :type boolean: bool
    :returns: Filtered DataFrame.
    :rtype: pd.DataFrame
    """
    df = dataframe.copy()
    df = df[df[by_column].isin(list_of_values) == boolean]
    df.reset_index(inplace=True, drop=True)

    return df

def remove_insertions(
    df: pd.DataFrame, keep: Literal["first", "last"] = "first"
) -> pd.DataFrame:
    """
    This function removes insertions from PDB DataFrames.

    Credit: Arian Jamasb, graphein (https://github.com/a-r-j/graphein)

    :param df: RNA Structure dataframe to remove insertions from.
    :type df: pd.DataFrame
    :param keep: Specifies which insertion to keep. Options are ``"first"`` or
        ``"last"``. Default is ``"first"``.
    :type keep: Literal["first", "last"]
    :return: RNA structure dataframe with insertions removed
    :rtype: pd.DataFrame
    """
    # Catches unnamed insertions
    duplicates = df.duplicated(
        subset=["chain_id", "residue_number", "atom_name", "alt_loc"],
        keep=keep,
    )
    df = df[~duplicates]

    return filter_dataframe(
        df, by_column="insertion", list_of_values=[""], boolean=True
    )

def filter_atoms(row):
    residue = row['residue_name']
    atom_name = row['atom_name']
    
    # 対応するアトムリストを取得
    valid_atoms = RESIDUE_TO_ATOMS.get(residue, [])
    
    # アトムが対応するリストに含まれているかを確認
    return atom_name in valid_atoms

def df_to_tensor(
    df: pd.DataFrame,
    fill_value: float = FILL_VALUE,
    center: bool = True
):
    """
    Transforms a DataFrame of an RNA structure into a
    ``length x num_atoms x 3`` tensor.

    Credit: Arian Jamasb, graphein (https://github.com/a-r-j/graphein)

    :param df: DataFrame of protein structure.
    :type df: pd.DataFrame
    :param atoms_to_keep: List of atom types to retain in the tensor.
    :type atoms_to_keep: List[str]
    :param fill_value: Value to fill missing entries with. Defaults to ``1e-5``.
    :type fill_value: float
    :param center: Whether to center the structure at the origin. Defaults to ``True``.
    :type center: bool
    :returns: ``Length x Num_Atoms (default 37) x 3`` tensor.
    :rtype: torch.Tensor
    """
    if center:
        df.x_coord -= df.x_coord.mean()
        df.y_coord -= df.y_coord.mean()
        df.z_coord -= df.z_coord.mean()

    seq = "".join(df.drop_duplicates(subset = ["residue_number", "residue_name"])["residue_name"].values.tolist())

    #num_residues = len([res for res in df.residue_number.unique()])
    num_residues = len(seq)
    #df = df.loc[df["atom_name"].isin(atoms_to_keep)]
    df = df[df.apply(filter_atoms, axis=1)]
    residue_indices = pd.factorize(np.array(df.residue_number))[0]
    #atom_indices = df["atom_name"].map(lambda x: atoms_to_keep.index(x)).values
    atom_indices = df["atom_name"].map(lambda x: ATOM_TO_IDX[x]).values

    positions = (
        torch.zeros((num_residues, KEEP_ATOM_NUM, 3)) + fill_value
    ).float()
    positions[residue_indices, atom_indices] = torch.tensor(
        df[["x_coord", "y_coord", "z_coord"]].values
    ).float()
    
    positions_idx = (
        torch.zeros((num_residues, KEEP_ATOM_NUM, 1))
    ).int()
    positions_idx[residue_indices, atom_indices] = 1
    
    return positions, positions_idx, seq

def get_info_from_pdb(pdb_path):
    df = cpdb.parse(pdb_path, df=True)
    df = remove_insertions(df)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
   
    coords, coords_mask, seq = df_to_tensor(df, center=False)
    
    return coords, coords_mask, seq


def get_feature(seq, coords, seq_bp_maps, init_bp_maps, maxlength, device):
    n = len(seq)

    R_matrices, T_vectors = [], []
    for i in range(len(coords)):
        R_matrice, T_vector = get_geopair_fe(coords[i])
        R_matrices.append(R_matrice)
        T_vectors.append(T_vector)

    R_matrices=torch.stack(R_matrices)
    T_vectors=torch.stack(T_vectors)

    seq_onehot = get_seq_fe(seq)
    pair_fe = []
    pair_fe.append(get_paring_feature(seq_onehot).unsqueeze(-1))
    pair_fe.append(get_paring_outerprot(seq_onehot))

    pair_fe.append(seq_bp_maps.permute(1,2,0))
    pair_fe.append(init_bp_maps.permute(1,2,0))
    pair_fe = torch.concat(pair_fe, dim = -1)
    
    if(n > maxlength):
        seq_onehot = seq_onehot[:maxlength, :]
        pair_fe = pair_fe[:maxlength, :maxlength, :]
        
        stack_map = stack_map[:maxlength, :maxlength]
        bp_map = bp_map[:maxlength, :maxlength]

        R_matrices = R_matrices[:,:maxlength, :maxlength]
        T_vectors = T_vectors[:,:maxlength, :maxlength]
        
    R_matrices = R_matrices.to(torch.float32).to(device).unsqueeze(0)
    T_vectors = T_vectors.to(torch.float32).to(device).unsqueeze(0)

    seq_onehot = seq_onehot.to(torch.float32).to(device).unsqueeze(0)
    pair_fe = pair_fe.to(torch.float32).to(device).unsqueeze(0)

    return R_matrices, T_vectors, seq_onehot, pair_fe
    
def loss_fn(y_pred, y_true, weight):
    #loss = torch.nn.BCEWithLogitsLoss(reduction='none')
    loss = torch.nn.BCELoss(reduction='none')
    return (loss(y_pred, y_true) * weight).sum()

def calc_measures(pred_, true_):
    pred = np.zeros_like(pred_)
    pred[pred_ > 0.5] = 1
    
    prc = precision_score(true_, pred)
    rcal = recall_score(true_, pred)
    f1 = f1_score(true_, pred)
    mcc = matthews_corrcoef(true_, pred)

    return pd.DataFrame([prc, rcal, f1, mcc], index = ["Precision", "Recall", "F1 score", "MCC"])

class predict():
    def __init__(self, modelpath, maxlength=800, device = "cuda:1"):
        self.model = Net(dim_seq_fe=4, dim_hidden = 3*16).to(device)
        self.model.load_state_dict(torch.load(modelpath, map_location=torch.device(device)))
        
        self.maxlength = maxlength
        self.device = device
        
    def __call__(self, coords, seq_ss_fes, str_ss_fes, seq, outpath):
        
        R_matrices, T_vectors, seq_onehot, pair_fe = get_feature(seq, coords, seq_ss_fes, str_ss_fes, self.maxlength, self.device)
            
        self.model.eval()
        out_raw = self.model(seq_onehot, pair_fe, R_matrices, T_vectors)
        #out = nn.Tanh()(out_raw)
        out = nn.Sigmoid()(out_raw)
        out = out.cpu().detach().clone().squeeze(0).numpy()
        out_raw = out_raw.cpu().detach().clone().squeeze(0).numpy()

        #out_ss = bp2ss(out_raw)
        return out 
        #return ((out-np.min(out))/(np.max(out)-np.min(out)))-0.5
        
        """
        out = nn.Sigmoid()(out).numpy()

        ss_fes = torch.concat([torch.tensor(seq_ss_fes), torch.tensor(str_ss_fes)], dim=0)
        ss_fes = torch.mean(ss_fes,dim = 0).detach().clone().numpy()

        out = out+ss_fes
        pairs = np.array(np.where(out >= 0.25))

        if(len(pairs) > 0):
            scores = out[pairs[0], pairs[1]][None,:]
            pairs = np.concatenate([pairs, scores])
            pairs = pairs.T.tolist()
            pairs = pd.DataFrame(pairs, columns = ["base_i", "base_j","scores"])
            pairs = pairs[pairs["base_i"] < pairs["base_j"]]
            pairs['base_i'] = pairs['base_i'].astype(int)
            pairs['base_j'] = pairs['base_j'].astype(int)
        else:
            pairs = pd.DataFrame([], columns = ["base_i", "base_j","scores"])

        pairs.to_csv(outpath, index=False)
        """

def get_init_ss(pdbfile, seq_len, CSSRPATH):
    command = [CSSRPATH, pdbfile, "-", "--outfmt", "1"]
    result = subprocess.run(command, capture_output=True, text=True)
    
    pairs, _ = get_ss_info(result.stdout, seq_len)
    pairs = np.array(pairs)

    return pairs

def out_ss(filepath, ss_info, score_mat):
    pairs = np.array(np.where(ss_info >= 0.5))
    
    if(len(pairs) > 0):
        scores = score_mat[pairs[0], pairs[1]][None,:]
        pairs = np.concatenate([pairs, scores])
        pairs = pairs.T.tolist()
        pairs = pd.DataFrame(pairs, columns = ["base_i", "base_j","scores"])
        pairs = pairs[pairs["base_i"] < pairs["base_j"]]
        pairs['base_i'] = pairs['base_i'].astype(int)
        pairs['base_j'] = pairs['base_j'].astype(int)
    else:
        pairs = pd.DataFrame([], columns = ["base_i", "base_j","scores"])

    pairs.to_csv(filepath, index=False)
    
def open_json(filepath):
    # JSONファイルを開いて読み込む
    with open(filepath, 'r', encoding='utf-8') as f:
        data = json.load(f)
        
    return data

#if __name__ == '__main__':
def bppredictmainpdb(initpdbpath, templatepdbpath, outpath, weight):
    #initpdbpath = sys.argv[1]
    #templatepdbpath = [sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]]
    
    #outpath = sys.argv[8]

    workingpath = os.path.dirname(os.path.abspath(__file__))
    modelpath = "{}/model.pt".format(workingpath)
    
    #_, _, init_seq = get_info_from_pdb(initpdbpath)
    #init_bp = main_get_struct_ss(workingpath, init_seq, initpdbpath)
    
    coords, str_ss_fes = [], []
    for i in range(len(templatepdbpath)):
        assert os.path.exists(templatepdbpath[i]), "Target pdb file doesn't exist."
        coord, _, seq = get_info_from_pdb(templatepdbpath[i])
        coords.append(coord)
  
        str_ss_fe = main_get_struct_ss(workingpath, seq, templatepdbpath[i])
        str_ss_fes.extend(str_ss_fe)

    str_ss_fes = torch.stack(str_ss_fes)
    seq_ss_fes = main_get_seq_ss(workingpath, seq)

    seq_ss_fes= torch.stack(seq_ss_fes)

    cls = predict(modelpath)
    predict_bp = cls(coords, seq_ss_fes, str_ss_fes, seq, outpath)
    
    _, _, init_seq = get_info_from_pdb(initpdbpath)
    init_bp = main_get_init_struct_ss(workingpath, init_seq, initpdbpath)
    
    weight_init = 1-weight
    weight_init_each = weight_init/6.0
  
    bp_mat = init_bp[0].numpy() * (weight_init_each*3) + init_bp[1].numpy() * weight_init_each + init_bp[2].numpy() * weight_init_each + init_bp[3].numpy() * weight_init_each + predict_bp * weight
    
    #bp_mat = init_bp[0].numpy() * 0.3 + init_bp[1].numpy() * 0.1 + init_bp[2].numpy() * 0.1 + init_bp[3].numpy() * 0.1 + predict_bp* 0.4
    bp_mat[bp_mat < 0.5] = -1

    final_ss_mat = bp2ss(bp_mat)
    
    outdir = os.path.dirname(outpath)
    outfile = os.path.basename(outpath)
    
    os.makedirs(f"{outdir}", exist_ok = True)
    out_ss(f"{outpath}", final_ss_mat, bp_mat)
    
def bppredictmainfasta(templatepdbpath, outpath):
    #initpdbpath = sys.argv[1]
    #templatepdbpath = [sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]]
    
    #outpath = sys.argv[8]

    workingpath = os.path.dirname(os.path.abspath(__file__))
    modelpath = "{}/model.pt".format(workingpath)
    
    #_, _, init_seq = get_info_from_pdb(initpdbpath)
    #init_bp = main_get_struct_ss(workingpath, init_seq, initpdbpath)
    
    coords, str_ss_fes = [], []
    for i in range(len(templatepdbpath)):
        assert os.path.exists(templatepdbpath[i]), "Target pdb file doesn't exist."
        coord, _, seq = get_info_from_pdb(templatepdbpath[i])
        coords.append(coord)
  
        str_ss_fe = main_get_struct_ss(workingpath, seq, templatepdbpath[i])
        str_ss_fes.extend(str_ss_fe)

    str_ss_fes = torch.stack(str_ss_fes)
    seq_ss_fes = main_get_seq_ss(workingpath, seq)

    seq_ss_fes= torch.stack(seq_ss_fes)

    cls = predict(modelpath)
    predict_bp = cls(coords, seq_ss_fes, str_ss_fes, seq, outpath)
    
    predict_bp[predict_bp < 0.5] = -1
    final_ss_mat = bp2ss(predict_bp)
    
    outdir = os.path.dirname(outpath)
    outfile = os.path.basename(outpath)
    
    os.makedirs(f"{outdir}", exist_ok = True)
    out_ss(f"{outpath}", final_ss_mat, predict_bp)
    
    
    
    









