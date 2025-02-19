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
from featurizer import get_geopair_fe, get_seq_fe, get_paring_feature, get_paring_outerprot, get_geoseq_fe
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, matthews_corrcoef
#from get_ss import main_get_seq_ss, main_get_struct_ss  
#from get_stack import get_bp_clarna_stack
from typing import Any, List, Literal, Optional
from Net4 import Net
#from ResNetBP import Net 
#from ResNetMAP import Net 
import subprocess
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


def get_feature(seq, coords, bp_list, maxlength, device):
    #coords = coords.unsqueeze(1)
    dist_map = torch.sqrt(((coords[:,None,:] - coords[None,:,:])**2).sum(-1))
    n = len(seq)
    seq_onehot = get_seq_fe(seq)
    pair_fe = []
    pair_fe.append(dist_map)
    pair_fe.append(get_paring_feature(seq_onehot).unsqueeze(-1))
    pair_fe.append(get_paring_outerprot(seq_onehot))
    
    bp_map = torch.zeros([n,n])
    if(len(bp_list) != 0):
        bp_map[bp_list[0], bp_list[1]] = 1
        bp_map[bp_list[1], bp_list[0]] = 1

    pair_fe.append(bp_map.unsqueeze(-1))
    pair_fe = torch.concat(pair_fe, dim = -1)
    if(n > maxlength):
        seq_onehot = seq_onehot[:maxlength, :]
        pair_fe = pair_fe[:maxlength, :maxlength, :]
        bp_map = bp_map[:maxlength, :maxlength]
        
    seq_onehot = seq_onehot.to(torch.float32).to(device).unsqueeze(0)
    pair_fe = pair_fe.to(torch.float32).to(device).unsqueeze(0)

    pair_fe = torch.permute(pair_fe, (0,3,1,2))
    
    return pair_fe
    
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
   
def loss_bcewithlogits(y_pred, y_true, weight):
    criterion = torch.nn.BCEWithLogitsLoss(reduction='none')
    #criterion = torch.nn.BCELoss(reduction='none')
    return (criterion(y_pred, y_true) * weight).sum()

def loss_ce(y_pred, y_true, weight):
    criterion = nn.CrossEntropyLoss(reduction='none')
    return (criterion(y_pred, y_true) * weight).sum()

class evaluate():
    def __init__(self, modelpath, maxlength=800, device = "cuda:1"):
        self.model = Net(dim_hidden = 3*16).to(device)
        self.model.load_state_dict(torch.load(modelpath, map_location=torch.device(device)))
        
        self.maxlength = maxlength
        self.device = device
        
    def __call__(self, coords, seq, bp_list):
        
        pair_fe = get_feature(seq, coords, bp_list, self.maxlength, self.device)
        
        #n = stack_map.size(0)    
        
        self.model.eval()
        #out_bp, out_stack, out_contact = self.model(seq_onehot, pair_fe, R_matrices, T_vectors)
        #out_bp, out_stack = self.model(seq_onehot, pair_fe, R_matrices, T_vectors)
        #out_bp = self.model(seq_onehot, pair_fe, R_matrices, T_vectors)
        out = self.model(pair_fe)
        
        out = nn.Sigmoid()(out).squeeze(0).cpu().detach().clone().numpy()
        pairs = np.array(np.where(out >= 0.5))
        
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
            
        return pairs, out
        

#def get_init_ss(pdbfile, seq_len):
    #command = [CSSRPATH, pdbfile, "-", "--outfmt", "1"]
    #result = subprocess.run(command, capture_output=True, text=True)
    
    #pairs, _ = get_ss_info(result.stdout, seq_len)
    #pairs = np.array(pairs)

    #return pairs


#if __name__ == '__main__':
    #pdbpath = sys.argv[1]
    #bppath = sys.argv[2]
    #outpath = sys.argv[3]
    
def stacking_prediction_main(pdbpath, bppath, outpath):

    workingpath = os.path.dirname(os.path.abspath(__file__))
    modelpath = "{}/model.pt".format(workingpath)

    coord, _, seq = get_info_from_pdb(pdbpath)
    bp_list = pd.read_csv(bppath)
    bp_list = (np.array(bp_list.loc[:, ["base_i", "base_j"]].values.tolist())).T

    #init_stack_pairs = get_bp_clarna_stack(pdbpath, seq)
    #init_stack_pairs["scores"]=1
    #init_stack_pairs["base_i"] = init_stack_pairs["base_i"] - 1
    #init_stack_pairs["base_j"] = init_stack_pairs["base_j"] - 1

    eval_cls = evaluate(modelpath)
    pairs, out = eval_cls(coord[:,0:2,:], seq, bp_list)

    #pairs = pd.concat([pred_pairs, init_stack_pairs]).drop_duplicates(subset=["base_i", "base_j"])

    pairs.to_csv(f"{outpath}/stack.list", index=False)
    #save_pickle(f"{outpath}/stack_predicted_matrix.pkl", out)








