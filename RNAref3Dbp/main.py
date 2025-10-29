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
import re
from pdbparser import PDBLoader
from affine import build_affine3d_from_coordinates
from torch.utils.data import DataLoader
from bppredictor import predict_bp_from_backbonestr, predict_bp_from_initbp
from stackpredictor import get_bpmap, predict_adjstack_from_str, predict_nonadjstack_from_str
import pickle

#sys.path.append(os.path.join(os.path.dirname(__file__), './bp'))

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

bp_pattern1 = re.compile(r'^([^,]+),([^,]+),(\d+),(\d+),(\d+\.\d+)$')
bp_pattern2 = re.compile(r'^([^,]+),([^,]+),(\d+),(\d+)')

def get_seq(pdbpath):
    # PDBパーサーの作成
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("RNA", pdbpath)
    
    # RNA配列を取得
    rna_sequences = []
    for model in structure:
        for chain in model:
            seq = []
            for residue in chain:
                res_name = residue.get_resname().strip()
                if res_name in rna_codes:
                    seq.append(rna_codes[res_name])
            if seq:
                seq_record = SeqRecord(Seq("".join(seq)), id=f"Chain_{chain.id}", description="")
                rna_sequences.append(seq_record)
                
    # FASTA形式で出力
    #with open("output.fasta", "w") as fasta_file:
        #write(rna_sequences, fasta_file, "fasta")
                
    return rna_sequences

def output_fasta(outpath, seq_id, sequences):
    with open(outpath, "w") as f:
        f.write("\n".join([f">{seq_id}", f"{sequences}"]))
        
def checking_weightval(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{x} is not a valid float")
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError(f"The weight must be between 0.0 and 1.0, but got {x}")
    return x

def remove_dir(dirpath):
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(Path(dirpath))
        
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
        
def output_int(df, outpath):
    df.loc[:, "chain_i"] = df.loc[:, "chain_i"].astype(str)
    df.loc[:, "chain_j"] = df.loc[:, "chain_j"].astype(str)
    df.loc[:, "base_i"] = df.loc[:, "base_i"].astype(str)
    df.loc[:, "base_j"] = df.loc[:, "base_j"].astype(str)
    df.loc[:, "score"] = df.loc[:, "score"].astype(str)
    interact = (df.loc[:, ["chain_i", "chain_j", "base_i", "base_j", "score"]]).values.tolist()

    interact = [",".join(interact[i]) for i in range(len(interact))]

    with open(outpath, "w") as f:
        f.write("\n".join(interact))

def load_basepairs(filepath):
    data = []
    with open(filepath, "r") as f:
        line = f.readline().rstrip()
        
        while line:
            
            match = bp_pattern1.match(line)
            
            if match:
                chain_i, chain_j, base_i, base_j, score = match.groups()
                data.append([chain_i, chain_j, int(base_i), int(base_j), float(score)])
            else:
                match = bp_pattern2.match(line)
                
                if match:
                    chain_i, chain_j, base_i, base_j = match.groups()
                    data.append([chain_i, chain_j, int(base_i), int(base_j), 1.0])
                else:
                    raise ValueError(f"Invalid line format: {line}")
                
            line = f.readline().rstrip()
            
    return pd.DataFrame(data, columns=["chain_i", "chain_j", "base_i", "base_j", "score"])
        
def preprocess(pdbpath):
    dist_atom = "C4'"
    bb_atoms = ["P","O5'","C5'","C4'","C3'","O3'"]
    
    pdbobj = PDBLoader(pdbpath, target_atoms = bb_atoms)
    coordinates_, info = pdbobj.get_structure()
    
    coordinates1 = torch.tensor([[coordinates_[j][atom] if atom in coordinates_[j] else [np.inf, np.inf, np.inf] for atom in frame1_atoms] for j in range(len(coordinates_))]).float().unsqueeze(0)
    coordinates2 = torch.tensor([[coordinates_[j][atom] if atom in coordinates_[j] else [np.inf, np.inf, np.inf] for atom in frame2_atoms] for j in range(len(coordinates_))]).float().unsqueeze(0)
    
    rots1, trans1, coord_mask1, maked_rots1, maked_trans1 = build_affine3d_from_coordinates(coordinates1)
    rots2, trans2, coord_mask2, maked_rots2, maked_trans2 = build_affine3d_from_coordinates(coordinates2)
    
    rots1 = torch.concat([rots1, maked_rots1.unsqueeze(1)], dim = 1)
    trans1 = torch.concat([trans1, maked_trans1.unsqueeze(1)], dim = 1)
    coord_mask1 = torch.concat([(coord_mask1.long()), torch.tensor([1]).unsqueeze(1)], dim = 1)
    
    rots2 = torch.concat([rots2, maked_rots2.unsqueeze(1)], dim = 1)
    trans2 = torch.concat([trans2, maked_trans2.unsqueeze(1)], dim = 1)
    coord_mask2 = torch.concat([(coord_mask2.long()), torch.tensor([1]).unsqueeze(1)], dim = 1)
            
    B, S, _, _ = coordinates1.size()
    
    coordinates_origin = coordinates1[:, :, 1, :]
    
    dist_map = ((coordinates_origin[..., None, :, :] - coordinates_origin[..., :, None, :]) ** 2).sum(-1)
    
    chains = [info[j]["chain_name"] for j in range(len(info))]
    pdb_idx = [info[j]["seq_idx"] for j in range(len(info))]
    pdb_nts = [info[j]["nt_name"] for j in range(len(info))]
    pdb_nt_mask = [info[j]["required_flag"] for j in range(len(info))]
    
    pdbidx2seqidx = dict(zip(pdb_idx, np.arange(S).tolist()))
    seqidx2pdbidx = dict(zip(np.arange(S).tolist(), pdb_idx))
    seqidx2nt = dict(zip(np.arange(S).tolist(), pdb_nts))
    
    return dist_map, info, [B, S], [rots1, trans1, coord_mask1], [rots2, trans2, coord_mask2], chains, pdb_idx, pdb_nts, pdb_nt_mask, pdbidx2seqidx, seqidx2pdbidx, seqidx2nt
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This package is for predicting base-pairing interactions.")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-p", "--pdb", help="Input pdb file", type=str)
    #parser.add_argument("--model", help="Directory of models", type=str)
    parser.add_argument("-bp", "--basepair", help="File of base-pair information", type=str)
    parser.add_argument("-o", "--output", help="Directory of output", required=True, type=str)
    
    maxlength=1024
    
    args = parser.parse_args()
    outpath=args.output
    bppath=args.basepair
    pdbpath = args.pdb
    
    modelpath = os.path.join(os.path.dirname(__file__), 'models')
    #args.model
    dist_map, str_info, str_size, str_fe_1, str_fe_2, chains, pdb_idx, pdb_nts, pdb_nt_mask, pdbidx2seqidx, seqidx2pdbidx, seqidx2nt = preprocess(pdbpath)
    
    if(len(list(set(chains))) > 1):
        print("Error: PDB file includes multi-chains. PDB file can be included only a chain.")
        sys.exit()
        
    chain_name = chains[0]
    seq = "".join(pdb_nts)
    
    if(bppath is None):
        bb_bp, basepairmaps, processflag = predict_bp_from_backbonestr(dist_map, str_info, str_size, str_fe_1, str_fe_2, pdb_nts, pdb_nt_mask, seqidx2pdbidx, seqidx2nt, modelpath, chain_name, f"{outpath}")
        
        output_int(bb_bp, f"{outpath}/RNAref3D_bp.txt")
        
        print("aaa")
    else:
        bb_bp = load_basepairs(bppath)
    
    #if((len(seq) < maxlength) and processflag):
        #bb_bp = predict_bp_from_initbp(basepairmaps, seq, dist_map, seqidx2pdbidx, seqidx2nt, modelpath, chain_name, outpath, maxlength=maxlength)
    #else:
        #if(os.path.isfile(f"{outpath}/RNAref3D_bp_bb.csv")):
            #os.makedirs(outpath, exist_ok = True)
            #shutil.copy2(f"{outpath}/RNAref3D_bp_bb.csv", f"{outpath}/RNAref3D_bp.csv")

    bpmap, nt2pnt = get_bpmap(bb_bp, pdbidx2seqidx, str_size[1])
    
    bb_adjstack, adjstackpairmaps = predict_adjstack_from_str(bpmap, nt2pnt, dist_map, str_info, str_size, str_fe_1, str_fe_2, pdb_nts, pdb_nt_mask, seqidx2pdbidx, seqidx2nt, modelpath, chain_name, outpath)
    bb_nonadjstack, nonadjstackpairmaps = predict_nonadjstack_from_str(bpmap, nt2pnt, dist_map, str_info, str_size, str_fe_1, str_fe_2, pdb_nts, pdb_nt_mask, seqidx2pdbidx, seqidx2nt, modelpath, chain_name, outpath)
    
    if((len(bb_adjstack) > 0) and (len(bb_nonadjstack) > 0)):
        bb_stack = pd.concat([bb_adjstack, bb_nonadjstack])
        bb_stack = bb_stack.reset_index(drop = True)
    elif(len(bb_adjstack) > 0):
        bb_stack = bb_adjstack
        bb_stack = bb_stack.reset_index(drop = True)
    elif(len(bb_nonadjstack) > 0):
        bb_stack = bb_nonadjstack
        bb_stack = bb_stack.reset_index(drop = True)
    else:
        bb_stack = pd.DataFrame([], columns = ["chain_i", "chain_j", "base_i", "base_j", "base_i_nt", "base_j_nt", "score"])
    
    output_int(bb_stack, f"{outpath}/RNAref3D_stack.txt")
    
    #os.makedirs(outpath, exist_ok = True)
    #bb_stack.to_csv(f"{outpath}/RNAref3D_stack.csv", index = False)











