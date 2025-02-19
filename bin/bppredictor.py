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

sys.path.append(os.path.join(os.path.dirname(__file__), '../bp'))
import bppredict
import npy2pdb

sys.path.append(os.path.join(os.path.dirname(__file__), '../bp/scripts'))
import Feature

sys.path.append(os.path.join(os.path.dirname(__file__), '../bp/DeepE2EPotential'))
import predict

rna_codes = {
    "A": "A", "C": "C", "G": "G", "U": "U", "T": "U",
    "DA": "A", "DC": "C", "DG": "G", "DT": "U", "DU": "U"
}

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This package is for predicting base-pairing interactions.")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pdb", help="Path to a pdb file", type=str)
    group.add_argument("--fasta", help="Path to a fasta file", type=str)
    parser.add_argument("-o", "--output", help="Path to the pdb file", required=True, type=str)
    parser.add_argument("-w", "--weight", help="Weight of the secondary structure predicted by the ResNet-based neural network", type=checking_weightval, default=0.4)
    
    args = parser.parse_args()
    outpath=args.output
    os.makedirs(f"{outpath}/bp_RNAoptimizer", exist_ok = True)
    
    weight=args.weight
    
    if(args.pdb is not None):
        pdbpath=args.pdb
        
        #getting sequence
        rna_sequences = get_seq(pdbpath)
        if(len(rna_sequences) > 1):
            raise ValueError("The pdb file containts multi chains. Specify the pdb file with a single chain structure.")
            sys.exit()
        seq_id = str(rna_sequences[0].id)
        seq = str(rna_sequences[0].seq)
        output_fasta(f"{outpath}/bp_RNAoptimizer/seq.fasta", seq_id, seq)
        
        fastapath=f"{outpath}/bp_RNAoptimizer/seq.fasta"
        
    elif(args.fasta is not None):
        fastapath=args.fasta
        
    #secondary structure prediction
    Feature.predict_ss(fastapath, f"{outpath}/bp_RNAoptimizer/seq_ss")
        
    #predict e2e structures
    predict.predict_e2e_str(fastapath, f"{outpath}/bp_RNAoptimizer/seq_ss.npy", f"{outpath}/bp_RNAoptimizer/e2e")
        
    #npy to pdb
    templatepdbpath=[]
    for i in range(6):
        npy2pdb.main_npy2pdb(f"{outpath}/bp_RNAoptimizer/e2e_{i}.npy", fastapath, f"{outpath}/bp_RNAoptimizer/e2e_{i}.pdb")
        templatepdbpath.append(f"{outpath}/bp_RNAoptimizer/e2e_{i}.pdb")
        
    if(args.pdb is not None):
        bppredict.bppredictmainpdb(pdbpath, templatepdbpath, f"{outpath}/bp.list", weight)
    elif(args.fasta is not None):
        bppredict.bppredictmainfasta(templatepdbpath, f"{outpath}/bp.list")

    remove_dir(f"{outpath}/bp_RNAoptimizer")
















