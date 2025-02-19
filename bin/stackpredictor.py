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

sys.path.append(os.path.join(os.path.dirname(__file__), '../stack'))
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
    parser.add_argument("--pdb", help="Path to a pdb file", required=True, type=str)
    parser.add_argument("--bp", help="Path to a file with base pair information", required=True, type=str)
    parser.add_argument("-o", "--output", help="Path to the pdb file", required=True, type=str)
    
    args = parser.parse_args()
    pdbpath=args.pdb
    bppath=args.bp
    outpath=args.output
    os.makedirs(outpath, exist_ok = True)
    
    predict.stacking_prediction_main(pdbpath, bppath, outpath)
















