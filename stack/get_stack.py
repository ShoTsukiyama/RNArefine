import os
import sys
import re
import numpy as np
import pandas as pd
import subprocess
import shutil

def get_bp_clarna_bp(pdbpath, seq):
    seqlen=len(seq)
    bin_path = shutil.which("rna_clarna_run.py")
    
    command = [bin_path, "-bps", "-ipdb", pdbpath]
    result = subprocess.run(command, capture_output=True, text=True)
    
    line = result.stdout.split("\n")
    pattern = r'([A-Za-z0-9])\s+(\d+)\s+([A-Za-z0-9])\s+(\d+)\s+bp\s+([A-Za-z])\s+([A-Za-z])\s+(\S+)\s+(\S+)'

    pairs=[]
    for i in range(len(line)):
        match = re.match(pattern, line[i])
        if match:
            chain_i, seqid1, chain_j, seqid2,nt1, nt2,bptype,score=match.groups()
            
            seqid1=int(seqid1)
            seqid2=int(seqid2) 

            if(seqid1 < seqid2): pairs.append([seqid1, seqid2])
            else: pairs.append([seqid2, seqid1])
            
    return pd.DataFrame(pairs, columns = ["base_i", "base_j"])

def get_bp_clarna_stack(pdbpath, seq):
    seqlen=len(seq)
    bin_path = shutil.which("rna_clarna_run.py")
    
    command = [bin_path, "-stack", "-ipdb", pdbpath]
    result = subprocess.run(command, capture_output=True, text=True)
    
    line = result.stdout.split("\n")
    pattern = r'([A-Za-z0-9])\s+(\d+)\s+([A-Za-z0-9])\s+(\d+)\s+bp\s+([A-Za-z])\s+([A-Za-z])\s+(\S+)\s+(\S+)'

    pairs=[]
    for i in range(len(line)):
        match = re.match(pattern, line[i])
        if match:
            chain_i, seqid1, chain2, seqid2,nt1, nt2,bptype,score=match.groups()
            
            seqid1=int(seqid1)
            seqid2=int(seqid2)

            if(seqid1 < seqid2): pairs.append([seqid1, seqid2])
            else: pairs.append([seqid2, seqid1])
            
    return pd.DataFrame(pairs, columns = ["base_i", "base_j"])





