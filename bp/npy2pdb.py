#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 17:45:09 2023

@author: tsukiyamashou
"""
import sys
import os
import copy
import pandas as pd
import numpy as np
import shutil
import string

def outpdb(coor_np,seq,savefile,start=0,end=10000,energystr=''):
    last_name=['P','C','N']
    L = len(seq)
    #wstr=[f'REMARK {str(energystr)}']
    wstr=[]
    templet='%6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s'
    count=1
    for i in range(L):
        if seq[i] in ['a','g','A','G']:
            Atom_name = [' P  '," C4'",' N9 ']
            #atoms = ['P','C4']

        elif seq[i] in ['c','u','C','U']:
            Atom_name = [' P  '," C4'",' N1 ']
        for j in range(coor_np.shape[1]):
            outs=('ATOM  ',count,Atom_name[j],seq[i],'A',i+1,coor_np[i][j][0],coor_np[i][j][1],coor_np[i][j][2],0,0,last_name[j],'')
            #outs=('ATOM  ',count,Atom_name[j],'ALA','A',i+1,coor_np[i][j][0],coor_np[i][j][1],coor_np[i][j][2],1.0,90,last_name[j],'')
            #print(outs)
            if i>=start-1 and i < end:
                wstr.append(templet % outs)

                count+=1

    wstr='\n'.join(wstr)
    wfile=open(savefile,'w')
    wfile.write(wstr)
    wfile.close()
    
#if __name__ == '__main__':
    #infile = sys.argv[1] 
    #fastafile = sys.argv[2] 
    #outfile = sys.argv[3]  

    #infile = "/Users/tsukiyamashou/Documents/Project_RNA3D/RNARefiner/RNArefiner/stacking_prediction/dataset/train_pdb_data/1a1tB/DRfold/e2e_0.npy"
    #fastafile = "/Users/tsukiyamashou/Documents/Project_RNA3D/RNARefiner/RNArefiner/stacking_prediction/dataset/train_pdb_data/1a1tB/seq.fa"
    #outfile = "/Users/tsukiyamashou/Documents/Project_RNA3D/RNARefiner/RNArefiner/stacking_prediction/dataset/train_pdb_data/1a1tB/DRfold/e2e_0.pdb"

    #coor_np = np.load(infile)
    #seq=open(fastafile).readlines()[1].strip()

    #assert len(seq) == coor_np.shape[0]
    
    #outpdb(coor_np,seq,outfile)
    
def main_npy2pdb(infile, fastafile, outfile):
    coor_np = np.load(infile)
    seq=open(fastafile).readlines()[1].strip()

    assert len(seq) == coor_np.shape[0]
    
    outpdb(coor_np,seq,outfile)
    











