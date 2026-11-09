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
import random
import shutil
import torch
import joblib
import pickle
import gemmi
from typing import Any, List, Literal, Optional

class PDBLoader():
    def __init__(
        self,
        pdbpath: str,
        target_model: Optional[int] = 0,
        target_atoms: Optional[dict[str, str]] = None
    ):
        
        if target_atoms is None:
            target_atoms = {
                "A": ["P","O5'","C5'","C4'","C3'","O3'"],
                "C": ["P","O5'","C5'","C4'","C3'","O3'"],
                "G": ["P","O5'","C5'","C4'","C3'","O3'"],
                "U": ["P","O5'","C5'","C4'","C3'","O3'"],
            }
            
        self.pdbpath = pdbpath
        self.target_model = target_model
        self.required = target_atoms
        
        self.coordinates, self.info = [], []
        
        self.structure = gemmi.read_structure(pdbpath)
        
        self.load_pdb()
        
    def normalize_atom_name(self, name: str):
        n = name.strip().upper().replace("*", "'")
        if n == "O1P": n = "OP1"
        if n == "O2P": n = "OP2"
        return n
    
    def is_rna_residue(self, res: gemmi.Residue):
        name = res.name.strip().upper()
        if name in {"A","C","G","U"}:
            return True
        atom_names = {self.normalize_atom_name(a.name) for a in res}
        return "O2'" in atom_names
    
    def residue_has_all_backbone_atoms(self, res: gemmi.Residue, required_set: set):
        atoms_present = {self.normalize_atom_name(a.name) for a in res}
        return required_set.issubset(atoms_present)
        
    def load_pdb(self):
        model = self.structure[self.target_model]
        
        info = []
        for ch in model:
            rna_res = [res for res in ch if res.het_flag == 'A' and self.is_rna_residue(res)]
            
            if not rna_res:
                continue
            
            for res in rna_res:
                ok = self.residue_has_all_backbone_atoms(res, set(self.required[res.name]))
                seqid = res.seqid.num
                
                self.coordinates.append({a.name: [a.pos.x, a.pos.y, a.pos.z] for a in res if(a.name in self.required[res.name])})
                
                self.info.append({"required_flag": ok, "chain_name": ch.name, "seq_idx": seqid, "nt_name": res.name})
                
    def get_structure(self):
        return self.coordinates, self.info
            
        































