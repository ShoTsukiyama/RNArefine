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

def fp32_autocast_context(device_type: str):
    """
    Returns an autocast context manager that disables downcasting by AMP.

    Args:
        device_type: The device type ('cpu' or 'cuda')

    Returns:
        An autocast context manager with the specified behavior.
    """
    if device_type == "cpu":
        return torch.amp.autocast(device_type, enabled=False)
    elif device_type == "cuda":
        return torch.amp.autocast(device_type, dtype=torch.float32)
    else:
        raise ValueError(f"Unsupported device type: {device_type}")

def graham_schmidt(
    neg_x_axis: torch.Tensor,
    origin: torch.Tensor,
    xy_plane: torch.Tensor,
    eps: float = 1e-10,):
    
    with fp32_autocast_context(origin.device.type):
        x_axis = origin - neg_x_axis
        xy_plane = xy_plane - origin
    
        e1 = xy_plane
        denom = torch.sqrt((x_axis**2).sum(dim=-1, keepdim=True) + eps)
        x_axis = x_axis / denom
        dot = (x_axis * e1).sum(dim=-1, keepdim=True)
        e1 = e1 - x_axis * dot
        denom = torch.sqrt((e1**2).sum(dim=-1, keepdim=True) + eps)
        e1 = e1 / denom
        e2 = torch.cross(x_axis, e1, dim=-1)

        rots = torch.stack([x_axis, e1, e2], dim=-1)
    
    return origin, rots

def build_affine3d_from_coordinates(
    coords: torch.Tensor,  # (N, CA, C).
):
    _MAX_SUPPORTED_DISTANCE = 1e6
    coord_mask = torch.all(
        torch.all(torch.isfinite(coords) & (coords < _MAX_SUPPORTED_DISTANCE), dim=-1),
        dim=-1,
    )
    
    def atom3_to_backbone_affine(bb_positions: torch.Tensor):
        P, C4, C3 = bb_positions.unbind(dim=-2)
        return graham_schmidt(P, C4, C3)
    
    coords = coords.clone().float()
    coords[~coord_mask] = 0
    
    #分子ごとに重心を計算 (N, L, 3, 3) -> (N, 3, 3)
    average_per_p_c4_c3 = coords.masked_fill(~coord_mask[..., None, None], 0).sum(1) / (
        coord_mask.sum(-1)[..., None, None] + 1e-8
    )
    
    trans, rots = atom3_to_backbone_affine(
        average_per_p_c4_c3.float()
    )
    
    B, S, _, _ = coords.shape
    
    assert isinstance(B, int)
    assert isinstance(S, int)
    affine_rot_mats = rots.flatten(-2)[..., None, :].expand(B, S, 9)
    affine_trans = trans[..., None, :].expand(B, S, 3)
    
    identity_rot = torch.eye(3, device = coords.device, dtype=torch.float32)
    identity_rot = identity_rot.view(*[1 for _ in range(2)], 3, 3)
    identity_rot = identity_rot.expand(B, S, -1, -1)
    
    affine_rot_mats = affine_rot_mats.where(
        coord_mask.any(-1)[..., None, None], identity_rot.flatten(-2)
    )
    
    trans_pos, rots_pos = atom3_to_backbone_affine(coords.float())
    
    final_rots  = torch.where(coord_mask[..., None, None], rots_pos,  affine_rot_mats.unflatten(-1, (3, 3)))  # (B,S,3,3)
    final_trans = torch.where(coord_mask[..., None],      trans_pos, affine_trans) 
    
    return final_rots, final_trans, coord_mask, rots, trans

def build_trans3d_from_coordinates(
    coords: torch.Tensor,  # (N, CA, C).
):
    
    B, S, _ = coords.shape
    device = coords.device
    
    _MAX_SUPPORTED_DISTANCE = 1e6
    coord_mask = torch.all(
        torch.isfinite(coords) & (coords < _MAX_SUPPORTED_DISTANCE), dim=-1
    )
    
    final_trans = coords.clone().float()
    
    tmp_coords = final_trans.masked_fill(~coord_mask.unsqueeze(-1), 0.0)
    valid_count = coord_mask.sum(dim=-1, keepdim=True) + 1e-8
    global_trans = tmp_coords.sum(dim=1) / valid_count # (B, 3)
    
    final_trans = torch.where(
        coord_mask.unsqueeze(-1), 
        final_trans, 
        global_trans.unsqueeze(1)
    )
    
    return final_trans, coord_mask, global_trans
    
    """
    P, C4, C3 = coords.unbind(dim=-2)
    
    minus_d = (P-C4).norm(dim = 2).unsqueeze(-1)
    x = torch.concat([-minus_d, torch.zeros_like(minus_d), torch.zeros_like(minus_d)], dim = 2)
    
    print(final_rots.unsqueeze(-3).size())
    print(x.unsqueeze(-1).unsqueeze(2).size())
    print(final_trans.unsqueeze(-2).size())
    
    x_global = (final_rots.unsqueeze(-3) @ x.unsqueeze(-1).unsqueeze(2)).squeeze(-1) + final_trans.unsqueeze(-2)
    
    Rinv = final_rots.transpose(-1, -2)
    
    print(Rinv.unsqueeze(-3).size())
    print(x_global.unsqueeze(-1).size())
    print(final_trans.unsqueeze(-2).size())
    
    attn_out = (Rinv.unsqueeze(-3) @ (x_global.unsqueeze(-1) - final_trans.unsqueeze(-2).unsqueeze(-1))).squeeze(-1)
    
    print(attn_out.squeeze(2).size())
    print(x.size())
    
    print((attn_out.squeeze(2) - x))
    
    sys.exit()
    return final_rots, final_trans, coord_mask, rots, trans
    """
    
    #P, C4, C3 = average_per_p_c4_c3.unbind(dim=-2)
    
    #origin, rots = affine_from_average
    
    # P, C4 は (B,3)
    #d = torch.linalg.norm(C4 - P, dim=-1, keepdim=True)  # (B,1)

    # ローカル座標（列ベクトル）：P は -x 方向
    #p_local = torch.stack([-d.squeeze(-1), torch.zeros_like(d.squeeze(-1)), torch.zeros_like(d.squeeze(-1))], dim=-1)  # (B,3)
    #p_local = p_local[..., None]  # (B,3,1)

    # 変換：world = origin + R @ local
    #p_world = rots @ p_local + origin[..., None]  # (B,3,1)
    #p_world = p_world.squeeze(-1)                 # (B,3)
    
    #print(p_world)
    #print(P)










































