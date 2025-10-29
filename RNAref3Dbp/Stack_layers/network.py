#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 17:45:09 2023

@author: tsukiyamashou
"""
import sys
import os
import copy
import functools
import einops
import math
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch
from typing import Any, List, Literal, Optional
import torch.utils.data as data
from torch.utils.data import DataLoader
from sklearn.model_selection import train_test_split
import torch
import torch.nn.functional as F
from Stack_layers.attention import MultiHeadAttention
from Stack_layers.geom_attention import GeometricReasoningOriginalImpl

class RelativePositionEmbedding(nn.Module):
    """
    Embedding layer for relative position embeddings. `bins` is the number of positions relative
    to the query position that are considered before clipping. For instance, if `bins=10`, then
    the relative position embedding will have 21 positions, [-10, 10].
    """

    def __init__(self, bins, embedding_dim, init_std=0.02):
        super().__init__()
        self.bins = bins

        self.embedding = torch.nn.Embedding(2 * bins + 3, embedding_dim)
        self.embedding.weight.data.normal_(0, init_std)

    def forward(self, diff, idx):
        """
        Input:
          query_residue_index: (B, ) tensor of source indices (dytpe=torch.long)
          key_residue_index: (B, L) tensor of target indices (dytpe=torch.long)
        Output:
          embeddings: B x L x embedding_dim tensor of embeddings
        """

        diff = diff.clamp(-self.bins, self.bins)
        diff = diff + self.bins + 1  # add 1 to adjust for padding index
        diff[idx == 0] = 0
        diff[idx == 2] = 2 * self.bins + 2
        
        output = self.embedding(diff.long())
        return output
    
class TokenEmbedding(nn.Module):
    def __init__(self, d_model):
        super().__init__()
        self.relative_positional_embedding = RelativePositionEmbedding(
            32, d_model, init_std=0.02
        )
        
        self.nt_embedding = torch.nn.Embedding(6, d_model)
        self.post_cat = torch.nn.Linear(2 * d_model, d_model)
        self.pos_scale = torch.nn.Parameter(torch.tensor(1.0))
        self.ln = torch.nn.LayerNorm(d_model)
        self.dropout = torch.nn.Dropout(p=0.1)
        
    def forward(self,
                seq_indices: torch.Tensor, 
                seqpos_rel_i: torch.Tensor, 
                seqpos_rel_j: torch.Tensor, 
                seqpos_idx_rel_i: torch.Tensor, 
                seqpos_idx_rel_j: torch.Tensor
    ):
        
        pos_i = self.relative_positional_embedding(seqpos_rel_i, seqpos_idx_rel_i)
        pos_j = self.relative_positional_embedding(seqpos_rel_j, seqpos_idx_rel_j)
        
        z = self.nt_embedding(seq_indices)
        
        z_i = z + self.pos_scale * pos_i
        z_j = z + self.pos_scale * pos_j
        
        zcat = torch.concat([z_i, z_j], dim = 2)
        z = self.post_cat(self.dropout(zcat))
        z = self.ln(z)
        
        return z
    
def swiglu_correction_fn(expansion_ratio: float, d_model: int) -> int:
    # set hidden dimesion to nearest multiple of 256 after expansion ratio
    return int(((expansion_ratio * d_model) + 255) // 256 * 256)


class SwiGLU(nn.Module):
    """
    SwiGLU activation function as an nn.Module, allowing it to be used within nn.Sequential.
    This module splits the input tensor along the last dimension and applies the SiLU (Swish)
    activation function to the first half, then multiplies it by the second half.
    """

    def __init__(self):
        super(SwiGLU, self).__init__()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x1, x2 = x.chunk(2, dim=-1)
        return F.silu(x1) * x2


def swiglu_ln_ffn(d_model: int, expansion_ratio: float, bias: bool):
    return nn.Sequential(
        nn.LayerNorm(d_model),
        nn.Linear(
            d_model, swiglu_correction_fn(expansion_ratio, d_model) * 2, bias=bias
        ),
        SwiGLU(),
        nn.Linear(swiglu_correction_fn(expansion_ratio, d_model), d_model, bias=bias),
    )


def gelu_ln_ffn(d_model: int, expansion_ratio: float, bias: bool):
    hidden_dim = int(expansion_ratio * d_model)
    return nn.Sequential(
        nn.LayerNorm(d_model),
        nn.Linear(d_model, hidden_dim, bias=bias),
        nn.GELU(),
        nn.Linear(hidden_dim, d_model, bias=bias),
    )
    
class UnifiedTransformerBlock(nn.Module):
    """
    A unified transformer block that can optionally incorporate geometric attention.

    This class defines a transformer block that can be configured to use geometric attention
    alongside the standard multi-head attention mechanism. It is designed to be a flexible
    component of transformer-based models, allowing for the integration of geometric reasoning.

    Parameters
    ----------
    d_model : int
        The dimensionality of the input and output features of the transformer block.
    n_heads : int
        The number of attention heads in the multi-head attention mechanism.
    n_layers : int
        The number of layers in the transformer block.
    use_geom_attn : bool, optional
        Whether to use geometric attention in addition to the standard multi-head attention. Defaults to False.
    v_heads : int, optional
        The number of heads to use for the geometric attention mechanism, if enabled. Must be specified if `use_geom_attn` is True.
    """

    def __init__(
        self,
        d_model: int,
        n_heads: int,
        use_geom_attn: bool = False,
        use_plain_attn: bool = True,
        use_flash_attn: bool = False,
        v_heads: Optional[int] = None,
        bias: bool = False,
        expansion_ratio: float = 4.0,
        residue_scaling_factor: float = 1,
        mask_and_zero_frameless: bool = False,
        qk_layernorm: bool = True,
        ffn_type: str = "swiglu",  # swiglu | gelu
    ):
        super().__init__()
        self.use_plain_attn = use_plain_attn
        if self.use_plain_attn:
            #if use_flash_attn:
                #self.attn = FlashMultiHeadAttention(
                    #d_model, n_heads, bias, qk_layernorm=qk_layernorm
                #)
            #else:
            self.attn = MultiHeadAttention(
                d_model, n_heads, bias, qk_layernorm=qk_layernorm
            )
            
            
        self.use_geom_attn = use_geom_attn
        if self.use_geom_attn:
            if v_heads is None:
                raise ValueError("v_heads must be specified when use_geom_attn is True")
            self.geom_attn = GeometricReasoningOriginalImpl(
                c_s=d_model,
                v_heads=v_heads,
                bias=bias,
                mask_and_zero_frameless=mask_and_zero_frameless,
            )
        if ffn_type == "swiglu":
            self.ffn = swiglu_ln_ffn(d_model, expansion_ratio, bias)
        elif ffn_type == "gelu":
            self.ffn = gelu_ln_ffn(d_model, expansion_ratio, bias)
        else:
            raise ValueError(f"Unknown ffn_type: {ffn_type}")
        self.scaling_factor = residue_scaling_factor

    def forward(
        self,
        x: torch.Tensor,
        sequence_id: torch.Tensor,
        rot: torch.Tensor,
        trans: torch.Tensor,
        bpmap: torch.Tensor,
        frames_mask: torch.Tensor,
        chain_id: torch.Tensor,
    ) -> torch.Tensor:
        """
        Forward pass for the UnifiedTransformerBlock.

        Parameters
        ----------
        x : torch.Tensor[float]
            Input tensor to the transformer block, typically the output from the previous layer.
        sequence_id : torch.Tensor[int]
            Tensor containing sequence IDs for each element in the batch, used for attention masking.
        frames : Affine3D
            Affine3D containing geometric frame information for geometric attention.
        frames_mask : torch.Tensor[bool]
            Boolean mask tensor indicating valid frames for geometric attention.
        chain_id : torch.Tensor[int]
            Tensor containing chain IDs for each element, used for attention masking in geometric attention.

        Returns
        -------
        torch.Tensor[float]
            The output tensor after applying the transformer block operations.
        """
        if self.use_plain_attn:
            r1 = self.attn(x, sequence_id)
            x = x + r1 / self.scaling_factor

        if self.use_geom_attn:
            r2 = self.geom_attn(x, rot, trans, bpmap, frames_mask, sequence_id, chain_id)
            x = x + r2 / self.scaling_factor

        r3 = self.ffn(x) / self.scaling_factor
        x = x + r3

        return x

    
class TransformerStack(nn.Module):
    """
    A stack of transformer blocks used in the ESM-3 model. Each block is a UnifiedTransformerBlock,
    which can either be geometric attention or standard multi-head attention.

    Args:
        d_model (int): The dimensionality of the input and output feature vectors.
        n_heads (int): The number of attention heads.
        v_heads (int): The number of voting heads.
        n_layers (int): The number of transformer blocks in the stack.
        n_layers_geom (int, optional): The number of transformer blocks that use geometric attention.
        scale_residue (bool, optional): Whether to scale the residue connections in each transformer block.
        mask_and_zero_frameless (bool, optional): Whether to mask and zero frameless positions in the input.
            Only applies in the geometric attention blocks, which is conditioned on the structure
    """

    def __init__(
        self,
        d_model: int,
        n_heads: int,
        v_heads: int,
        n_layers: int,
        n_layers_geom: int = 1,
        scale_residue: bool = True,
        mask_and_zero_frameless: bool = False,
        bias: bool = False,
        qk_layernorm: bool = True,
        ffn_type: str = "swiglu",  # swiglu | gelu
        expansion_ratio: float = 8 / 3,
        use_flash_attn: bool = False,
    ):
        super().__init__()
        self.blocks = nn.ModuleList(
            [
                UnifiedTransformerBlock(
                    d_model,
                    n_heads,
                    v_heads=v_heads,
                    use_geom_attn=i < n_layers_geom,
                    use_flash_attn=use_flash_attn,
                    residue_scaling_factor=(
                        math.sqrt(n_layers / 36) if scale_residue else 1.0
                    ),
                    expansion_ratio=expansion_ratio,
                    mask_and_zero_frameless=mask_and_zero_frameless,
                    bias=bias,
                    qk_layernorm=qk_layernorm,
                    ffn_type=ffn_type,
                )
                for i in range(n_layers)
            ]
        )
        self.norm = nn.LayerNorm(d_model, bias=False)

    def forward(
        self,
        x: torch.Tensor,
        sequence_id: Optional[torch.Tensor] = None,
        rot: Optional[torch.Tensor] = None,
        trans: Optional[torch.Tensor] = None,
        bpmap: Optional[torch.Tensor] = None,
        affine_mask: Optional[torch.Tensor] = None,
        chain_id: Optional[torch.Tensor] = None,
    ) -> tuple[torch.Tensor, torch.Tensor, list[torch.Tensor]]:
        """
        Forward pass of the TransformerStack.

        Args:
            x (torch.Tensor): The input tensor of shape (batch_size, sequence_length, d_model).
            sequence_id (torch.Tensor): The sequence ID tensor of shape (batch_size, sequence_length).
            affine (Affine3D | None): The affine transformation tensor or None.
            affine_mask (torch.Tensor | None): The affine mask tensor or None.
            chain_id (torch.Tensor): The protein chain tensor of shape (batch_size, sequence_length).
                Only used in geometric attention.

        Returns:
            post_norm: The output tensor of shape (batch_size, sequence_length, d_model).
            pre_norm: The embedding of shape (batch_size, sequence_length, d_model).
        """
        *batch_dims, _ = x.shape
        if chain_id is None:
            chain_id = torch.ones(size=batch_dims, dtype=torch.int64, device=x.device)
        hiddens = []
        for block in self.blocks:
            x = block(x, sequence_id, rot, trans, bpmap, affine_mask, chain_id)
            hiddens.append(x)
        return self.norm(x), x, hiddens

        
class GeometricEncoderStack(TransformerStack):
    def __init__(self, d_model, n_heads, v_heads, n_layers):
        super().__init__(d_model, n_heads, v_heads, 0)
        self.blocks = nn.ModuleList(
            [
                UnifiedTransformerBlock(
                    d_model,
                    n_heads,
                    v_heads=v_heads,
                    use_geom_attn=True,
                    use_plain_attn=False,
                    expansion_ratio=4,
                    bias=True,
                )
                for i in range(n_layers)
            ]
        )
        self.norm = nn.Identity()
        
class Outputlayer(nn.Module):
    def __init__(self, in_dim=2048, hidden=None, p_drop=0.1, use_ln=True):
        super().__init__()
        layers = []
        if hidden: 
            layers += [nn.Linear(in_dim, hidden)]
            if use_ln:
                layers += [nn.LayerNorm(hidden)] 
            layers += [nn.GELU()]
            if p_drop > 0:
                layers += [nn.Dropout(p_drop)]
            in_dim = hidden
        self.body = nn.Sequential(*layers)
        self.out = nn.Linear(in_dim, 1)
        nn.init.zeros_(self.out.weight)
        nn.init.zeros_(self.out.bias)

    def forward(self, x):       
        x = self.body(x)             
        y = self.out(x)              
        return y

class BPDetectNet(nn.Module):
    def __init__(self, d_model = 256, n_heads = 1, v_heads = 128, n_layers = 5, d_out = 128, n_codes = 4096):
        super().__init__()
        self.token_emb = TokenEmbedding(d_model)
        
        self.transformer = GeometricEncoderStack(d_model, n_heads, v_heads, n_layers)
        
        self.outlayer = Outputlayer(d_model*2, hidden=d_model)
        
    def forward(
        self, 
        seq_indices: torch.Tensor, 
        rot: torch.Tensor, 
        trans: torch.Tensor, 
        bpmap: torch.Tensor, 
        seqpos_rel_i: torch.Tensor, 
        seqpos_rel_j: torch.Tensor, 
        mask: torch.Tensor, 
        seqpos_idx_rel_i: torch.Tensor, 
        seqpos_idx_rel_j: torch.Tensor
    ):
        B, L = seq_indices.size()
        z = self.token_emb(seq_indices, seqpos_rel_i, seqpos_rel_j, seqpos_idx_rel_i, seqpos_idx_rel_j)
        z, _, _ = self.transformer(z, None, rot, trans, bpmap, mask, None)
        z = z[:, :2, :]
        z = z.view(B, -1)
        out = self.outlayer(z)
        
        return out
        
        
        
        














    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
