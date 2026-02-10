from math import sqrt
import sys
import torch
from einops import rearrange
from torch import nn
from torch.nn import functional as F


class GeometricReasoningOriginalImpl(nn.Module):
    def __init__(
        self,
        c_s: int,
        v_heads: int,
        num_vector_messages: int = 1,
        mask_and_zero_frameless: bool = True,
        divide_residual_by_depth: bool = False,
        bias: bool = False,
    ):
        super().__init__()
        self.c_s = c_s
        self.v_heads = v_heads
        self.num_vector_messages = num_vector_messages
        self.mask_and_zero_frameless = mask_and_zero_frameless

        self.s_norm = nn.LayerNorm(c_s, bias=bias)
        dim_proj = (
            2 * self.v_heads * 3 + self.v_heads * 3 * self.num_vector_messages
        )
        self.proj = nn.Linear(c_s, dim_proj, bias=bias)
        channels_out = self.v_heads * 3 * self.num_vector_messages
        self.out_proj = nn.Linear(channels_out, c_s, bias=bias)

        self.bp_proj = nn.Embedding(2, self.v_heads)
        # The basic idea is for some attention heads to pay more or less attention to rotation versus distance,
        # as well as to control the sharpness of the softmax (i.e., should this head only attend to those residues
        # very nearby or should there be shallower dropoff in attention weight?)
        self.distance_scale_per_head = nn.Parameter(torch.zeros((self.v_heads)))
        #self.rotation_scale_per_head = nn.Parameter(torch.zeros((self.v_heads)))
        self.basepair_scale_per_head = nn.Parameter(torch.zeros((self.v_heads)))

        initial_centers = torch.linspace(0.0, 12.0, steps=self.v_heads)
        self.centers = nn.Parameter(initial_centers.view(1, self.v_heads, 1, 1))
        self.sigmas = nn.Parameter(torch.ones(self.v_heads).view(1, self.v_heads, 1, 1) * 1.5)

    def forward(self, s, trans, bpmap, affine_mask, sequence_id, chain_id):
        if sequence_id is None:
            sequence_id = torch.zeros_like(s[..., 0], dtype=torch.int64)
        attn_bias = sequence_id.unsqueeze(-1) == sequence_id.unsqueeze(-2)
        attn_bias = attn_bias.unsqueeze(1).float()
        attn_bias = attn_bias.masked_fill(
            ~affine_mask[:, None, None, :], torch.finfo(attn_bias.dtype).min
        )
        chain_id_mask = chain_id.unsqueeze(1) != chain_id.unsqueeze(2)
        attn_bias = attn_bias.masked_fill(
            chain_id_mask.unsqueeze(1), torch.finfo(s.dtype).min
        )

        ns = self.s_norm(s)
        vec_dist = self.proj(ns)

        _p_dist = rearrange(vec_dist, "... (h c) -> ... h c", c=3)

        q_dist_local, k_dist_local, v_local = _p_dist.split(
            [self.v_heads, self.v_heads, self.v_heads * self.num_vector_messages],
            dim=-2,
        )
        q_dist_local = q_dist_local + trans.unsqueeze(-2)
        k_dist_local = k_dist_local + trans.unsqueeze(-2)

        q_dist_local = rearrange(q_dist_local, "b s h d -> b h s 1 d")
        k_dist_local = rearrange(k_dist_local, "b s h d -> b h 1 s d")
        v_local = rearrange(
            v_local, "b s (h m) d -> b h s (m d)", m=self.num_vector_messages
        )

        #distance = (q_dist_local - k_dist_local).norm(dim=-1) / sqrt(3)
        #diff = distance - self.centers.view(1, self.v_heads, 1, 1)
        #distance_term = - (diff ** 2) / (2 * (F.softplus(self.sigmas) + 1e-5) ** 2)
        
        distance = (q_dist_local - k_dist_local).norm(dim=-1) / sqrt(3)
        diff = (distance - self.centers.view(1, self.v_heads, 1, 1)).abs()
        scale = F.softplus(self.sigmas) + 1e-5
        distance_term = - diff / scale

        distance_term_weight = rearrange(
            F.softplus(self.distance_scale_per_head), "h -> h 1 1"
        )
        
        basepair_term_weight = rearrange(
            F.softplus(self.basepair_scale_per_head), "h -> h 1 1"
        )
        
        bp_term = self.bp_proj(bpmap).permute(0, 3, 1, 2)

        attn_weight = (
            bp_term * basepair_term_weight + distance_term * distance_term_weight
        )

        if attn_bias is not None:
            attn_weight = attn_weight + attn_bias

        attn_weight = torch.softmax(attn_weight, dim=-1)

        attn_out = attn_weight.matmul(v_local)

        _v_g = rearrange(attn_out, "b h s (m d) -> b s (h m) d", m=self.num_vector_messages)

        #attn_out = _v_g - trans.unsqueeze(-2)
        attn_out = _v_g

        attn_out = rearrange(
            attn_out, "b s (h m) d -> b s (h m d)", m=self.num_vector_messages
        )
        if self.mask_and_zero_frameless:
            attn_out = attn_out.masked_fill(~affine_mask[..., None], 0.0)
        s = self.out_proj(attn_out)

        return s

