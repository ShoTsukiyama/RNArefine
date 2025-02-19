import torch
import torch.nn as nn
import numpy as np
from torch.nn import functional as F
from collections import OrderedDict

class EmbedSeq(nn.Module):
    def __init__(self,seq_dim,m_dim,z_dim):
        super(EmbedSeq, self).__init__()
        self.qlinear = nn.Linear(seq_dim,z_dim)
        self.klinear = nn.Linear(seq_dim,z_dim)
        self.slinear = nn.Linear(seq_dim,m_dim)
        #self.pos = self.compute_pos().float()
        #self.pos1d=self.compute_apos()
        #self.poslinear = nn.Linear(65,z_dim)
        #self.poslinear2 = nn.Linear(14,m_dim)

    #def tocuda(self,device):
        #self.to(device)
        #self.pos.to(device)

    #def compute_apos(self,maxL=2000):
        #d = torch.arange(maxL)
        #m = 14
        #d =(((d[:,None] & (1 << np.arange(m)))) > 0).float()
        #return d

    #def compute_pos(self,maxL=2000):
        #a = torch.arange(maxL)
        #b = (a[None,:]-a[:,None]).clamp(-32,32)
        #return F.one_hot(b+32,65)

    def forward(self,seq):
        #if self.pos.device != seq.device:
            #self.pos = self.pos.to(seq.device)
        #if self.pos1d.device != seq.device:
            #self.pos1d = self.pos1d.to(seq.device)
        # msa N L D, seq L D

        N,L,D=seq.shape
        #s = self.slinear(seq)
        m = self.slinear(seq)
        #p = self.poslinear2(self.pos1d[:L])

        #m = s + p

        sq=self.qlinear(seq)
        sk=self.klinear(seq)
        z=sq[:,None,:,:]+sk[:,:,None,:]
        
        #z = z + self.poslinear( self.pos[None,:L,:L,:]  )
        return m,z








