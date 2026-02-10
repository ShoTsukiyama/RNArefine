from numpy import select
import torch
from torch import nn
from torch.nn import functional as F
import math

attdrop = True

class Linear(nn.Module):
    def __init__(self,dim_in,dim_out):
        super(Linear,self).__init__()
        self.linear = nn.Linear(dim_in,dim_out)
    def forward(self,x):
        x = self.linear(x)
        return x


class LinearNoBias(nn.Module):
    def __init__(self,dim_in,dim_out):
        super(LinearNoBias,self).__init__()
        self.linear = nn.Linear(dim_in,dim_out,bias=False)
    def forward(self,x):
        x = self.linear(x)
        return x

def DropAtt(x,dim,droprate = 0.25):
    shapes = x.shape
    L = shapes[dim]
    num_of_dim = len(shapes)
    randmask = torch.rand(L).to(x.device)
    themask = randmask<droprate
    randmask[themask] = 0
    randmask[~themask] = 1
    for i in range(dim):
        randmask = randmask[None]
    for i in range(dim,num_of_dim-1):
        randmask = randmask[...,None]
    x = x * randmask
    x = x /(x.sum(dim=dim,keepdims=True) + 1e-06)
    return x

class TriOut(nn.Module):
    def __init__(self,z_dim,c=32):
        super(TriOut,self).__init__()
        self.z_dim = z_dim
        self.norm =nn.LayerNorm(z_dim)
        self.onorm =nn.LayerNorm(c)
        self.alinear = Linear(z_dim,c)
        self.blinear = Linear(z_dim,c)
        self.aglinear = Linear(z_dim,c)
        self.bglinear = Linear(z_dim,c)
        self.glinear = Linear(z_dim,z_dim)
        self.olinear = Linear(c,z_dim)

    def forward(self,z):
        z = self.norm(z)
        a = self.alinear(z) * torch.sigmoid(self.aglinear(z))
        b = self.blinear(z) * torch.sigmoid(self.bglinear(z))
        o = torch.einsum('ilc,jlc->ijc',a,b)
        o = self.onorm(o)
        o = self.olinear(o)
        o = o * torch.sigmoid(self.glinear(z))
        return o

class TriIn(nn.Module):
    def __init__(self,z_dim,c=32):
        super(TriIn,self).__init__()
        self.z_dim = z_dim
        self.norm =nn.LayerNorm(z_dim)
        self.onorm =nn.LayerNorm(c)
        self.alinear = Linear(z_dim,c)
        self.blinear = Linear(z_dim,c)
        self.aglinear = Linear(z_dim,c)
        self.bglinear = Linear(z_dim,c)
        self.glinear = Linear(z_dim,z_dim)
        self.olinear = Linear(c,z_dim)

    def forward(self,z):
        z = self.norm(z)
        a = self.alinear(z) * torch.sigmoid(self.aglinear(z))
        b = self.blinear(z) * torch.sigmoid(self.bglinear(z))
        o = torch.einsum('lic,ljc->ijc',a,b)
        o = self.onorm(o)
        o = self.olinear(o)
        o = o * torch.sigmoid(self.glinear(z))
        return o

class TriAttStart(nn.Module):
    def __init__(self,z_dim,N_head=4,c=8):
        super(TriAttStart,self).__init__()
        self.z_dim = z_dim
        self.N_head = N_head
        self.c = c
        self.sq_c = 1/math.sqrt(c)
        self.norm=nn.LayerNorm(z_dim)
        self.qlinear = Linear(z_dim,c*N_head)
        self.klinear = Linear(z_dim,c*N_head)
        self.vlinear = Linear(z_dim,c*N_head)
        self.blinear = Linear(z_dim,N_head)
        self.glinear = Linear(z_dim,c*N_head)
        self.olinear = Linear(c*N_head,z_dim)
        
        #self.bplinear = Linear(2,N_head)

    def forward(self,z_, bp):
        L1,L2,D=z_.shape
        z = self.norm(z_)
        q = self.qlinear(z).reshape(L1,L2,self.N_head,self.c)
        k = self.klinear(z).reshape(L1,L2,self.N_head,self.c)
        v = self.vlinear(z).reshape(L1,L2,self.N_head,self.c)
        b = self.blinear(z)
        att = torch.einsum('blhc,bkhc->blkh',q,k)*self.sq_c + b[None,:,:,:] #+ self.bplinear(bp)[None,:,:,:]
        att = F.softmax(SignedSqrt(att),dim=2)
        if attdrop:
            if self.training:
                att = DropAtt(att,dim=2)
        o = torch.einsum('blkh,bkhc->blhc',att,v)
        o = (torch.sigmoid(self.glinear(z).reshape(L1,L2,self.N_head,self.c)) * o).reshape(L1,L2,-1)
        o = self.olinear(o)
        return o
    
class TriAttEnd(nn.Module):
    def __init__(self,z_dim,N_head=4,c=8):
        super(TriAttEnd,self).__init__()
        self.z_dim = z_dim
        self.N_head = N_head
        self.c = c
        self.sq_c = 1/math.sqrt(c)
        self.norm=nn.LayerNorm(z_dim)
        self.qlinear = Linear(z_dim,c*N_head)
        self.klinear = Linear(z_dim,c*N_head)
        self.vlinear = Linear(z_dim,c*N_head)
        self.blinear = Linear(z_dim,N_head)
        self.glinear = Linear(z_dim,c*N_head)
        self.olinear = Linear(c*N_head,z_dim)
        
        #self.bplinear = Linear(2,N_head)

    def forward(self,z_, bp):
        L1,L2,D=z_.shape
        z = self.norm(z_)
        q = self.qlinear(z).reshape(L1,L2,self.N_head,self.c)
        k = self.klinear(z).reshape(L1,L2,self.N_head,self.c)
        v = self.vlinear(z).reshape(L1,L2,self.N_head,self.c)
        b = self.blinear(z)
        att = torch.einsum('blhc,kbhc->blkh',q,k)*self.sq_c + b[None,:,:,:].permute(0,2,1,3) #+ self.bplinear(bp)[None,:,:,:]
        att = F.softmax(SignedSqrt(att),dim=2)
        if attdrop:
            if self.training:
                att = DropAtt(att,dim=2)
        o = torch.einsum('blkh,klhc->blhc',att,v)
        o = (torch.sigmoid(self.glinear(z).reshape(L1,L2,self.N_head,self.c)) * o).reshape(L1,L2,-1)
        o = self.olinear(o)
        return o
    def forward2(self,z_):
        z = z_.permute(1,0,2)
        L1,L2,D=z_.shape
        z = self.norm(z_)
        q = self.qlinear(z).reshape(L1,L2,self.N_head,self.c)
        k = self.klinear(z).reshape(L1,L2,self.N_head,self.c)
        v = self.vlinear(z).reshape(L1,L2,self.N_head,self.c)
        b = self.blinear(z)
        att = torch.einsum('blhc,bkhc->blkh',q,k)*self.sq_c + b[None,:,:,:]
        att = F.softmax(att,dim=2)
        if attdrop:
            if self.training:
                att = DropAtt(att,dim=2)
        o = torch.einsum('blkh,bkhc->blhc',att,v)
        o = (torch.sigmoid(self.glinear(z).reshape(L1,L2,self.N_head,self.c)) * o).reshape(L1,L2,-1)
        o = self.olinear(o)
        o = o.permute(1,0,2)
        return o
    
class PairTrans(nn.Module):
    def __init__(self,z_dim,c_expand=2):
        super(PairTrans,self).__init__()
        self.z_dim=z_dim
        self.c_expand=c_expand
        self.norm = nn.LayerNorm(z_dim)
        self.linear1 = Linear(z_dim,z_dim*c_expand)
        self.linear2 = Linear(z_dim*c_expand,z_dim)
    def forward(self,z):
        a = self.linear1(self.norm(z))
        a = self.linear2(F.relu(a))
        return a 

def SignedSqrt( x):
    x = torch.sqrt(torch.relu(x)) - torch.sqrt(torch.relu(-x))
    return x

class MSARow(nn.Module):
    def __init__(self,m_dim,z_dim,N_head=8,c=8):
        super(MSARow,self).__init__()
        self.N_head = N_head
        self.c = c
        self.sq_c = 1/math.sqrt(c)
        self.norm1=nn.LayerNorm(m_dim)
        self.qlinear = LinearNoBias(m_dim,N_head*c)
        self.klinear = LinearNoBias(m_dim,N_head*c)
        self.vlinear = LinearNoBias(m_dim,N_head*c)
        self.norm_z  = nn.LayerNorm(z_dim)
        self.zlinear = LinearNoBias(z_dim,N_head)
        self.glinear = Linear(m_dim,N_head*c)
        self.olinear = Linear(N_head*c,m_dim)
        
        #self.bplinear = Linear(1,N_head)

    def forward(self,m,z,bp):
        N,L,D = m.shape  
        m = self.norm1(m)
        q = self.qlinear(m).reshape(N,L,self.N_head,self.c) #s rq h c 
        k = self.klinear(m).reshape(N,L,self.N_head,self.c) #s rv h c 
        v = self.vlinear(m).reshape(N,L,self.N_head,self.c)
        b = self.zlinear(self.norm_z(z))
        g = torch.sigmoid(self.glinear(m)).reshape(N,L,self.N_head,self.c)
         
        att=torch.einsum('bqhc,bvhc->bqvh',q,k) * (self.sq_c) + b[None,:,:,:] #+ self.bplinear(bp)[None,:,:,:]# rq rv h
        att=F.softmax(SignedSqrt(att),dim=2)
        if attdrop:
            if self.training:
                att = DropAtt(att,dim=2)
        o = torch.einsum('bqvh,bvhc->bqhc',att,v) * g
        m_ = self.olinear(o.reshape(N,L,-1))
        return m_

class MSACol(nn.Module):
    def __init__(self,m_dim,N_head=8,c=8):
        super(MSACol,self).__init__()
        self.N_head = N_head
        self.c = c
        self.sq_c = 1/math.sqrt(c)
        self.norm1=nn.LayerNorm(m_dim)
        self.qlinear = LinearNoBias(m_dim,N_head*c)
        self.klinear = LinearNoBias(m_dim,N_head*c)
        self.vlinear = LinearNoBias(m_dim,N_head*c)

        self.glinear = Linear(m_dim,N_head*c)
        self.olinear = Linear(N_head*c,m_dim)
        
        #self.bplinear = Linear(1,N_head)

    def forward(self,m,bp):
        # m : N L 32  
        N,L,D = m.shape  
        m = self.norm1(m)
        q = self.qlinear(m).reshape(N,L,self.N_head,self.c) #s rq h c 
        k = self.klinear(m).reshape(N,L,self.N_head,self.c) #s rv h c 
        v = self.vlinear(m).reshape(N,L,self.N_head,self.c)

        g = torch.sigmoid(self.glinear(m)).reshape(N,L,self.N_head,self.c)

        att=torch.einsum('slhc,tlhc->stlh',q,k) * (self.sq_c) #+ self.bplinear(bp)[None,:,:,:]  # rq rv h
        att=F.softmax(SignedSqrt(att),dim=1)
        if attdrop:
            if self.training:
                att = DropAtt(att,dim=1)
        o = torch.einsum('stlh,tlhc->slhc',att,v) * g
        m_ = self.olinear(o.reshape(N,L,-1))
        return m_

class MSATrans(nn.Module):
    def __init__(self,m_dim,c_expand=2):
        super(MSATrans,self).__init__()
        self.c_expand=4
        self.m_dim=m_dim
        self.norm=nn.LayerNorm(m_dim)
        self.linear1 = Linear(m_dim,m_dim*c_expand)
        self.linear2 = Linear(m_dim*c_expand,m_dim)
    def forward(self,m):
        m = self.norm(m)
        m = self.linear1(m)
        m = self.linear2(F.relu(m))
        return m

class MSAOPM(nn.Module):
    def __init__(self,m_dim,z_dim,c=12):
        super(MSAOPM,self).__init__()
        self.m_dim=m_dim
        self.c=c
        self.norm=nn.LayerNorm(m_dim)
        self.linear1=Linear(m_dim,c)
        self.linear2=Linear(m_dim,c)
        self.linear3=Linear(c*c,z_dim)
    def forward(self,m):
        N,L,D=m.shape
        o=self.norm(m)
        a=self.linear2(o)
        b=self.linear1(o)
        o = torch.einsum('nia,njb->nijab',a,b).mean(dim=0)
        o = self.linear3(o.reshape(L,L,-1))
        return o

class EvoBlock(nn.Module):
    def __init__(self,m_dim,z_dim,docheck=False):
        super(EvoBlock,self).__init__()
        
        self.msa_row=MSARow(m_dim,z_dim)
        self.msa_col=MSACol(m_dim)
        self.msa_trans=MSATrans(m_dim)

        self.msa_opm=MSAOPM(m_dim,z_dim)

        self.pair_triout= TriOut(z_dim)
        self.pair_triin = TriIn(z_dim)
        self.pair_tristart= TriAttStart(z_dim)
        self.pair_triend  = TriAttEnd(z_dim)
        self.pair_trans = PairTrans(z_dim)
        #self.docheck=docheck
        #if docheck:
            #print('will do checkpoint')

    def layerfunc_msa_row(self,m,z):
        return self.msa_row(m,z) + m
    def layerfunc_msa_col(self,m):
        return self.msa_col(m) + m
    def layerfunc_msa_trans(self,m):
        return self.msa_trans(m) + m
    def layerfunc_msa_opm(self,m,z):
        return self.msa_opm(m) + z

    def layerfunc_pair_triout(self,z):
        return self.pair_triout(z) + z
    def layerfunc_pair_triin(self,z):
        return self.pair_triin(z) + z
    def layerfunc_pair_tristart(self,z):
        return self.pair_tristart(z) + z
    def layerfunc_pair_triend(self,z):
        return self.pair_triend(z) + z      
    def layerfunc_pair_trans(self,z):
        return self.pair_trans(z) + z  
    def forward(self,m,z,bp):
        if True:
            m = m + self.msa_row(m,z,bp)
            m = m + self.msa_col(m,bp)
            m = m + self.msa_trans(m)
            z = z + self.msa_opm(m)
            z = z + self.pair_triout(z)
            z = z + self.pair_triin(z)
            z = z + self.pair_tristart(z, bp)
            z = z + self.pair_triend(z, bp)
            z = z + self.pair_trans(z)
            return m,z
        
def one_d(idx_, d, max_len=2056*8):
    idx = idx_[None]
    K = torch.arange(d//2).to(idx.device)
    sin_e = torch.sin(idx[..., None] * math.pi / (max_len**(2*K[None]/d))).to(idx.device)
    cos_e = torch.cos(idx[..., None] * math.pi / (max_len**(2*K[None]/d))).to(idx.device)
    return torch.cat([sin_e, cos_e], axis=-1)[0]

class Evoformer(nn.Module):
    def __init__(self,m_dim,z_dim,docheck=False):
        super(Evoformer,self).__init__()
        self.layers=[9]
        self.docheck=docheck
        if docheck:
            pass
            #print('will do checkpoint')
        self.evos=nn.ModuleList([EvoBlock(m_dim,z_dim,True) for i in range(self.layers[0])])
        self.olayer = Linear(z_dim,1)
        
        self.poslinear=Linear(64,z_dim)
        self.lm_layer_m = Linear(512,m_dim)
        self.lm_layer_z = Linear(128,z_dim)
        self.lm_layer_bp = Linear(1,z_dim)

    def layerfunc(self,layermodule,m,z):
        m_,z_=layermodule(m,z)
        return m_,z_

    def forward(self, m, z, bp, seq_idx):
        
        N,L,D = m.shape
        
        seq_idx = seq_idx[None]
        relative_pos = seq_idx[:, :, None] - seq_idx[:, None, :]
        relative_pos = relative_pos.reshape([1, L**2])
        relative_pos = one_d(relative_pos,64)

        z = self.lm_layer_z(z) + self.poslinear(relative_pos.reshape([1, L, L, -1])[0]) + self.lm_layer_bp(bp)
        m = self.lm_layer_m(m)
        
        for i in range(len(self.evos)):
            m,z=self.evos[i](m,z,bp)
        out = self.olayer(z)
        
        return out.squeeze(-1)

    #def forward_n(self,m,z,starti,endi):
        #for i in range(starti,endi):
            #print(i)
            #m,z=self.evos[i](m,z)
        #return m,z
    #def forward(self,m,z):
        

        #m,z = checkpoint(self.forward_n,m,z,0,3)
        #m,z = checkpoint(self.forward_n,m,z,3,6)
        #m,z = checkpoint(self.forward_n,m,z,6,10)
        #m,z = checkpoint(self.forward_n,m,z,10,13)
        #m,z = checkpoint(self.forward_n,m,z,13,16)
        #return m,z






















