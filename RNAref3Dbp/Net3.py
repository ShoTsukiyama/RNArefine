import sys
import torch
import torch.nn as nn
from torch.nn import functional as F
from collections import OrderedDict

def conv3x3(in_planes: int, out_planes: int, stride: int = 1, groups: int = 1, dilation: int = 1) -> nn.Conv2d:
    return nn.Conv2d(in_planes, out_planes, kernel_size=3, stride=stride,
                     padding=dilation, groups=groups, bias=False, dilation=dilation)

def conv1x1(in_planes: int, out_planes: int, stride: int = 1) -> nn.Conv2d:
    return nn.Conv2d(in_planes, out_planes, kernel_size=1, stride=stride, bias=False)

def pick_groups(C, max_groups=32):
    for g in [32, 16, 8, 4, 2, 1]:
        if g <= max_groups and C % g == 0:
            return g
    return 1

class BaseResBlock(nn.Module):
    def __init__(self, inplanes, planes, stride = 1, dilation = 1):
        super(BaseResBlock, self).__init__()
        
        g1 = pick_groups(inplanes, 32)
        self.norm1 = nn.GroupNorm(g1, inplanes)

        self.bn1 = nn.BatchNorm2d(inplanes)
        self.relu1 = nn.ReLU(inplace=False)
        self.conv1 = conv3x3(inplanes, planes, stride)
        self.dropout = nn.Dropout(p=0.2)
        
        g2 = pick_groups(planes, 32)
        self.norm2 = nn.GroupNorm(g2, planes)

        self.relu2 = nn.ReLU(inplace=False)
        self.conv2 = conv3x3(planes, planes, dilation=dilation)

    def forward(self, x):
        identity = x

        #out = self.bn1(x)
        out = self.relu1(self.norm1(x))
        #out = self.relu1(out)
        out = self.dropout(out)
        out = self.conv1(out)
        out = self.relu2(self.norm2(out))
        #out = self.relu2(out)
        out = self.dropout(out)
        out = self.conv2(out)
        out = out + identity

        return out

"""
class StrNet(nn.Module):
    def __init__(self, dim_seq_fe, dim_hidden=96, maxL=1200):
        super(StrNet, self).__init__()
        #self.zlinear1 = nn.Linear(147, dim_hidden)
        self.zlinear1 = nn.Linear(19, dim_hidden)
        self.zlinear2 = nn.Linear(dim_hidden, dim_hidden)
        self.relu = nn.ReLU()
        
        #self.dropout = nn.Dropout(p=0.3)
        self.pos = self.compute_pos().float()
        self.poslinear = nn.Linear(65,dim_hidden)

    def compute_pos(self,maxL=2000):
        a = torch.arange(maxL)
        b = (a[None,:]-a[:,None]).clamp(-32,32)
        return F.one_hot(b+32,65)

    def forward(self, init_bp_, lmz, pfe, outerprot):
        if self.pos.device != init_bp_.device:
            self.pos = self.pos.to(init_bp_.device)

        #z = torch.concat([init_bp_, lmz, pfe, outerprot], dim = 3)
        z = torch.concat([init_bp_, pfe, outerprot], dim = 3)
        
        B,L,L,C = z.size()
        
        s_z = self.zlinear1(z)
        s_z = s_z  + self.poslinear( self.pos[None,:L,:L,:])
        s_z = self.relu(s_z)
        s_z = self.zlinear2(s_z)

        return s_z
"""

class RelPos(nn.Module):
    def __init__(self, dim_hidden=96, pos_range=32, maxL=2048):
        super().__init__()
        self.pos_emb = nn.Embedding(2*pos_range + 1, dim_hidden)
        self.pos_range = pos_range
        self.register_buffer("idx", torch.arange(maxL), persistent=False)
        self.scale = nn.Parameter(torch.tensor(1e-2))

    def forward(self, L):
        a = self.idx[:L]
        rel = (a[None, :] - a[:, None]).clamp(-self.pos_range, self.pos_range) + self.pos_range
        return self.scale * self.pos_emb(rel)

class StrNet(nn.Module):
    def __init__(self, dim_hidden=96, maxL=1200):
        super(StrNet, self).__init__()
        #self.zlinear1 = nn.Linear(147, dim_hidden)
        #self.zlinear2 = nn.Linear(dim_hidden, dim_hidden)
        
        self.enc_init  = nn.Linear(1,   dim_hidden)
        #self.enc_lmz   = nn.Linear(128, dim_hidden)
        #self.enc_pfe   = nn.Linear(1,   dim_hidden)
        self.enc_outer = nn.Linear(10,  dim_hidden)
        self.enc_distmap = nn.Linear(1, dim_hidden)
        #self.enc_distmap = nn.Embedding(1, dim_hidden)
        
        self.act   = nn.ReLU(inplace=False)
        #self.drop  = nn.Dropout(p=0.1)

        self.mix1  = nn.Linear(dim_hidden, dim_hidden)
        self.mix2  = nn.Linear(dim_hidden, dim_hidden)

        self.pos   = RelPos(dim_hidden=dim_hidden, pos_range=32)

    def compute_pos(self,maxL=1032):
        a = torch.arange(maxL)
        b = (a[None,:]-a[:,None]).clamp(-32,32)
        return F.one_hot(b+32,65)

    def forward(self, init_bp_, dist_map, pfe, outerprot):
        #if self.pos.device != init_bp_.device:
            #self.pos = self.pos.to(init_bp_.device)

        #z = torch.concat([init_bp_, lmz, pfe, outerprot], dim = 3)
        
        z = (
            self.enc_init(init_bp_.float()) +
            #self.enc_lmz(lmz.float()) +
            #self.enc_pfe(pfe.float()) +
            self.enc_outer(outerprot.float()) + 
            self.enc_distmap(dist_map.float())
        ) 
        
        B,L,_,H = z.shape
        z = z + self.pos(L)[None, :, :, :]  # 位置加算

        z = self.act(self.mix1(z))
        #z = self.drop(z)
        z = self.mix2(z)

        return z


class ResNet(nn.Module):
    def __init__(self, init_dim_hidden, dim_hidden=96, layer_num=5, maxL=200):
        super().__init__() 
        first_layer = nn.Sequential(
            nn.Conv2d(init_dim_hidden, dim_hidden, kernel_size=1),
        )
        res_layers = []
        for i in range(layer_num):
            dilation = pow(2, (i % 3))
            res_layers.append(BaseResBlock(inplanes=dim_hidden, planes=dim_hidden, dilation=dilation))
        res_layers = nn.Sequential(*res_layers)
        final_layer = nn.Conv2d(dim_hidden, 1, kernel_size=1) 

        layers = OrderedDict()
        layers["first"] = first_layer
        layers["resnet"] = res_layers
        #layers["final"] = final_layer

        self.proj = nn.Sequential(layers)

    def forward(self, z):
        output = self.proj(z)
        return output

class PairLinear(nn.Module):
     def __init__(self, dim_init, dim_hidden=96):
         super().__init__()
         self.linear1 = nn.Linear(dim_init, dim_hidden)
         self.linear2 = nn.Linear(dim_hidden, dim_hidden)
         self.relu = nn.ReLU()

     def forward(self, z):
         z = self.linear1(z)
         z = self.relu(z)
         z = self.linear2(z)

         return z

class Net(nn.Module):
    def __init__(self, dim_hidden=96, maxL=1024):
        super().__init__()
        self.maxL = maxL
        self.strnet = StrNet(dim_hidden)
        #self.pairnet = PairLinear(17, dim_hidden)

        self.resnet = ResNet(dim_hidden, dim_hidden, 32)
        self.outlinear = nn.Linear(dim_hidden, int(dim_hidden/2))
        self.outlinear2 = nn.Linear(int(dim_hidden/2), 1)
        #self.outlinear_stack = nn.Linear(dim_hidden, 1)
        #self.outlinear_map = nn.Linear(dim_hidden, 30)
        #self.resnet = ResNet(dim_hidden+17, dim_hidden, 20)
        #self.softmax = nn.Softmax(dim=-1)

    def forward(self, init_bp_, dist_map, pfe, outerprot):
        z = self.strnet(init_bp_, dist_map, pfe, outerprot) 
        
        z = z.permute(0, 3, 1, 2)
        
        out = self.resnet(z)
        out = out.permute(0, 2, 3, 1)
        
        #out = out + out.permute(0, 2, 1, 3)

        out = self.outlinear(out)
        out = self.outlinear2(out).squeeze(-1)
        
        out = 0.5 * (out + out.permute(0, 2, 1))

        return out
        #out_bp = self.outlinear_bp(out)
        #out_stack = self.outlinear_stack(out)
        #out_maps = self.outlinear_map(out)

        #out_maps = out_maps.reshape([out_maps.size(0), out_maps.size(1), out_maps.size(2), 3, 10])
        #out_maps = self.softmax(out_maps)

        #return out_bp.squeeze(-1), out_stack.squeeze(-1), out_maps
        #return out_bp.squeeze(-1), out_stack.squeeze(-1)








