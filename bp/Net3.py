import torch
import torch.nn as nn
from torch.nn import functional as F
from collections import OrderedDict
from EmbSeq import EmbedSeq

def conv3x3(in_planes: int, out_planes: int, stride: int = 1, groups: int = 1, dilation: int = 1) -> nn.Conv2d:
    return nn.Conv2d(in_planes, out_planes, kernel_size=3, stride=stride,
                     padding=dilation, groups=groups, bias=False, dilation=dilation)

def conv1x1(in_planes: int, out_planes: int, stride: int = 1) -> nn.Conv2d:
    return nn.Conv2d(in_planes, out_planes, kernel_size=1, stride=stride, bias=False)

class BaseResBlock(nn.Module):
    def __init__(self, inplanes, planes, stride = 1, dilation = 1):
        super(BaseResBlock, self).__init__()

        self.bn1 = nn.BatchNorm2d(inplanes)
        self.relu1 = nn.ReLU(inplace=True)
        self.conv1 = conv3x3(inplanes, planes, stride)
        self.dropout = nn.Dropout(p=0.3)
        self.relu2 = nn.ReLU(inplace=True)
        self.conv2 = conv3x3(planes, planes, dilation=dilation)

    def forward(self, x):
        identity = x

        out = self.bn1(x)
        out = self.relu1(out)
        out = self.dropout(out)
        out = self.conv1(out)
        out = self.relu2(out)
        out = self.dropout(out)
        out = self.conv2(out)
        out += identity

        return out

class StrNet(nn.Module):
    def __init__(self, dim_seq_fe, dim_hidden=96, maxL=200):
        super(StrNet, self).__init__()
        self.embedpairseq = EmbedSeq(dim_seq_fe, dim_hidden, dim_hidden)
        self.zlinear = nn.ModuleList([nn.Linear(dim_hidden, dim_hidden) for _ in range(6)])
        self.slinear = nn.Linear(dim_hidden, dim_hidden)
        self.szlinear = nn.Linear(dim_hidden, dim_hidden)
        self.relu = nn.ReLU()
        
        self.dropout = nn.Dropout(p=0.3)
        self.pos = self.compute_pos().float()
        self.poslinear = nn.Linear(65,dim_hidden)

    def compute_pos(self,maxL=2000):
        a = torch.arange(maxL)
        b = (a[None,:]-a[:,None]).clamp(-32,32)
        return F.one_hot(b+32,65)

    def forward(self, seq_, R_matrices, T_vector):
        if self.pos.device != seq_.device:
            self.pos = self.pos.to(seq_.device)

        m, s = self.embedpairseq(seq_)
        B,L,C = m.size()
        C_times = int(C/3)

        m = m.reshape(B,L,C_times,3)
        m = m[:,None,:,None,:,:].repeat(1,T_vector.size(1),1,L,1,1)

        z = torch.einsum('blijmn,blijkm->blijkn',R_matrices,m+T_vector.unsqueeze(4).expand(-1,-1,-1,-1,C_times,-1))
        
        _,N,_,_,_,_ = z.size()
        z = z.reshape(B,N,L,L,C)
        
        z = torch.stack([linear(z[:,i,:,:,:]) for i,linear in enumerate(self.zlinear)])
        s = self.slinear(s)

        z = self.dropout(z)
        s = self.dropout(s)
        
        s_z = self.relu(z * s[None,:,:,:,:])
        #s_z = self.relu(s_z)
        s_z = self.szlinear(s_z)
        s_z = s_z.permute(1,2,3,0,4)
        s_z = s_z  + self.poslinear( self.pos[None,:L,:L,None,:])
        s_z = s_z.reshape(B,L,L,N*C)
        #s_z = s_z.mean(-2) + self.poslinear( self.pos[None,:L,:L,None,:])

        return s_z

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
    def __init__(self, dim_seq_fe, dim_hidden=96, maxL=200):
        super().__init__()
        self.maxL = maxL
        self.strnet = StrNet(dim_seq_fe, dim_hidden)
        #self.pairnet = PairLinear(17, dim_hidden)

        self.resnet = ResNet(dim_hidden*6+45, dim_hidden, 40)
        self.outlinear = nn.Linear(dim_hidden, 1)
        #self.outlinear_stack = nn.Linear(dim_hidden, 1)
        #self.outlinear_map = nn.Linear(dim_hidden, 30)
        #self.resnet = ResNet(dim_hidden+17, dim_hidden, 20)
        self.softmax = nn.Softmax(dim=-1)

    def forward(self, seq_, pair_fe, R_matrices, T_vector):
        z = self.strnet(seq_, R_matrices, T_vector) 
        #pz = self.pairnet(pair_fe)

        z = torch.concat([z, pair_fe], dim = -1)
        z = z.permute(0,3,1,2)
        
        out = self.resnet(z)
        out = out.permute(0, 2, 3, 1)
        out = out + out.permute(0, 2, 1, 3)

        out = self.outlinear(out)
        
        return out.squeeze(-1)
        #out_bp = self.outlinear_bp(out)
        #out_stack = self.outlinear_stack(out)
        #out_maps = self.outlinear_map(out)

        #out_maps = out_maps.reshape([out_maps.size(0), out_maps.size(1), out_maps.size(2), 3, 10])
        #out_maps = self.softmax(out_maps)

        #return out_bp.squeeze(-1), out_stack.squeeze(-1), out_maps
        #return out_bp.squeeze(-1), out_stack.squeeze(-1)








