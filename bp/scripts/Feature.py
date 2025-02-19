import os,sys
import numpy as np 
from pathlib import Path
import json

def ss2matrix(ssstr):
    L=len(ssstr)
    ssmatrix = np.zeros( [L,L])
    slist=[]
    for i in range(L):
        if ssstr[i] in ['(','<','[','{']:
            slist.append(i)
        elif ssstr[i] in [')','>',']','}']:
            j = slist[-1]
            ssmatrix[i,j]=1
            ssmatrix[j,i]=1
            slist.pop(-1)
        elif ssstr[i] not in  ['.','-']:
            print('unknown ss state',ssstr[i], i)
            assert False
    return ssmatrix
def PETfold_runner(fastafile,saveprefix,petfoldbin):
    #petfoldbin =f'/nfs/amino-home/liyangum/eRNA/programs/PETfold/bin/PETfold'
    petfold_savefile=saveprefix+'.petfold'
    ppfile          =saveprefix+'.petfoldrr'
    fasta_ = os.path.abspath(fastafile)
    cmd = f'{petfoldbin} -f {fasta_} -r {ppfile} >{petfold_savefile}'
    os.system(cmd)
    petfoldline=open(petfold_savefile).readlines()[-2].strip()
    ss1=ss2matrix(petfoldline.split()[-1])
    ss2=np.genfromtxt(ppfile,skip_header=1,skip_footer=1)
    seq = open(fastafile).readlines()[1].strip()
    L = len(seq)
    L1,L2 = ss2.shape
    L3 = ss1.shape[0]
    os.remove(petfold_savefile)
    os.remove(ppfile)
    if L==L1 and L==L2 and L==L3:
        #np.save(saveprefix,np.stack(ss1,ss2))
        return np.concatenate([ss1[...,None],ss2[...,None]],axis=-1)
    else:
        print('PETfold lengths does not id',L,L1,L2,L3)
        assert False


def ViennaRNA_runner(fastafile,saveprefix_,ViennaRNAbin):
    #ViennaRNAbin='/nfs/amino-home/liyangum/eRNA/programs/ViennaRNA-2.4.18-bin/bin/RNAfold'
    saveprefix = os.path.abspath(saveprefix_)
    path = Path(saveprefix)
    seq=open(fastafile).readlines()[1].strip()
    workingdir = path.parent.absolute()
    name = os.path.basename(saveprefix)
    newfastafile = saveprefix+'.tmp.fasta'
    lines=[f'>{name}_ViennaRNA_name',seq]
    wfile =open(newfastafile,'w')
    wfile.write('\n'.join(lines))
    wfile.close()
    os.chdir(workingdir)
    cmd = f'{ViennaRNAbin} --outfile={name}.ViennaRNAss -p  --noPS {os.path.abspath(newfastafile)}'
    os.system(cmd)
    ss_str = open(f'{saveprefix}.ViennaRNAss').readlines()[2].strip().split()[0]
    ss1=ss2matrix(ss_str)
    L=len(seq)
    pslines=open(os.path.join(workingdir,f'{name}_ViennaRNA_name_dp.ps')).readlines()
    pslines=[aline.strip() for aline in pslines]
    pslines=[aline.strip() for aline in pslines if aline.endswith('ubox') and len(aline.split())==4 and aline[0].isdigit()]

    ss2 = np.zeros([L,L])
    for aline in pslines:
        words = aline.split()
        i,j,score = int(words[0])-1,int(words[1])-1,float(words[2])
        ss2[i,j]=score
        ss2[j,i]=score
    os.remove(os.path.join(workingdir,f'{name}_ViennaRNA_name_dp.ps'))
    os.remove(newfastafile)
    os.remove(f'{saveprefix}.ViennaRNAss')
    if L==ss1.shape[0]:
        return np.concatenate([ss1[...,None],ss2[...,None]],axis=-1)
    else:
        print('ViennaRNA lengths does not id',L,ss1.sape)
        assert False


def Extraxt_ss(fastafile,saveprefix,petfoldbin,ViennaRNAbin):
    savepath = os.path.abspath(saveprefix) +''
    ss1 = PETfold_runner(fastafile,saveprefix+'_pet',petfoldbin)
    ss2 = ViennaRNA_runner(fastafile,saveprefix+'_vie',ViennaRNAbin)
    ss = np.concatenate([ss1,ss2],axis=-1)
    np.save(savepath,ss)
    return ss

def load_json(filepath):
    with open(filepath, 'r') as file:
        data = json.load(file)
    return data    

#if __name__ == '__main__':
def predict_ss(fastafile, saveprefix):
    #fastafile=sys.argv[1]
    #saveprefix=sys.argv[2]

    current_directory = os.path.dirname(os.path.abspath(__file__))
    config = load_json(f"{current_directory}/../../config.json")
    petfoldpath = config["PETfoldPATH"]
    rnafoldpath = config["RNAfoldPATH"]
    ss = Extraxt_ss(fastafile,saveprefix,petfoldpath,rnafoldpath)
    return ss
