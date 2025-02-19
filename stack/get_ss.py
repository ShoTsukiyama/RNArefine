
import sys
import os
import subprocess
import numpy as np
import pandas as pd
import torch

PythonPATH="/home/dl-box/anaconda3/envs/drfold2/bin/python"
RNAfoldPATH="/usr/local/bin/RNAfold"
CSSRPATH="/RNArefinement/program/src/single/bp_stack_contact/CSSR-Linux64"
CT2DBPATH="/usr/local/bin/ct2db"
PETfoldPATH="/RNArefinement/DRfold/PETfold/bin/PETfold"
LinearfoldPATH="/RNArefinement/program/src/single/bp_stack_contact/LinearFold-master/linearfold"

def load_txt(filepath):
    with open(filepath) as f:
        data = f.readlines()
    return data

def out_txt(filepath, data):
    with open(filepath, "w") as f:
        f.write(data)

def exec_command(command):
    results = subprocess.run(command, capture_output=True, text=True)

    return results.stdout, results.stderr 

def exec_command_pipe(command1, command2):
    try:
        with subprocess.Popen(command1, stdout=subprocess.PIPE, text=True) as proc1:
            with subprocess.Popen(command2, stdin=proc1.stdout, stdout=subprocess.PIPE, text=True) as proc2:
                output = proc2.communicate()[0]
                
        return output
    except subprocess.CalledProcessError as e:
        print(f"Command failed")
        print("Output:", e.output)
        print("Error:", e.stderr)
        sys.exit()

def exec_command_stdin(command_list, content):
    process = subprocess.Popen(
                command_list,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
                )
    stdout, stderr = process.communicate(input=content)

    return stdout, stderr

def get_ss_info(brk, seqnum, inv = True):
    pairs, pair_types = [], []
    brk1, brk2, brk3, brk4 = 0, 0, 0, 0
    brk1_p, brk2_p, brk3_p, brk4_p = [], [], [], []
    pairs_mat = np.zeros([seqnum, seqnum])

    for i in range(len(brk)):
        if(brk[i] == "("):
            brk1 += 1
            brk1_p.append(i)

        elif(brk[i] == "["):
            brk2 += 1
            brk2_p.append(i)
        elif(brk[i] == "<"):
            brk3 += 1
            brk3_p.append(i)
        elif(brk[i] == "{"):
            brk4 += 1
            brk4_p.append(i)
        if(brk[i] == ")"):
            brk1 -= 1
            pairs.append((brk1_p[-1], i))
            pair_types.append(0)
            if(inv): pairs.append((i, brk1_p[-1]))
            pairs_mat[brk1_p[-1], i] = 1
            pairs_mat[i, brk1_p[-1]] = 1
            brk1_p.pop()
        elif(brk[i] == "]"):
            brk2 -= 1
            pairs.append((brk2_p[-1], i))
            pair_types.append(1)
            if(inv): pairs.append((i, brk2_p[-1]))
            pairs_mat[brk2_p[-1], i] = 1
            pairs_mat[i, brk2_p[-1]] = 1
            brk2_p.pop()
        elif(brk[i] == ">"):
            brk3 -= 1
            pairs.append((brk3_p[-1], i))
            pair_types.append(2)
            if(inv): pairs.append((i, brk3_p[-1]))
            pairs_mat[brk3_p[-1], i] = 1
            pairs_mat[i, brk3_p[-1]] = 1
            brk3_p.pop()
        elif(brk[i] == "}"):
            brk4 -= 1
            pairs.append((brk4_p[-1], i))
            pair_types.append(3)
            if(inv): pairs.append((i, brk4_p[-1]))
            pairs_mat[brk4_p[-1], i] = 1
            pairs_mat[i, brk4_p[-1]] = 1
            brk4_p.pop()

    return pairs, pair_types, pairs_mat

def get_ct(pairs, seq):
    pairs_ = np.array(pairs).T
    if(len(pairs_) > 0):
        pairs_dict = {**dict(zip(pairs_[0], pairs_[1])), **dict(zip(pairs_[1], pairs_[0]))}
    else:
        pairs_dict = {}
    seqlen = len(seq)
    seqidx = np.arange(1,seqlen+1,1).tolist()
    seqnextidx = [0] + seqidx[:seqlen-1]
    seqprevidx = seqidx[1:] + [0]

    ct=["{}".format(seqlen)]
    for i in range(seqlen):
        if(i in pairs_dict):
            ct.append("{} {}     {}     {}     {}  {}".format(seqidx[i], seq[i], seqnextidx[i], seqprevidx[i], pairs_dict[i]+1,seqidx[i]))
        else:
            ct.append("{} {}     {}     {}     {}  {}".format(seqidx[i], seq[i], seqnextidx[i], seqprevidx[i], 0,seqidx[i]))

    return "\n".join(ct)

def get_ss_linearfold(basepath, seq, model_type = None):
    command1 = ["echo", seq]
    if(model_type is None):
        command2 = [LinearfoldPATH]
    else:
        command2 = [LinearfoldPATH, "-V"]

    res = exec_command_pipe(command1, command2)

    seq_len = len(seq)
    res = res.split("\n")
    if(len(res) >= 2):
        if(len(res[1]) >= seq_len):
            ss = res[1][:seq_len]
        else:
            print(res)
            sys.exit()
    else:
        print(res)
        sys.exit()

    return ss

def get_ss_linearfold_with_str(basepath, seq, pairs, pair_types, model_type = None):
    seq_len = len(seq)

    pseudopairs = []
    constraints = list("?"*len(seq))
    for i in range(len(pairs)):
        if(pair_types[i] == 0):
            if((seq[pairs[i][0]] == "A" and seq[pairs[i][1]] == "U") or
               (seq[pairs[i][0]] == "U" and seq[pairs[i][1]] == "A") or
               (seq[pairs[i][0]] == "G" and seq[pairs[i][1]] == "C") or
               (seq[pairs[i][0]] == "C" and seq[pairs[i][1]] == "G") or
               (seq[pairs[i][0]] == "G" and seq[pairs[i][1]] == "U") or
               (seq[pairs[i][0]] == "U" and seq[pairs[i][1]] == "G")):
                if(pairs[i][0] < pairs[i][1]):
                    constraints[pairs[i][0]] = "("
                    constraints[pairs[i][1]] = ")"
                else:
                    constraints[pairs[i][0]] = ")"
                    constraints[pairs[i][1]] = "("
        else:
            pseudopairs.append((pairs[i][0], pairs[i][1]))
            constraints[pairs[i][0]] = "."
            constraints[pairs[i][1]] = "."

    constraints = "".join(constraints)

    command1 = ["echo", "{}\n{}".format(seq, constraints)]

    if(model_type is None):
        command2 = [LinearfoldPATH, "--constraints"]
    else:
        command2 = [LinearfoldPATH, "-V", "--constraints"]

    res = exec_command_pipe(command1, command2)
    res = res.split("\n")
    if(len(res) >= 3):
        if(len(res[2]) >= seq_len):
            new_pairs, _, _ = get_ss_info(res[2][:seq_len], len(seq), False)
        else:
            print(res)
            sys.exit()
    else:
        print(res)
        sys.exit()

    pairs = np.array(new_pairs + pseudopairs)
    ct = get_ct(pairs, seq)
    out_txt("{}/pairs.ct".format(basepath), ct)
    command = [CT2DBPATH, "{}/pairs.ct".format(basepath)]
    stdout, stderr = exec_command(command)
    res = stdout.split()

    if(len(res) >= 3):
        if(len(res[2]) >= seq_len):
            ss = res[2][:seq_len]
        else:
            print(stdout)
            print(stderr)
            sys.exit()
    else:
        print(stdout)
        print(stderr)
        sys.exit()

    return ss

def get_ss_petfold(basepath, seq):
    fasta_ = "{}/sequence.fa".format(basepath)
    ppfile          ="{}/pairs.petfoldrr".format(basepath)

    command = [PETfoldPATH, "-f", fasta_, "-r", ppfile]
    stdout, stderr = exec_command(command)
    
    pfold, petfold = None, None
    if("Pfold RNA structure:" in stdout):
        pfold = stdout.split("\n")[0].split("Pfold RNA structure:")[1].strip()

    if("PETfold RNA structure:" in stdout):
        petfold = stdout.split("\n")[2].split("PETfold RNA structure:")[1].strip()

    return pfold, petfold

def get_ss_rnafold(seq):
    seq_len = len(seq)
    command = [RNAfoldPATH]

    stdout, stderr = exec_command_stdin(command, ">seq\n{}".format(seq))
    res = stdout.split("\n")
    if(len(res) >= 3):
        if(len(res[2]) >= seq_len):
            ss = res[2][:seq_len]
            return ss
        else:
            print(stdout)
            print(stderr)
            sys.exit()
    else:
        print(stdout)
        print(stderr)
        sys.exit()

    return ss

def get_ss_rnafold_with_str(basepath, seq, pairs, pair_types):
    seq_len = len(seq)
    constraints = list("."*seq_len)
    pseudopairs = []

    for i in range(len(pairs)):
        if(pair_types[i] == 0):
            #constraints.append("F {} {} 1".format(pairs[i][0]+1, pairs[i][1]+1))
            if(pairs[i][0] < pairs[i][1]):
                constraints[pairs[i][0]] = "("
                constraints[pairs[i][1]] = ")"
            else:
                constraints[pairs[i][0]] = ")"
                constraints[pairs[i][1]] = "("
        else:
            #constraints.append("P {} {} 1".format(pairs[i][0]+1, pairs[i][1]+1))
            constraints[pairs[i][0]] = "x"
            constraints[pairs[i][1]] = "x"
            pseudopairs.append((pairs[i][0], pairs[i][1]))

    constraints="\n".join(constraints)
    #out_txt("{}/constraints_rnafold.txt".format(basepath), constraints)

    #command = [RNAfoldPATH, "--commands={}/constraints_rnafold.txt".format(basepath)]
    command = [RNAfoldPATH, "--enforceConstraint", "-C"]
    stdout, stderr = exec_command_stdin(command, ">seq\n{}\n{}".format(seq, constraints))
    res = stdout.split("\n")
    if(len(res) >= 3):
        if(len(res[2]) >= seq_len):
            new_pairs, _, _ = get_ss_info(res[2][:seq_len], len(seq), False)
        else:
            print(stdout)
            print(stderr)
            sys.exit()
    else:
        print(stdout)
        print(stderr)
        sys.exit()

    pairs = np.array(new_pairs + pseudopairs)
    ct = get_ct(pairs, seq)
    out_txt("{}/pairs.ct".format(basepath), ct)
    command = [CT2DBPATH, "{}/pairs.ct".format(basepath)]
    stdout, stderr = exec_command(command)
    res = stdout.split()

    if(len(res) >= 3):
        if(len(res[2]) >= seq_len):
            ss = res[2][:seq_len]
            return ss
        else:
            print(stdout)
            print(stderr)
            sys.exit()
    else:
        print(stdout)
        print(stderr)
        sys.exit()

    print(stdout)
    print(stderr)
    sys.exit()

def get_cssr(pdbfile):
    command = [CSSRPATH, pdbfile, "-", "--outfmt", "1"]
    result = subprocess.run(command, capture_output=True, text=True)
    brk = result.stdout

    return brk

def main_get_seq_ss(basepath, seq):
    out_txt("{}/sequence.fa".format(basepath), ">seq\n{}".format(seq))

    ss_linearfoldv = get_ss_linearfold(basepath, seq, "V")
    ss_linearfold = get_ss_linearfold(basepath, seq)
    _, ss_petfold = get_ss_petfold(basepath, seq)
    ss_rnafold = get_ss_rnafold(seq)

    _, _, pairs_mat_linearfoldv = get_ss_info(ss_linearfoldv, len(seq), False)
    _, _, pairs_mat_linearfold = get_ss_info(ss_linearfold, len(seq), False)
    _, _, pairs_mat_petfold = get_ss_info(ss_petfold, len(seq), False)
    _, _, pairs_mat_rnafold = get_ss_info(ss_rnafold, len(seq), False)

    return [torch.tensor(pairs_mat_linearfoldv), torch.tensor(pairs_mat_linearfold), torch.tensor(pairs_mat_petfold), torch.tensor(pairs_mat_rnafold)]

def main_get_struct_ss(basepath, seq, pdbfile):
    out_txt("{}/sequence.fa".format(basepath), ">seq\n{}".format(seq))

    brk = get_cssr(pdbfile)
    pairs, pair_types, pairs_mat_cssr = get_ss_info(brk, len(seq), False)

    ss_linearfoldv_str = get_ss_linearfold_with_str(basepath, seq, pairs, pair_types, "-V")
    ss_linearfold_str = get_ss_linearfold_with_str(basepath, seq, pairs, pair_types)
    ss_rnafold_str = get_ss_rnafold_with_str(basepath, seq, pairs, pair_types)

    _, _, pairs_mat_linearfoldv_str = get_ss_info(ss_linearfoldv_str, len(seq), False)
    _, _, pairs_mat_linearfold_str = get_ss_info(ss_linearfold_str, len(seq), False)
    _, _, pairs_mat_rnafold_str = get_ss_info(ss_rnafold_str, len(seq), False)

    return [torch.tensor(pairs_mat_cssr), torch.tensor(pairs_mat_linearfoldv_str), torch.tensor(pairs_mat_linearfold_str), torch.tensor(pairs_mat_rnafold_str)]




















