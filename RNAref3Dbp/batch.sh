#!/bin/bash
#base="/home/RNArefiner/program/get_params"
#WDIR="/home/RNArefiner/dataset"
base="/Users/sho/Documents/project/RNArefiner/RNArefiner/stacking_prediction/RNAref3Dbp2"
Model="/Users/sho/Documents/project/RNArefiner/RNArefiner/stacking_prediction/RNAref3Dbp2/models"
WDIR="/Users/sho/Documents/project/RNArefiner/RNArefiner/stacking_prediction/test_new/test_data"

#base="/Users/sho/Documents/project/RNArefiner/RNArefiner/stacking_prediction/program_new3/bp"
#Model="/Users/sho/Documents/project/RNArefiner/RNArefiner/stacking_prediction/program_new3/bp/models"
#INDIR="/Users/sho/Documents/project/RNArefiner/RNArefiner/dataset/train_pdb_data"
#WDIR="/Users/sho/Documents/project/RNArefiner/RNArefiner/stacking_prediction/dataset/train_data"

PythonBase="python"

#for f in $INDIR/*; do
for f in $WDIR/*; do
    file=$(basename "$f")
    echo $file
    
    $PythonBase $base/main.py --pdb $WDIR/$file/pred_Rhofold/unrelaxed_model.pdb --model $Model -o $WDIR/$file/pred_Rhofold
    $PythonBase $base/main.py --pdb $WDIR/$file/pred_DRfold2/unrelaxed_model.pdb --model $Model -o $WDIR/$file/pred_DRfold2
    $PythonBase $base/main.py --pdb $WDIR/$file/pred_DRfold/DPRcg.pdb --model $Model -o $WDIR/$file/pred_DRfold
    $PythonBase $base/main.py --pdb $WDIR/$file/native/native.pdb --model $Model -o $WDIR/$file/native
    
    #$PythonBase $base/main.py --pdb $WDIR/$file/native/native.pdb --model $Model -o $WDIR/$file/native
    #mkdir -p $WDIR/$file
    #cp $INDIR/$file/pred.pdb $WDIR/$file
    #cp $INDIR/$file/seq.fa $WDIR/$file
    #cp $INDIR/$file/clarna_base_pair_info_pred.csv $WDIR/$file
    #cp $INDIR/$file/mc_base_pair_info_pred.csv $WDIR/$file
    #cp $INDIR/$file/cssr_base_pair_info_pred.csv $WDIR/$file
    #cp $INDIR/$file/rnaview_base_pair_info_pred.csv $WDIR/$file
    #cp $INDIR/$file/dssr_base_pair_info_pred.csv $WDIR/$file

    #if [ ! -f $WDIR/$file/RNAref3D_bp_bb.csv ]; then
        #echo $file
        #$PythonBase $base/bppredictor.py --pdb $INDIR/$file/pred.pdb --model $Model -o $WDIR/$file 
    #fi
done
































































