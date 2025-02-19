#!/bin/bash

#Basepath="/home/user/デスクトップ/data1/refinement/globalRNArefiner"
#Basepath="/Users/tsukiyamashou/Documents/Project_RNA3D/RNARefiner/RNArefiner/bbRNArefiner"
#Basepath="/Users/tsukiyamashou/Documents/Project_RNA3D/RNARefiner/RNArefiner/bbRNArefiner"
Basepath=$(dirname "$(realpath "$0")")

#PYTHON="/home/dl-box/anaconda3/envs/drfold2/bin/python"
#PYTHON="/Users/tsukiyamashou/anaconda3/bin/python"
#PYTHON="/home/dl-box/anaconda3/envs/drfold2/bin/python"

#export PETFOLDBIN="/RNArefinement/DRfold/PETfold/bin"

#if [ ! -d $2/RNArefiner/deep ]; then
    #mkdir -p $2/RNArefiner/deep
    #mkdir -p $2/RNArefiner/template
#fi

PETFOLDDIR=$(awk -F'"' '/"PETfoldPATH":/ {print $4}' $Basepath/config.json)
PETFOLDBIN=$(dirname "$PETFOLDDIR")
export PETFOLDBIN="$PETFOLDBIN"
ArenaDIR=$(awk -F'"' '/"ArenaPATH":/ {print $4}' $Basepath/config.json)
PYTHON=$(awk -F'"' '/"PythonPATH":/ {print $4}' $Basepath/config.json)

INSTR="$1"  #initialstructure.pdb
IN="$2"                # input.fasta
WDIR="$3"/RNArefiner  # working folder

INREAL=$(realpath $IN)

mkdir -p $WDIR/bp/e2e
mkdir -p $WDIR/stack
mkdir -p $WDIR/refinement

if [ ! -s $WDIR/bp/e2e/seq_ss.npy ]
then
    echo "Generating ss features"
    $PYTHON $Basepath/e2e/scripts/Feature.py $INREAL $WDIR/bp/e2e/seq_ss 
fi


if [ ! -s $WDIR/bp/e2e/e2e_5.npy ]
then
    echo "Predicting e2e models"
    $PYTHON $Basepath/e2e/DeepE2EPotential/predict.py $INREAL $WDIR/bp/e2e/seq_ss.npy $WDIR/bp/e2e/e2e
fi

echo "Generating template structures"
for((i=0; i<=5; i++))
do
    if [ ! -s $WDIR/bp/e2e/e2e_$i.pdb ]
    then
        $PYTHON $Basepath/e2e/npy2pdb.py $WDIR/bp/e2e/e2e_$i.npy $INREAL $WDIR/bp/e2e/e2e_$i.pdb
    fi
done

if [ ! -s $WDIR/refinement/initial_structure.pdb ]
then
    $ArenaDIR $INSTR $WDIR/refinement/initial_structure.pdb
fi

if [ ! -s $WDIR/bp/bp.list ]
then
    echo "Predicting base pairs"
    $PYTHON $Basepath/bp/predict.py $WDIR/refinement/initial_structure.pdb $WDIR/bp/e2e/e2e_0.pdb $WDIR/bp/e2e/e2e_1.pdb $WDIR/bp/e2e/e2e_2.pdb $WDIR/bp/e2e/e2e_3.pdb $WDIR/bp/e2e/e2e_4.pdb $WDIR/bp/e2e/e2e_5.pdb $WDIR/bp/bp.list 
fi

if [ ! -s $WDIR/stack/stack.list ]
then
    echo "Predicting stacking pairs"
    $PYTHON $Basepath/stack/predict.py $WDIR/refinement/initial_structure.pdb $WDIR/bp/bp.list  $WDIR/stack/stack
fi

#echo "Start refinement"
mkdir -p $WDIR/refinement/e2e_$model_idx/base
echo "$Basepath/src/baserefinement/baserefinement -i $WDIR/refinement/initial_structure.pdb -o $WDIR/refinement -bp $WDIR/bp/bp.list -st $WDIR/stack/stack.list -time 10"
$Basepath/src/baserefinement/baserefinement -i $WDIR/refinement/initial_structure.pdb -o $WDIR/refinement -bp $WDIR/bp/bp.list -st $WDIR/stack/stack.list -time 10


#if [ ! -s $WDIR/contact/contact_20.list ]
#then
    #echo "Predicting contact map"
    #$PYTHON $Basepath/contact/predict.py $WDIR/contact $WDIR/e2e/e2e_0.pdb $WDIR/e2e/e2e_1.pdb $WDIR/e2e/e2e_2.pdb $WDIR/e2e/e2e_3.pdb $WDIR/e2e/e2e_4.pdb $WDIR/e2e/e2e_5.pdb $WDIR/contact/contact_8.list 8
    #$PYTHON $Basepath/contact/predict.py $WDIR/contact $WDIR/e2e/e2e_0.pdb $WDIR/e2e/e2e_1.pdb $WDIR/e2e/e2e_2.pdb $WDIR/e2e/e2e_3.pdb $WDIR/e2e/e2e_4.pdb $WDIR/e2e/e2e_5.pdb $WDIR/contact/contact_20.list 20
#fi

#model_idx=0
#backbone refinement
#mkdir -p $WDIR/refinement/e2e_$model_idx/bb
#$Basepath/src/bbrefinement/bbrefinement -i $WDIR/e2e/e2e_$model_idx"_bb.pdb" -o $WDIR/refinement/e2e_$model_idx"/bb" -t $WDIR/template -bp $WDIR/bp/bp.list -st $WDIR/stack/stack.list -ct1 $WDIR/contact/contact_8.list -ct2 $WDIR/contact/contact_20.list -r 3 -c2 10

#refinement of backbone + sugar ring + parts of base atoms
#mkdir -p $WDIR/refinement/e2e_$model_idx/sugar
#$Basepath/src/Attachatoms/InitStr -i $WDIR/refinement/e2e_$model_idx"/bb/final_structure.pdb" -o $WDIR/refinement/e2e_$model_idx"/sugar/e2e_"$model_idx"_sugar.pdb" -t 2
#$Basepath/src/sugarrefinement/sugarrefinement -i $WDIR"/refinement/e2e_"$model_idx"/sugar/e2e_"$model_idx"_sugar.pdb" -o $WDIR"/refinement/e2e_"$model_idx"/sugar" -t $WDIR/template -bp $WDIR/bp/bp.list -st $WDIR/stack/stack.list -ct1 $WDIR/contact/contact_8.list -ct2 $WDIR/contact/contact_20.list -r 3 -c2 10

#$Basepath/src/Attachatoms/InitStr -i $WDIR/e2e/e2e_$i.pdb -o $WDIR/e2e/e2e_$model_idx"_bb.pdb" -t 2
#$Basepath/src/sugarrefinement/sugarrefinement -i $WDIR/e2e/e2e_$model_idx"_bb.pdb" -o $WDIR"/refinement/e2e_"$model_idx"/sugar" -t $WDIR/template -bp $WDIR/bp/bp.list -st $WDIR/stack/stack.list -ct1 $WDIR/contact/contact_8.list -ct2 $WDIR/contact/contact_20.list -r 3 -c2 10

#full atomic refinement
#mkdir -p $WDIR/refinement/e2e_$model_idx/base
#$Basepath/src/Attachatoms/InitStr -i $WDIR/refinement/e2e_$model_idx"/sugar/final_structure.pdb" -o $WDIR/refinement/e2e_$model_idx"/base/e2e_$model_idx.pdb" -t 1
#$Basepath/src/baserefinement/baserefinement -i $WDIR/refinement/e2e_$model_idx"/base/e2e_$model_idx.pdb" -o $WDIR/refinement/e2e_$model_idx"/base" -t $WDIR/template -bp $WDIR/bp/bp.list -st $WDIR/stack/stack.list -ct1 $WDIR/contact/contact_8.list -ct2 $WDIR/contact/contact_20.list -r 3 -c2 10

#$Basepath/src/baserefinement/baserefinement -i $WDIR"/refinement/initial_structure.pdb" -o $WDIR/refinement/e2e_$model_idx"/base" -t $WDIR/template -bp $WDIR/bp/bp.list -st $WDIR/stack/stack.list -ct1 $WDIR/contact/contact_8.list -ct2 $WDIR/contact/contact_20.list -r 3 -c2 10



#for((i=0; i<=5; i++))
#do
    #mkdir -p $WDIR/RNArefiner/structure_$i
    #cp -r $WDIR/RNArefiner/template $WDIR/RNArefiner/structure_$i
#$Basepath/src/RNArefiner -i $WDIR/structures/structure_$i.pdb -o $WDIR/RNArefiner/structure_$i -bp $WDIR/bp/bp.list -st $WDIR/stack/stack.list -t 0 -c1 40 -c2 6
#$Basepath/src/RNArefiner -i $WDIR/structures/structure_$i.pdb -o $WDIR/RNArefiner/structure_$i -bp $WDIR/bp/bp.list -st $WDIR/stack/stack.list -t 3 -c2 5
#done

#$Basepath/src/RNArefiner -i $WDIR/structures/structure_0.pdb -o $WDIR/RNArefiner/structure_0 -bp $WDIR/bp/bp.list -st $WDIR/stack/stack.list -ct1 $WDIR/contact/contact_8.list -ct2 $WDIR/contact/contact_20.list -t 0 -c1 10 -c2 5
#$Basepath/src/RNArefiner -i $WDIR/structures/structure_0.pdb -o $WDIR/RNArefiner/structure_0 -bp $WDIR/bp/bp.list -st $WDIR/stack/stack.list -ct1 $WDIR/contact/contact_8.list -ct2 $WDIR/contact/contact_20.list -t 0 -c1 40 -c2 6


