#!/bin/bash

Basepath=$(dirname "$(realpath "$0")")
PETFOLDDIR=$(awk -F'"' '/"PETfoldPATH":/ {print $4}' $Basepath/config.json)
PETFOLDBIN=$(dirname "$PETFOLDDIR")
export PETFOLDBIN="$PETFOLDBIN"

INSTR="$1"  #initialstructure.pdb
WDIR="$2"  # working folder

PYTHONPATH=$(awk -F'"' '/"PythonPATH":/ {print $4}' $Basepath/config.json)

#echo "Predicting base-pair interactions..."
#$PYTHONPATH $Basepath/bppredictor.py --pdb $INSTR -o $WDIR

#echo "Predicting base-stack interactions..."
#$PYTHONPATH $Basepath/stackpredictor.py --pdb $INSTR --bp $WDIR/rnaref3d.bpseq -o $WDIR

#$Basepath/bin/RNAref3D -i $INSTR -o $WDIR/Prerefined_superimposed_str.pdb -bp $WDIR/bp.list -st $WDIR/stack.list -refmode 2

#$Basepath/bin/RNAref3D -i $INSTR -o $WDIR/EM_str_1_temp.pdb -bp $WDIR/bp.list -st $WDIR/stack.list -refmode 4 

#python $Basepath/RNAref3Dbp/main.py -p $INSTR -o $WDIR/RNAref3D
$Basepath/bin/RNAref3D -i $INSTR -o $WDIR/RNAref3D -bp $WDIR/RNAref3D/RNAref3D_bp.txt -st $WDIR/RNAref3D/RNAref3D_stack.txt
#$Basepath/bin/RNAref3D -i $INSTR -o $WDIR/Pre_ENE_1.pdb -bp $WDIR/bp.list -st $WDIR/stack.list -refmode 4


