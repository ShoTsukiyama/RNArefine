<img width="1088" alt="Image" src="https://github.com/user-attachments/assets/b9ce2f46-959a-408e-85bf-170babdf951a" />
This package refines RNA 3D structures and predicts base-pairing and base-stacking interactions.  
RNAoptimizer can be used under The MIT License.  
(This project was accomplished through a collaboration between the Yang Zhang lab in the National University of Singapore and the Hiroyuki Kurata lab in the Kyushu Institute of Technology.)  
  
# 1. Installation
## 1-1. Setting of third-party software
RNAoptimizer used the following third-party software.  
Install the software in your environment and configure the path in the RNArefiner/config.json file.  
| Software  | Description |
| :------------- | :------------- |
| [PETfold](https://rth.dk/resources/petfold/download.php)  | RNA secondary structure prediction tool |
| [RNAfold](https://www.tbi.univie.ac.at/RNA/)  | RNA secondary structure prediction tool  |
| [Linearfold](https://github.com/LinearFold/LinearFold)  | RNA secondary structure prediction tool |
| [CSSR](https://github.com/pylelab/CSSR)  | Base-pair identification tool  |
| [CT2DB](https://www.tbi.univie.ac.at/RNA/ct2db.1.html)  | Tool to produce dot-bracket notation of an RNA secondary structure  |
| [Arena](https://github.com/pylelab/Arena)  | Full-atomic reconstruction of RNA 3d structures  |
  
*Set the PETfold path as follows.
```
export PETFOLDBIN="path/to/PETfold/bin"
```
*In addition to those packages, configure the Python path in the RNArefiner/config.json file. 
  
## 1-2. Downloading of pre-defined parameters
(Not necessary for base-paring and base-stacking prediction)
Download pre-defined parameters from [here](https://1drv.ms/f/c/7340de942a07a197/EkuUvze_mWpJqsx22H2cz1wBSbRjelLEU5s3J6M0EC_mVg?e=d1WekU), and locate them on the /RNArefiner/param.  
  
## 1-3. Setup of refinement package
(Not necessary for base-paring and base-stacking prediction)
Compile on the RNArefiner/src
```
cd RNArefiner/src
make
```

# 2. RNA refinement
This package offers two methods for structure refinement and the prediction of base-pairing and base-stacking interactions: an automated process and a customizable process.  
The automated process enables structure refinement without detailed parameter settings,  
while the customizable process allows for user-defined adjustments to each process.  

## 2-1. Automated process
RNAoptimizer refines the input structure based on predicted base-pairing and base-stacking interactions.  
If you want to run the whole process without any customizations, run the following bash script.
```
bash rnarefiner.sh [Path to the initial structure] [Path to the output directory]
```
  
## 2-2. Customizable process
If you want to run a specific process with any customized processes rather than the whole process of RNAoptimizer, see the following sections.

  
### 2-1. Base-pairing prediction
RNAoptimizer predicts an optimized RNA secondary structure by integrating the input 3D structure with a secondary structure predicted by a ResNet-based neural network as follows.
```
python bp_predict.py --pdb or --fasta [Path to pdb or fasta file] -o [Output directory] [option]
```
(Other options)
| Option  | Description | Default |
| :-------------: | :-------------: | :-------------: |
| -w  | Weight of the secondary structure predicted by the ResNet-based neural network. (see the original paper) | 0.4 |
  
ex)
If you want to obtain the optimized secondary structure of the input structure, run the following command.  
```
python RNArefiner/bin/bppredictor.py --pdb input.pdb  -o ./output
```
  
If you only want to obtain the secondary structure predicted by the ResNet-based neural network from the sequence, run the following command.  
```
python RNArefiner/bin/bppredictor.py --fasta input.fa  -o ./output
```
As a result, bp.list file is outputted under the directory specified by -o.
  
### 2-2. Base-stacking prediction
RNAoptimizer predicts optimized RNA base-stacking interactions by the ResNet-based neural network from the sequence, input structure, and secondary structure.
```
python RNArefiner/bin/bppredictor.py --pdb [pdb file] --bp [bp.list file] -o [Output directory]
```
ex)
```
python RNArefiner/bin/bppredictor.py --pdb input.pdb --bp bp.list -o ./output
```
As a result, stack.list file is outputted under the directory specified by -o.
  
### 2-3. RNA 3d structure refinement
The following command executes refinements of RNA 3d structures.
```
/RNArefiner/bin/rnaoptimizer -i [path of input PDB] -o [directory to output] -bp [path of base pairing information] -st [path of stacking information] [other options]
```
Here, base-paring and base-stacking interaction information must be inputted as follows.  
(Other options)  
| Option | Description | Default |
| :-------------: | :-------------: | :-------------: |
| -libdir | Directory of user-defined template structures (All PDB files in the directory are loaded) | None |
| -libwg | Weights of an user-defined template structures-based energy (-libdir must be specified) | 10 if -libdir is specified, otherwise None |
| -libcwg | Weights of each template structure (-libdir must be specified) | 1/(template number) if -libdir is specified, otherwise None |
| -time | Simulation time (hours) | 10 |

The energy based on the user-defined template structures is calculated as follows.


Specify the directory including the customized template structures by the option of "-libdir".


In the energy calculation, the distance maps between the template structures and the refining structure.
If you emphasize the effects of template structure, set a higher weight by the option of "-libwg".
In addition, you can assign the different 

$$ E_{\text{template}} = w_{\text{template}} \times {\sum_{i=1}^{n} u_{\text{template},i} \times E_{\text{template},i}} $$
$$ x = {-b \pm \sqrt{b^2-4ac} \over 2a} $$
  
  





