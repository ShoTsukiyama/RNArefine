# RNAoptimizer
This package refines RNA 3D structures and predicts base-pairing and base-stacking interactions.  
RNAoptimizer can be used under The MIT License.
(This project was accomplished through a collaboration between the Yang Zhang lab in the National University of Singapore and the Kurata Hiroyuki lab in the Kyushu Institute of Technology.)  
  
## 1. Installation
### 1-1. Setting of third-party software
RNAoptimizer used the following third-party software.  
Install the software in your environment and configure the path in the RNArefiner/config.json file.  
| Software  | Description |
| :------------- | :------------- |
| [PETfold](https://rth.dk/resources/petfold/download.php)  | RNA secondary structure prediction tool |
| [RNAfold](https://www.tbi.univie.ac.at/RNA/)  | RNA secondary structure prediction tool  |
| [Linearfold](https://github.com/LinearFold/LinearFold)  | RNA secondary structure prediction tool |
| [CSSR](https://github.com/pylelab/CSSR)  | Base-pair identification tool  |
| [CT2DB](https://www.tbi.univie.ac.at/RNA/ct2db.1.html)  | Tool to produce dot bracket notation of an RNA secondary structure  |
| [Arena](https://github.com/pylelab/Arena)  | Full-atomic reconstruction of RNA 3d structures  |
In addtion to those packages, configure the python path in the RNArefiner/config.json file. 
  
### 1-2. Downloading of pre-defined parameters
(Not necessary for base-paring and base-stacking prediction)
Download pre-defined parameters from [here](https://1drv.ms/f/c/7340de942a07a197/EkuUvze_mWpJqsx22H2cz1wBSbRjelLEU5s3J6M0EC_mVg?e=d1WekU), and locate them on the /RNArefiner/param.  
  
### 1-3. Setup of refinement package
(Not necessary for base-paring and base-stacking prediction)
Compile on the RNArefiner/src
```
cd RNArefiner/src
make
```

## 2. Base-pairing prediction
A ResNet-based neural network model predicts base-paring (secondary structure) as follows.

  
## 3. Base-stacking prediction
A ResNet-based neural network model predicts base-stacking as follows.

  
## 4. RNA 3d structure refinement

  
## Citation





