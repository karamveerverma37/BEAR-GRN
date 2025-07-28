<img width="300" height="300" alt="Untitled design" src="https://github.com/user-attachments/assets/fc550e11-42ee-46e8-804c-904d90231fc7" />

Here we present, BEAR-GRN (Benchmarking Evaluation and Assessment Resource for GRNs), a framework for bechmarking multi-omics based GRN inference methods. The overall workflow includes, calculation of Accuracy, Stability and Scalability of GRN inference methods on different datasets.

#Installation
BEAR-GRN can be installed as the conda environment:

`conda intall`

OR

Using git clone : 
`git clone`

Input: User need a inferred GRN from any new method using the input data provided in `INPUT.DATA` dir.
# Accuracy Tutorial
##Step 1: use add_grn.py to add the infeered GRN to the BEAR-GRN.
`python add_grn.py --dataset K562 --method_name New_method --grn_path /home/jkl/inferred_grn.csv`

##Step 2: Run BEAR-GRN
`Rscript BEAR-GRN.R --metric accuracy --out_dir output`

##Requirements:
###Python 3.9
###R/4.3.2


