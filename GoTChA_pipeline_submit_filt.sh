#!/bin/bash
#BSUB -J gotcha_filt
#BSUB -P acc_MDS 
#BSUB -q premium
#BSUB -n 1
#BSUB -W 3:00
#BSUB -R rusage[mem=130GB]
#BSUB -R span[hosts=1] 
#BSUB -oo %J.stdout
#BSUB -eo %J.stderr
#BSUB -R "select[osmajor!=CENT7]"
#BSUB -L /bin/bash


module load R/4.2.0
module load nlopt

Rscript GoTChA_FastqFilt.R

