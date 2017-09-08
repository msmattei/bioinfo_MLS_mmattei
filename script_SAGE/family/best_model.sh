#BSUB -L /bin/bash
#BSUB -o /best_method.out
#BSUB -e /best_method.err
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -R "rusage[mem=40000]"
#BSUB -M 40194304

module add R/3.3.2


Rscript best_method.R
