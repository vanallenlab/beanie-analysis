#!/bin/bash
#$ -N cells100_large            # job name
#$ -l h_rt=8:00:00               # run time
#$ -l h_vmem=60G                 # virtual memory

#$ -o /home/unix/sjohri/valab_sjohri/projects/beanie_private/profiling/errors/100cells_large.out
#$ -e /home/unix/sjohri/valab_sjohri/projects/beanie_private/profiling/errors/100cells_large.err

source /broad/software/scripts/useuse
reuse Anaconda3
source activate /home/unix/sjohri/valab_sjohri/envs/beanie_revisions
python /home/unix/sjohri/valab_sjohri/projects/beanie_private/profiling/scripts/100cells_large.py