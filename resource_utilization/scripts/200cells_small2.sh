#!/bin/bash
#$ -N cells200_small2            # job name
#$ -l h_rt=6:00:00               # run time
#$ -l h_vmem=60G                 # virtual memory

#$ -o /home/unix/sjohri/valab_sjohri/projects/beanie_private/profiling/errors/200cells_small2.out
#$ -e /home/unix/sjohri/valab_sjohri/projects/beanie_private/profiling/errors/200cells_small2.err

source /broad/software/scripts/useuse
reuse Anaconda3
source activate /home/unix/sjohri/valab_sjohri/envs/beanie_revisions
python /home/unix/sjohri/valab_sjohri/projects/beanie_private/profiling/scripts/200cells_small2.py