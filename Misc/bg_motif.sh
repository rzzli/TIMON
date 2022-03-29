#!/bin/bash
# This script run the motif counts for ENCODE samples, used as background
# Dendritic cell samples were not included in this script, analyzed differently

source /gpfs/data01/glasslab/home/zhl022/anaconda3/etc/profile.d/conda.sh
conda activate rickli
pscript=/home/zhl022/daima/projects/Motif_coocc_Mar2022/scripts/run_celltype20_ncore_fac.py
for cell in $(ls -d /home/zhl022/daima/projects/Motif_coocc_Mar2022/background_atac_tissue/enhancer_promoter_split_fa/*)
do 
    bcell=$(basename $cell)
    #echo $bcell
    for fa in $(ls $cell/*fa)
    do
        echo $fa
        mkdir -p /home/zhl022/daima/projects/Motif_coocc_Mar2022/background_atac_tissue/npy1/${bcell}/
        /home/zhl022/.conda/envs/rickli/bin/python $pscript $fa $bcell /home/zhl022/daima/projects/Motif_coocc_Mar2022/background_atac_tissue/npy1/ 20 3
        echo "done"
    done
    
done