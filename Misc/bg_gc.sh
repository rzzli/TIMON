#!/bin/bash


source /gpfs/data01/glasslab/home/zhl022/anaconda3/etc/profile.d/conda.sh
conda activate rickli
pscript=/home/zhl022/daima/projects/Motif_coocc_Mar2022/scripts/calcGC.py
for cell in $(ls -d /home/zhl022/daima/projects/Motif_coocc_Mar2022/background_atac_tissue/enhancer_promoter_split_fa/*)
do 
    bcell=$(basename $cell)
    #echo $bcell
    for fa in $(ls $cell/*fa)
    do
        echo $fa
        mkdir -p /home/zhl022/daima/projects/Motif_coocc_Mar2022/background_atac_tissue/npy1_GC/${bcell}/
        /home/zhl022/.conda/envs/rickli/bin/python $pscript $fa $bcell /home/zhl022/daima/projects/Motif_coocc_Mar2022/background_atac_tissue/npy1_GC/ 3
        echo "done"
    done
    
done