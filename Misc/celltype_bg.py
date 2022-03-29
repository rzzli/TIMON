#!/home/zhl022/.conda/envs/rickli/bin/python
import sys
import os
import numpy as np
sys.path.append("/home/zhl022/daima/projects/Motif_coocc_Mar2022/scripts/")
from motifFreq import read_fasta, motifFreq
if __name__ == "__main__":
    """
    This script runs the motifFreq file and modules within to compute motif frequencies and save
        the file into a npy format file. 
    """

     
    fasta_file = sys.argv[1]
    cell_type = sys.argv[2]
    out_dir = sys.argv[3]
    n_core = int(sys.argv[4])
    sub_fac= int(sys.argv[5])
    basename=os.path.basename(fasta_file).replace(".fa","")
    
    ## run this
    myobj=motifFreq(fasta_file,ncore=n_core,subset_factor=sub_fac)
    mtx=myobj.all_peaks_motifs_score()
    #print("mtx shape:",mtx.shape)
    save_path=out_dir + '/' + cell_type  + '/' + basename + ".npy"
    print("saving to: " + save_path)
    with open(save_path, 'wb') as f:
        np.save(f,mtx)
    