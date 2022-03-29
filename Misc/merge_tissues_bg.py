import pandas as pd
import numpy as np
import glob
from os.path import exists
import matplotlib.pyplot as plt
import random
"""
This script merge motif counts in peaks from different tissues into one npy file.

"""
jpgFilenamesList1 = glob.glob('/home/zhl022/daima/projects/Motif_coocc_Mar2022/background_atac_tissue_DC/npy1/*/*enhancer*npy')
jpgFilenamesList2 = glob.glob('/home/zhl022/daima/projects/Motif_coocc_Mar2022/background_atac_tissue/npy1/*/*enhancer*npy')
jpgFilenamesList= jpgFilenamesList1 + jpgFilenamesList2

peak_one_hot_array=np.empty(shape=(0,401))
peak_gc_arr=np.array([])
for peak in jpgFilenamesList:
    gc_peak_file=peak.replace("npy1","npy1_GC").replace(".npy","_GC_perc.npy")
    if (exists(peak)) and (exists(gc_peak_file)):
        with open(peak, 'rb') as f:
            peak_load = np.load(f)
        peak_one_hot_array=np.vstack((peak_one_hot_array,peak_load))
        with open(gc_peak_file, 'rb') as f:
            gc_load = np.load(f)
        peak_gc_arr=np.hstack((peak_gc_arr,gc_load)) 

with open('../all_bg/bg_all_one_hot_all.npy', 'wb') as f:
    np.save(f, peak_one_hot_array)
    np.save(f, bin1_idx)
    np.save(f, bin2_idx)
    np.save(f, bin3_idx)
    np.save(f, bin4_idx)
    np.save(f, bin5_idx)
    np.save(f, bin6_idx)