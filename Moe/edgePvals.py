import pandas as pd
import numpy as np
import random
from scipy import stats
import multiprocessing
from itertools import combinations
import sys
import re
import os
"""
This module takes two types of input:
1. a background motif co-occurence file, usually I use the pre-computed background peak data from ENCODE
2. data to analyze, for example, fetal microglia peak motif occurence (not cooccurence yet, will compute here)


The module will return p values for each pair of TFs, if using 401 TFs, there will be 80200 combinations. 
"""
class sigEdges:
    def __init__(self, data_motif_mtx, data_gc_arr, bg_path, n_iter=10000,n_core=0,gc_match= True):
        self.bg_path=bg_path
        self.peak_one_hot_array=None
        self.bin1_idx=None
        self.bin2_idx=None
        self.bin3_idx=None
        self.bin4_idx=None
        self.bin5_idx=None
        self.bin6_idx=None

        with open(self.bg_path, 'rb') as f:
            self.peak_one_hot_array=np.load(f)
            self.bin1_idx=np.load(f )
            self.bin2_idx=np.load(f )
            self.bin3_idx=np.load(f )
            self.bin4_idx=np.load(f )
            self.bin5_idx=np.load(f )
            self.bin6_idx=np.load(f )
        self.data_gc_arr=data_gc_arr
        self.data_motif_mtx=data_motif_mtx
        self.bg_coo_mtx=None
        self.data_coo_mtx=None
        self.n_iter=n_iter
        self.n_core = n_core
        self.gc_match = gc_match
        self.tfpairs=None
        
        
    def genIdxGc(self,gc_arr,n_iter):
        """
        Input
        gc_arr: gc_arr computed from data
        n_ter: numbers of iteration needed, default: 10000

        Output
        final_mtx: (n_iter, data_peak_num)

        Note
        This function returns GC matched indencies for randomly 
            draw peaks from background files, for n iterations.

        The output matrix can be used to sample co-occur using drawCo function.
        """

        gc_bin1_num=(gc_arr<30).sum()
        gc_bin2_num=((gc_arr>=30) & (gc_arr<40)).sum()
        gc_bin3_num=((gc_arr>=40) & (gc_arr<50)).sum()
        gc_bin4_num=((gc_arr>=50) & (gc_arr<60)).sum()
        gc_bin5_num=((gc_arr>=60) & (gc_arr<70)).sum()
        gc_bin6_num=(gc_arr>=70).sum()
        bg_bin1_dim=self.bin1_idx.shape[0]
        bg_bin2_dim=self.bin2_idx.shape[0]
        bg_bin3_dim=self.bin3_idx.shape[0]
        bg_bin4_dim=self.bin4_idx.shape[0]
        bg_bin5_dim=self.bin5_idx.shape[0]
        bg_bin6_dim=self.bin6_idx.shape[0]


        rng = np.random.default_rng(12345)
        bin1_idx_mtx= [self.bin1_idx[rng.integers(low=0, high=bg_bin1_dim, size=gc_bin1_num)] for _ in range(self.n_iter) ]
        bin2_idx_mtx= [self.bin2_idx[rng.integers(low=0, high=bg_bin2_dim, size=gc_bin2_num)] for _ in range(self.n_iter) ]
        bin3_idx_mtx= [self.bin3_idx[rng.integers(low=0, high=bg_bin3_dim, size=gc_bin3_num)] for _ in range(self.n_iter) ]
        bin4_idx_mtx= [self.bin4_idx[rng.integers(low=0, high=bg_bin4_dim, size=gc_bin4_num)] for _ in range(self.n_iter) ]
        bin5_idx_mtx= [self.bin5_idx[rng.integers(low=0, high=bg_bin5_dim, size=gc_bin5_num)] for _ in range(self.n_iter) ]
        bin6_idx_mtx= [self.bin6_idx[rng.integers(low=0, high=bg_bin6_dim, size=gc_bin6_num)] for _ in range(self.n_iter) ]

        mtx= [bin1_idx_mtx,bin2_idx_mtx,bin3_idx_mtx,bin4_idx_mtx,bin5_idx_mtx,bin6_idx_mtx]
        final_mtx=np.empty(shape=(self.n_iter, gc_arr.shape[0]))
        for i in range(self.n_iter):
            final_mtx[i,:] = np.concatenate((mtx[0][i],mtx[1][i],mtx[2][i],mtx[3][i],mtx[4][i],mtx[5][i])) 
        return final_mtx.astype(int)
    
    def genIdxNoGc(self,gc_arr,n_iter):
        rng = np.random.default_rng(12345)
        return np.array([rng.integers(low=0, high=self.peak_one_hot_array.shape[0], size=gc_arr.shape[0]) for _ in range(self.n_iter) ])

        
    def coOcc(self,sub_mtx):
        """
        Input
        sub_mtx: (n_peak,n_tf) matrix

        Output
        coocc: (n_tf,n_tf) diag matrix with diag=0
        
        Note
        Computes co-occurrence    
        """
        coocc= np.dot(np.transpose(sub_mtx),sub_mtx)
        np.fill_diagonal(coocc , 0)
        return coocc
    def drawCo(self,sample_iter):
        return np.array([self.coOcc(self.peak_one_hot_array[sample_iter[i,:]]) for i in range(sample_iter.shape[0])])
    def pValEdges(self,tu):
        node1=tu[0]
        node2=tu[1]
        return stats.percentileofscore(self.bg_coo_mtx[:,node1,node2],self.data_coo_mtx[node1,node2])
    
    def run_edge(self):

        self.data_coo_mtx = self.coOcc(self.data_motif_mtx)
        if self.gc_match:
            bg_sample_iter= self.genIdxGc(self.data_gc_arr,self.n_iter)
        else:
            bg_sample_iter= self.genIdxNoGc(self.data_gc_arr,self.n_iter)
        self.bg_coo_mtx = self.drawCo(bg_sample_iter)
        
        n_tf= self.data_motif_mtx.shape[1]
        self.tfpairs= list(combinations(range(n_tf),2))

        # if n_core is less than 1, do not use multi process 
        if self.n_core >=1:
            with multiprocessing.Pool(self.n_core) as P:
                p_list_fetal=list(P.map(self.pValEdges,self.tfpairs))
        else:
            p_list_fetal=list(map(self.pValEdges,self.tfpairs))
        p_list_fetal=np.array(p_list_fetal)
        return p_list_fetal