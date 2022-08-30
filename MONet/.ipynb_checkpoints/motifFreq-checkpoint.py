import numpy as np
import pandas as pd
import sys
import re
import os
from Bio import motifs, SeqIO, Seq
from itertools import combinations
import multiprocessing
from scipy import stats
import Bio
import random
from Bio.SeqUtils import GC


"""
This module takes input from a fa file, and compute the co-occurence matrix (n_peak, n_tf).
The input can be background peaks (all tissue cell types) or peaks of interest (e.g. fetal microglia). 

"""
def readFasta(fasta_file, skip_duplicate=True, fmt='fasta'):
    ############
    # Read in sequences
    ############
    alphabet = Bio.Seq.IUPAC.Alphabet.IUPAC.IUPACUnambiguousDNA() # need to use this alphabet for motif score calculation
    id_seq_dict = {} # {sequenceID: fastq sequence}
    duplicate_keys = []
    for seq_record in SeqIO.parse(fasta_file, fmt):  
        seq_record.seq.alphabet = alphabet
        if seq_record.id in id_seq_dict.keys():
            duplicate_keys.append(seq_record.id)
        else:
            id_seq_dict[seq_record.id] = seq_record.seq
    # delete duplicate keys
    if skip_duplicate:
        for dk in duplicate_keys:
            del id_seq_dict[dk]
        if len(duplicate_keys) > 0:
            print('Ignore duplicate keys in %s: %s' % (fasta_file, duplicate_keys))
    return id_seq_dict

class uniqueTF:
    ############
    # This class computes non-overlapping motifs for ONE peak
    # Input mtx: np matrix e.g. (#tf, peak_size)  e.g.(401,300)
    # Input tfsizes: dict of transcription sizes e.g. {SPI1: 15....}
    # Input motif_cutoff: dict of motif cutoff scores e.g. {SPI1:7.82....}
    # Input tfnames: tfnames to use, list, can be subseted from full list e.g. 401 tfs to 70tfs .... 
    # To run: myobj = uniqueTF(mtx,tfsizes,motif_cutoff,tfnames)
    #         peak_onehot_short=peak_onehot_obj.repeatFillGaps()
    # Return: np array of 0 or 1 to indicate if TF is bind to the current peak
    # It should be noted the peak_onehot_short is (#tfsubset,1)
    ############
    def __init__(self, mtx, tfsizes, motif_cutoff,tfnames):
        self.mtx=mtx
        self.tfsizes=tfsizes
        self.motif_cutoff=motif_cutoff
        self.tfnames=tfnames
        self.peak_size=self.mtx.shape[1]
        self.track_array =  np.empty((self.peak_size,4))
        self.track_array[:] = np.nan
        self.tf_onehot = None
    def structure(self,A):
        """
        check for np.nan in array, return starting and length of continuous np.nan 
        https://stackoverflow.com/questions/64266229/fast-way-to-find-length-and-start-index-of-repeated-elements-in-array
        """
        A=A*1
        ZA = np.concatenate(([0], A, [0]))
        indices = np.flatnonzero( ZA[1:] != ZA[:-1] )
        counts = indices[1:] - indices[:-1]

        #return np.array([indices[::2], counts[::2]])
        return list(zip(indices[::2], counts[::2]))
        #return list(zip(indices[::2], counts[::2])), np.array([list(a) for a in zip(indices[::2], counts[::2])])
    def checkTfMax(self,mtx_sub, tfi, pos,tf_indexs_good_size):
        # input: mtx_sub: a matrix of given bp range, with all tfs mtx_sub = mtx[:,bp_start:bp_end+1]
        # input: tfi: index of tf amount 401, a single int
        # input pos: relative position, pos_300= pos + bp_start 
        # input tf_indexs_good_size: [idxs..] of tfs with size meet the requirement
        # outif T/F whether the tf is the best tf for the tf list with good sizes
        # pos >= 1 and <=299
        #tf_sizes=sizes[tf_indexs_good_size]
        tf_size=self.tfsizes[tfi]  # get tf size
        if pos-tf_size > -1:
            #legit_other_tf_max = [mtx_sub[tf_indexs_good_size[i],pos-sizes[tf_indexs_good_size[i]] +1 : pos+1] for i in tf_indexs_good_size if pos-sizes[tf_indexs_good_size[i]]>1] 
            legit_other_tf_max = [mtx_sub[i,pos-self.tfsizes[i] +1 : pos+1] for i in tf_indexs_good_size if pos-self.tfsizes[i]> -1] 
            #legit_other_tf_max=[]
            #for i in tf_indexs_good_size:
                #if pos - self.tfsizes[i]> -1:
                    #legit_other_tf_max.append(mtx_sub[i,pos-self.tfsizes[i] +1 : pos+1])
            legit_other_tf_max_np=np.concatenate(legit_other_tf_max).max()

            return (True if mtx_sub[tfi,pos]>=legit_other_tf_max_np else False)
        else:
            return False
            #max_pre=mtx_1[:,: pos].max()
            #return (True if self.mtx[tfi,pos]>max_pre else False)
    def fillGaps(self,tf_index,gap_tup):
        ## need to know which tfs meet the size requirements
        ## tf_indexs [1,2,3....] tfs that fits the size requirement
        ## gap_tup (gap_pos,gap_length)
        # 1. if tf_index is empty, means no tf less than gap size, return nothing 
        # 2. 
        if len(tf_index) ==0:
            return
        else:
            sub_mtx=self.mtx[:,gap_tup[0]:gap_tup[0]+gap_tup[1]]  # only tf meet size and gaps [:,gap_start:gap_end]
            max_argsort=np.argsort(sub_mtx,axis=0)   # (gap,)
            for i in range(max_argsort.shape[1]):
                tfi=[idx for idx in max_argsort[:,i] if idx in tf_index][-1]
                current_300_pos = i + gap_tup[0]
                if self.checkTfMax(sub_mtx,tfi,i,tf_index):
                    tf_selected_size= self.tfsizes[tfi]
                    if np.all(np.isnan(self.track_array[current_300_pos-tf_selected_size +1 ,:] )) :  # empty, simply replace the whole thing
                        motif_score=self.mtx[tfi,current_300_pos]
                        motif_in=np.hstack((np.repeat([[motif_score,tfi,tf_selected_size]],tf_selected_size,axis=0),np.arange(tf_selected_size)[:, np.newaxis]))
                        self.track_array[current_300_pos-tf_selected_size +1: current_300_pos+1,:]=motif_in
                        #track_array[current_300_pos-tf_selected_size +1: current_300_pos+1,:] = np.repeat([[motif_score,tfi,tf_selected_size]],tf_selected_size,axis=0)
                    else:
                        need_to_na_size= int(self.track_array[current_300_pos-tf_selected_size +1,3])
                        self.track_array[current_300_pos-tf_selected_size +1 - need_to_na_size :current_300_pos-tf_selected_size +1,:] =np.nan
                        motif_score=self.mtx[tfi,current_300_pos]
                        motif_in=np.hstack((np.repeat([[motif_score,tfi,tf_selected_size]],tf_selected_size,axis=0),np.arange(tf_selected_size)[:, np.newaxis]))
                        self.track_array[current_300_pos-tf_selected_size +1: current_300_pos+1,:]=motif_in
    def repeatFillGaps(self):
        gap_tu=self.structure(np.isnan(self.track_array[:,0])) # init with list [(300,0)]
        gap_tu_pre=list() # init with empty list 
        gap_tu_diff=[tu for tu in gap_tu if tu not in gap_tu_pre]
        while (len(gap_tu_diff)>0) and (max([i[1] for i in gap_tu_diff])>=self.tfsizes.min()):
            gap_tu_diff_sub= [tu for tu in gap_tu_diff if tu[1]>self.tfsizes.min()]
            for tu in gap_tu_diff_sub:
                tf_sub=np.where(self.tfsizes<=tu[1])[0]
                self.fillGaps(tf_sub,tu)
            gap_tu_pre=gap_tu.copy()
            gap_tu=self.structure(np.isnan(self.track_array[:,0]))
            gap_tu_diff=[tu for tu in gap_tu if tu not in gap_tu_pre]
        #unq=list(set(np.unique(self.track_array[~np.isnan(self.track_array[:,0]),:],axis=0)[:,1].astype(int)))
        #self.tf_onehot=np.zeros(len(self.tfsizes))
        #self.tf_onehot[unq]=1
        ta_sub=self.track_array[:,:3]
        unq=np.unique(ta_sub[~np.isnan(ta_sub[:,0]),:],axis=0)   #only keep unique rows
        tf_unq=np.array(self.tfnames)[unq[:,1].astype(int)]
        cutoff_unq=[self.motif_cutoff[tf] for tf in tf_unq]
        unq_pass_cutoff=unq[np.where(unq[:,0]-cutoff_unq >0),:][0,...]
        self.tf_onehot=np.zeros(len(self.tfnames))
        self.tf_onehot[unq_pass_cutoff[:,1].astype(int)]=1
        return np.array(self.tf_onehot) 

class motifFreq:
    def __init__(self,fa_path_A, ncore=30, subset_factor=1, tfexp=None, 
                 thres_path="../data/HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt",
                 homer_path="../data/HOCOMOCOv11_core_HUMAN_mono_homer_format_0.0001.motif"):
        self.fa_path_A=fa_path_A
        self.thres_path=thres_path
        self.homer_path=homer_path
        self.subset_factor=int(subset_factor)
        self.tfexp=tfexp 
        self.tfnames=None
        self.motif_pssm=None
        self.motif_cutoff=None
        self.seq_dict_A= None
        self.peak_id_A=None
        self.mtxA=None
        self.pairs= None
        self.ncore=ncore
        self.tfsizes=None
        self.mtx=None
        """
        full TF 
        """
        self.tfnames_full=None # sorted list of TF names etc [AP1,.....]
        self.tfsizes_full=None # np array, index correspond to tfname sorted, array[10,12....]
        
        self.subset_tf_index=None
        
        
    def loadMotifLib(self):
        thres_df=pd.read_csv(self.thres_path,sep='\t')   # df containing different cutoff
        thres_df_tf=thres_df.apply(lambda x: x[0].split('_')[0],axis=1) # pd series containing tf names 401
        thres_df_thres=np.array(thres_df.iloc[:,-1])  # last columns of threshold df containing cutoffs for
        self.motif_cutoff=dict(zip(thres_df_tf,thres_df_thres))
        self.motif_pssm={}
        with open(self.homer_path) as handle:
            for m in motifs.parse(handle, "pfm-four-columns"):
                tfname=m.name.split()[1].split("_")[0]
                pssm=m.pssm
                self.motif_pssm[tfname]=pssm
        self.tfnames_full=sorted(list(self.motif_cutoff.keys()))
        self.tfsizes_full=np.array([len(self.motif_pssm[tf]['G']) for tf in self.tfnames_full])
        if self.tfexp is None:
            self.tfnames=self.tfnames_full
            self.tfsizes=self.tfsizes_full
            self.subset_tf_index= list(range(len(self.tfnames)))
            #self.motif_pssm=self.motif_pssm_full
            #self.motif_cutoff=self.motif_cutoff_full
        else:
            if isinstance(self.tfexp, (np.ndarray, list))==False:
                raise TypeError("TF list is expected to be a list")
            elif len(self.tfexp)==0:
                raise ValueError("TF list length is zero")
            else:
                self.tfexp=sorted(list(self.tfexp))
                self.subset_tf_index = [self.tfnames_full.index(tf_selected) for tf_selected in self.tfexp if tf_selected in self.tfnames_full]
                self.tfnames = [self.tfnames_full[idx] for idx in self.subset_tf_index]
                self.tfsizes = np.array([self.tfsizes_full[idx] for idx in self.subset_tf_index])
                #self.motif_pssm=self.motif_pssm_full
                #self.motif_cutoff=self.motif_cutoff_full
                unmapped_tf = [tf_unselected for tf_unselected in self.tfexp if tf_unselected not in self.tfnames_full]
                print("%d / %d selcted TF symbol not mapped" % (len(unmapped_tf),len(self.tfexp)))
                print("unmapped TF symbols are:")
                print(", ".join(unmapped_tf))
            
            

        
    def loadSeq(self):
        self.seq_dict_A= readFasta(self.fa_path_A)  
        #subset peaks to 1/3
        temp_peak_id_A=sorted([key for key, value in self.seq_dict_A.items() ])
        sub_idx=[i for i in range(len(temp_peak_id_A)) if i%self.subset_factor==0]
        self.peak_id_A=[temp_peak_id_A[i] for i in range(len(temp_peak_id_A)) if i in sub_idx]
        
        #self.peak_id_A=self.peak_id_A[:20]

    def scoreMotifBoth(self,peakID,motif_name):
        '''
        compute highest score for current motif subtract by cutoff
        input: one peak id , one motif 
        output: np.array shape (300,) position-wise max score for given motif scores
        '''

        seq=self.seq_dict_A[peakID]
        fwd_pssm = self.motif_pssm[motif_name]
        rev_pssm = fwd_pssm.reverse_complement()
        # calculate motif size -1
        m_size1 = len(fwd_pssm['G']) -1


        fwd_scores = fwd_pssm.calculate(seq) # scores for forward orientation
        score_to_insert_fwd= fwd_scores[0]
        fwd_scores300=np.insert(fwd_scores,0,[score_to_insert_fwd]*m_size1)

        rev_scores = rev_pssm.calculate(seq) # scores for reverse orientation
        score_to_insert_rev= rev_scores[0]
        rev_scores300=np.insert(rev_scores,0,[score_to_insert_rev]*m_size1)

        max_scores= np.amax([fwd_scores300,rev_scores300],axis=0)

        return max_scores
    def allMotifScores(self,peakID):
        '''
        for each peak, compute all 401 motifs status
        return a (401,) array
        '''
        #all_motif_one_seq_scores=[self.scoreMotifBoth(peakID,motif_name) for motif_name in self.tfnames]
        self.mtx=np.array([self.scoreMotifBoth(peakID,motif_name) for motif_name in self.tfnames])
        peak_onehot_obj=uniqueTF(self.mtx,self.tfsizes,self.motif_cutoff,self.tfnames)
        peak_onehot_short=peak_onehot_obj.repeatFillGaps()
        # convert short one hot (# selected tf,1) to full one hot (401,1)
        full_tf_onehot=np.zeros((len(self.tfnames_full,)))
        whereto_one=np.array(self.subset_tf_index)[np.where(peak_onehot_short==1)[0]]  # index out of 401 where there is a co-occurence event
        full_tf_onehot[whereto_one]=1                                
        return np.array(full_tf_onehot) 
 
    
    def allPeaksMotifsScore(self):
        '''
        compute all peaks, motif score, returns npeak x 401 motif array
        
        '''
        self.loadMotifLib()
        self.loadSeq()
        with multiprocessing.Pool(self.ncore) as pool:
            self.mtxA=np.array(list(pool.map(self.allMotifScores,self.peak_id_A)))
        return self.mtxA


class motifFreq_overlap:
    def __init__(self,fa_path_A, ncore=30, subset_factor=1, tfexp=None, 
                 thres_path="/home/zhl022/daima/projects/Motif_coocc_Mar2022/motif_lib/HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt",
                 homer_path="/home/zhl022/daima/projects/Motif_coocc_Mar2022/motif_lib/HOCOMOCOv11_core_HUMAN_mono_homer_format_0.0001.motif"):
        self.fa_path_A=fa_path_A
        self.thres_path=thres_path
        self.homer_path=homer_path
        self.subset_factor=int(subset_factor)
        self.tfexp=tfexp 
        self.tfnames=None
        self.motif_pssm=None
        self.motif_cutoff=None
        self.seq_dict_A= None
        self.peak_id_A=None
        self.mtxA=None
        self.pairs= None
        self.ncore=ncore
        self.tfsizes=None
        self.mtx=None
        """
        full TF 
        """
        self.tfnames_full=None # sorted list of TF names etc [AP1,.....]
        self.tfsizes_full=None # np array, index correspond to tfname sorted, array[10,12....]
        
        self.subset_tf_index=None
        
        
    def loadMotifLib_overlap(self):
        thres_df=pd.read_csv(self.thres_path,sep='\t')   # df containing different cutoff
        thres_df_tf=thres_df.apply(lambda x: x[0].split('_')[0],axis=1) # pd series containing tf names 401
        thres_df_thres=np.array(thres_df.iloc[:,-1])  # last columns of threshold df containing cutoffs for
        self.motif_cutoff=dict(zip(thres_df_tf,thres_df_thres))
        self.motif_pssm={}
        with open(self.homer_path) as handle:
            for m in motifs.parse(handle, "pfm-four-columns"):
                tfname=m.name.split()[1].split("_")[0]
                pssm=m.pssm
                self.motif_pssm[tfname]=pssm
        self.tfnames_full=sorted(list(self.motif_cutoff.keys()))
        self.tfsizes_full=np.array([len(self.motif_pssm[tf]['G']) for tf in self.tfnames_full])
        if self.tfexp is None:
            self.tfnames=self.tfnames_full
            self.tfsizes=self.tfsizes_full
            self.subset_tf_index= list(range(len(self.tfnames)))
            #self.motif_pssm=self.motif_pssm_full
            #self.motif_cutoff=self.motif_cutoff_full
        else:
            if isinstance(self.tfexp, (np.ndarray, list))==False:
                raise TypeError("TF list is expected to be a list")
            elif len(self.tfexp)==0:
                raise ValueError("TF list length is zero")
            else:
                self.tfexp=sorted(list(self.tfexp))
                self.subset_tf_index = [self.tfnames_full.index(tf_selected) for tf_selected in self.tfexp if tf_selected in self.tfnames_full]
                self.tfnames = [self.tfnames_full[idx] for idx in self.subset_tf_index]
                self.tfsizes = np.array([self.tfsizes_full[idx] for idx in self.subset_tf_index])
                #self.motif_pssm=self.motif_pssm_full
                #self.motif_cutoff=self.motif_cutoff_full
                unmapped_tf = [tf_unselected for tf_unselected in self.tfexp if tf_unselected not in self.tfnames_full]
                print("%d / %d selcted TF symbol not mapped" % (len(unmapped_tf),len(self.tfexp)))
                print("unmapped TF symbols are:")
                print(", ".join(unmapped_tf))
    def loadSeq_overlap(self):
        self.seq_dict_A= readFasta(self.fa_path_A)  
        #subset peaks to 1/3
        temp_peak_id_A=sorted([key for key, value in self.seq_dict_A.items() ])
        sub_idx=[i for i in range(len(temp_peak_id_A)) if i%self.subset_factor==0]
        self.peak_id_A=[temp_peak_id_A[i] for i in range(len(temp_peak_id_A)) if i in sub_idx]
    def scoreMotif_overlap(self,peakID,motif_name):
        '''
        for each peak and each TF, check if motif is above threshold
        output: 1 if motif is present in the peak, 0 if not
        '''
        #self.loadMotifLib()
        #self.loadSeq()
        

        seq_dict=self.seq_dict_A

        seq=seq_dict[peakID]
        fwd_pssm = self.motif_pssm[motif_name]
        rev_pssm = fwd_pssm.reverse_complement()
        fwd_scores = np.array(fwd_pssm.calculate(seq)) # scores for forward orientation
        rev_scores = np.array(rev_pssm.calculate(seq)) # scores for reverse orientation

        if_motif= 1 if np.max(np.concatenate((fwd_scores, rev_scores)))>=self.motif_cutoff[motif_name] else 0
        return if_motif
    def allMotifScore_overlap(self,peakID):
        '''
        for each peak, compute all 401 motifs status
        return a (401,) array
        '''
        all_motif_one_seq_scores=[self.scoreMotif_overlap(peakID,motif_name) for motif_name in self.tfnames]
        return np.array(all_motif_one_seq_scores) 
 
    
    def allPeaksMotifsScore_overlap(self):
        '''
        compute all peaks, motif score, returns npeak x 401 motif array
        
        '''
        self.loadMotifLib_overlap()
        self.loadSeq_overlap()
        with multiprocessing.Pool(self.ncore) as pool:
            self.mtxA=np.array(list(pool.map(self.allMotifScore_overlap,self.peak_id_A)))
        return self.mtxA


"""
This module takes input from a fa file, and compute the GC ratios out of 100,
    and output a numpy array object (.npy). 

"""
class calcGC:
    def __init__(self,fa_path_A, subset_factor=1 ):
        
        self.fa_path_A=fa_path_A
        self.subset_factor=subset_factor
        self.peak_id_A=None
        self.gc_100=None
        
    def computeGC(self):
        self.seq_dict_A= readFasta(self.fa_path_A)  
        #subset peaks to 1/3
        temp_peak_id_A=sorted([key for key, value in self.seq_dict_A.items() ])
        sub_idx=[i for i in range(len(temp_peak_id_A)) if i%self.subset_factor==0]
        self.peak_id_A=[temp_peak_id_A[i] for i in range(len(temp_peak_id_A)) if i in sub_idx]
        self.gc_100 = np.array([GC(self.seq_dict_A[peak ] )for peak in self.peak_id_A])
        return self.gc_100


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
    
    def runEdge(self):

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