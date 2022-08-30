import sys
import os
import numpy as np
from Bio import motifs, SeqIO, Seq
from Bio.SeqUtils import GC


"""
This module takes input from a fa file, and compute the GC ratios out of 100,
    and output a numpy array object (.npy). 



"""

def read_fasta(fasta_file, skip_duplicate=True, fmt='fasta'):
    '''
    Read in sequences
    '''
    alphabet = Seq.IUPAC.Alphabet.IUPAC.IUPACUnambiguousDNA() # need to use this alphabet for motif score calculation
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

class calcGC:
    def __init__(self,fa_path_A, subset_factor=1 ):
        
        self.fa_path_A=fa_path_A
        self.subset_factor=subset_factor
        self.peak_id_A=None
        self.gc_100=None
        
    def computeGC(self):
        self.seq_dict_A= read_fasta(self.fa_path_A)  
        #subset peaks to 1/3
        temp_peak_id_A=sorted([key for key, value in self.seq_dict_A.items() ])
        sub_idx=[i for i in range(len(temp_peak_id_A)) if i%self.subset_factor==0]
        self.peak_id_A=[temp_peak_id_A[i] for i in range(len(temp_peak_id_A)) if i in sub_idx]
        self.gc_100 = np.array([GC(self.seq_dict_A[peak ] )for peak in self.peak_id_A])
        return self.gc_100

#fa=calcGC('/home/zhl022/daima/archived/proj/encode_download/Jul16_allHepG2/out_fa_misc/CEBPD/CEBPD_1_1.fa',subset_factor=3)
#xx=fa.computeGC()

if __name__ == "__main__":

     
    fasta_file = sys.argv[1]
    cell_type = sys.argv[2]
    out_dir = sys.argv[3]
    sub_fac= int(sys.argv[4])
    basename=os.path.basename(fasta_file).replace(".fa","")
    
    ## run this
    myobj=calcGC(fasta_file,subset_factor=sub_fac)
    gc_array=myobj.computeGC()
    save_path=out_dir + '/' + cell_type  + '/' + basename +  "_GC_perc.npy"
    print("saving to: " + save_path)
    with open(save_path, 'wb') as f:
        np.save(f,gc_array)