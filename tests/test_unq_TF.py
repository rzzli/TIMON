#!/usr/bin/env python

import sys
sys.path.append('../')
from Moe import motifFreq, readFasta, uniqueTF
import numpy as np
import unittest

class TestNonOverlap(unittest.TestCase):
    def setUp(self):
        self.obj=motifFreq("./fixtures/fixture_10seq.fa",ncore=10,subset_factor=1)
        self.obj.loadMotifLib()
        self.obj.loadSeq()
        self.test_peak_id = self.obj.peak_id_A[9]
    def test_gap(self):
        self.mtx=np.array([self.obj.scoreMotifBoth(self.test_peak_id,motif_name) for motif_name in self.obj.tfnames])
        peak_onehot_obj=uniqueTF(self.mtx,self.obj.tfsizes,self.obj.motif_cutoff,self.obj.tfnames)
        peak_onehot_short=peak_onehot_obj.repeatFillGaps()
        gap_tu=peak_onehot_obj.gap_tu # list of tuples of gaps e.g.[(0, 28), (48, 1), (63, 9)...]
        track=peak_onehot_obj.track_array[:,:3]  # a (300bp,3) array, first column motif score, second column tfidx, third column tfsize
        for tu in gap_tu:
            #track
            self.assertTrue(np.all(np.isnan(track[tu[0]:tu[0]+tu[1],0])))

 
if __name__ == '__main__':
    unittest.main()