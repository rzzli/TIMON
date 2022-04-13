#!/usr/bin/env python

import sys
sys.path.append('../')
from Moe import motifFreq
import numpy as np
import unittest

class TestNonOverlap(unittest.TestCase):
    def setUp(self):
        self.obj=motifFreq("./fixtures/fixture_10seq.fa",ncore=10,subset_factor=1)
        self.motif_mtx=self.obj.allPeaksMotifsScore()
    def test_output_shape(self):
        self.assertEqual(self.motif_mtx.shape,(10,401))
    def test_output_content(self):
        # should be either 0 and 1
        self.assertTrue(np.array_equal(np.unique(self.motif_mtx),np.array([0., 1.])))
        #self.assertEqual(5,6)
    def test_load_lib(self):
        self.assertEqual(self.obj.tfsizes_full.shape,(401, ))
    
if __name__ == '__main__':
    unittest.main()

