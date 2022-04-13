#!/usr/bin/env python

import sys
sys.path.append('../')
from Moe import calcGC
import numpy as np
import unittest

class TestNonOverlap(unittest.TestCase):
    def setUp(self):
        self.obj=calcGC("./fixtures/fixture_10seq.fa",subset_factor=1)
        self.gc_array=self.obj.computeGC()
    def test_gc_output_shape(self):
        self.assertEqual(self.gc_array.shape,(10,))
 
    
if __name__ == '__main__':
    unittest.main()