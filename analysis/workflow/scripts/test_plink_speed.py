#!/usr/bin/env python

import sys
import time
import numpy as np
from pgenlib import PgenReader

variant_ct_start = 3220890
variant_ct_end = 3253426
variant_ct = variant_ct_end - variant_ct_start
pgen = PgenReader(bytes("/projects/ps-gymreklab/resources/datasets/ukbiobank/array_imputed/pfile_converted/chr1.pgen", "utf8"), sample_subset=np.arange(int(sys.argv[1]), dtype=np.uint32))
sample_ct = pgen.get_raw_sample_ct()
data = np.empty((variant_ct, sample_ct * 2), dtype=np.int32)
pgen.read_alleles_range(variant_ct_start, variant_ct_end, data)
