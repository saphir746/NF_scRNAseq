#!/usr/bin/env python3

import scanpy as sc
import scrublet as scr
import scipy
import numpy as np
import os
import re
import sys

dir='/camp/stp/babs/working/schneid/projects/sahaie/giovanni.giangreco/Characterisation_of_CAF_in_HPV_cancer_scrnaseq/Data_interim_files/'

#Files = os.listdir(dir)
#Files = [f for f in Files if f.endswith('_raw.h5ad')]

f=sys.argv[1]

def scrub_doublets(f):
	#filename=dir+f
	adata = sc.read(f) #filename)
	count_matrix=adata.X.tocsc()
	scrub = scr.Scrublet(count_matrix, expected_doublet_rate=0.06)
	doublet_scores, predicted_doublets = scrub.scrub_doublets()
	adata.obs["scrublet_score"]=doublet_scores
	adata.obs["scrublet_pred"]=predicted_doublets
	filename_new=re.sub('raw', 'scrubletted', f)
	adata.write(filename_new)

scrub_doublets(f)

#[scrub_doublets(f) for f in Files]
