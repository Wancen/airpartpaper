#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 00:45:20 2021

@author: wancen
"""

import sys
sys.path.insert(0, '..')
import numpy as np
from scdali import run_scdali
from scdali.utils.stats import apply_fdr_bh


def practice(ase, cts, x, ncores):
# =============================================================================
#   pvals_join = run_scdali(A=ase, D=cts, cell_state=x, model='scDALI-Joint', base_rate=.5)['pvalues']
#   pvals_corrected = apply_fdr_bh(pvals_join)
#   significant=np.where(pvals_corrected<0.1)[0]
# =============================================================================
  gp_results2 = run_scdali(A=ase, D=cts, cell_state=x, model='GP',base_rate=0.5,gp_kernel='RBF',n_cores=ncores)
  mu = gp_results2['posterior_mean']
  sd = np.sqrt(gp_results2['posterior_var'])
  return(mu,sd)

def testsig(ase, cts, x, ncores):
    pvals_join = run_scdali(A=ase, D=cts, cell_state=x, model='scDALI-Het', base_rate=.5)['pvalues']
    pvals_corrected = apply_fdr_bh(pvals_join)
    #significant=np.where(pvals_corrected<0.1)[0]
    return(pvals_corrected)