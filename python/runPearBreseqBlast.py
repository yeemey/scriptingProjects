#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 16:48:53 2016

python 3.5.2
BioPython 1.68

@author: ymseah
"""

from blastUnmapped import run_software, parse_results
import re, os, subprocess

run_object = run_software()
run_object.batch_run_pear('/Users/ymseah/scriptingProjects/python/test/', '/Users/ymseah/scriptingProjects/')


'''
  
#test_object.run_breseq('/pear_output_dir/', 'Ancestor', '/Users/ymseah/scriptingProjects/python/')
#test_object.fastq_to_fasta("sic_M1_pear.assembled.unmatched.fastq","sic_M1_pear.assembled.unmatched.fasta")
#test_object.parse_blast("sic_M1_pear.assembled.unmatchedBLAST.out")
'''