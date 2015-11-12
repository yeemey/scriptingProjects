# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 12:59:47 2015

@author: ymseah
"""

import numpy as np
import matplotlib.pyplot as plt
import itertools # for zipping, less memory consumption than zip

filepath = raw_input("Enter path to SDT matrix file: ")
filepath2 = raw_input("Enter path to TreePuzzle distance output file: ")

# function to get sequence IDs from SDT matrix output file
def getSDTid(filepath):
    with open(filepath) as sdtfile:
        seqid = []
        for lines in sdtfile:
            eachline = lines.split(',')
            seqid.append(eachline[0])
    return seqid

# function to calculate distance from similarity, and vice-versa (complement)
def comp(num):
    compnum = 1-num
    return compnum

# function to list all SDT distances
def getSDTdist(filepath):
    with open(filepath) as sdtsims:
        sdtdist = []
        for lines in sdtsims:
            simsline = lines.split(',')
            similarities = simsline[1:(len(simsline)-1)]
            for item in similarities:
                item = item.strip('\n')
                item = float(item)
                sdtdist.append(comp(item))
    return sdtdist    

# function to get pairwise uncorrected distances from SDT output file
# key (tuples of sequence IDs): value (comp(SDT identities))
def getPUD(filepath):
    SDTid = getSDTid(filepath)
    SDTdist = getSDTdist(filepath) # distance values
    seqpairs = [] # tuple keys
    i = 1
    while i<len(SDTid):
        j = 0
        while j < i:
            eachseqpair = (SDTid[i], SDTid[j])
            seqpairs.append(eachseqpair)
            j += 1
        i += 1
    PUDS = dict(itertools.izip(seqpairs, SDTdist))
    return PUDS

# function to get a list of lists comprising sequence ID and TreePuzzle distances
def getTPZiddist(filepath2):
    with open(filepath2) as tpzfile:
        numseq = tpzfile.readline()
        numseq = numseq.strip('\n')
        numseq = int(numseq)
        seqid2 = []
        seqid2element = []
        tpzbyid = []
        for lines in tpzfile:
            eachline2 = lines.split(' ')
            for item in eachline2:
                item = item.strip('\n')
                if item != '':
                    tpzbyid.append(item)
        while len(seqid2) < numseq:
            seqid2element = tpzbyid[:numseq+1]
            seqid2.append(seqid2element)
            del tpzbyid[:numseq+1]
    return seqid2

# function to get pairwise evolutionary distances from TreePuzzle
# key (tuples of sequence IDs): value (TreePuzzle identities))
def getPED(filepath2):
    TPZiddist = getTPZiddist(filepath2)
    TPZdistperseq = [] #list of lists of distances for each sequence
    TPZdistances = []
    seqpairs2 = [] # tuple keys
    i = 0
    while i < len(TPZiddist):
        j = 0
        while j < len(TPZiddist):
            eachseqpair2 = (TPZiddist[i][0], TPZiddist[j][0])
            seqpairs2.append(eachseqpair2)
            j += 1
        i += 1
    for entries in TPZiddist:
        TPZdistperseq.append(entries[1:])
    for distances in TPZdistperseq:
        ii = 0
        while ii < len(TPZdistperseq):
            distances[ii] = float(distances[ii])
            ii += 1
        TPZdistances = TPZdistances + distances
    
    PEDSinit = dict(itertools.izip(seqpairs2, TPZdistances))
    PEDS = {}
    for key,value in PEDSinit.iteritems():
        if key[0] != key[1]:
            PEDS[key] = value
    return PEDS

PUD = getPUD(filepath)
PED = getPED(filepath2)

