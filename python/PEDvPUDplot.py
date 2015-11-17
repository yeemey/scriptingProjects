# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 12:59:47 2015

Make scatterplot of pairwise uncorrected vs evolutionary distances.

Pairwise uncorrected distances are obtained from SDT 1.0 (filename_mat.txt), while pairwise maximum likelihood distances are obtained from Tree-Puzzle 5.2 (outdist.txt). Both files are read to extract sequence ID and distance values.

"""

import numpy as np
import matplotlib.pyplot as plt
import itertools

filepath = raw_input("Enter path to SDT matrix file: ")
filepath2 = raw_input("Enter path to TreePuzzle distance output file: ")

def getSDTid(filepath):
    """
    Reads SDT 1.0 matrix output (filename_mat.txt) and returns list of sequence IDs.
    
    ::param filepath::  Path to SDT file (input).
    ::param sdtfile::   Input file.
    ::param seqid::     List of sequence ID strings.
    ::param eachseqid:: Individual sequence ID string.
    ::param line::      Line from input file.
    ::param eachline::  List of elements from param line.
    """
    with open(filepath) as sdtfile:
        seqid = []
        eachseqid = ""
        for line in sdtfile:
            eachline = line.split(',')
            eachseqid = eachline[0]
            if len(eachseqid) > 10:
                eachseqid = eachseqid[:10]
                seqid.append(eachseqid)
            else:
                seqid.append(eachseqid)
    return seqid

def getSDTdist(filepath):
    """
    Reads SDT 1.0 matrix output (filename_mat.txt) and returns list of pairwise distances as float numbers. 
    
    Order of list elements corresponds to order of distances (left to right) in each matrix row (top to bottom) in the SDT file. Pairwise identities of 0 between the same pairwise sequences are excluded. Pairwise distances are calculated as the complement of identities.
    
    ::param filepath::     Path to SDT file (input).
    ::param sdtsims::      Input file.
    ::param sdtdist::      List of pairwise distances (float).
    ::param line::         Line from input file.
    ::param simsline::     List of elements from param line.
    ::param similarities:: Slice of param simsline to exclude sequence ID & distance of 0 (pairwie comparison with self).
    ::param item::         Individual distance measure.
    """
    with open(filepath) as sdtsims:
        sdtdist = []
        for line in sdtsims:
            simsline = line.split(',')
            similarities = simsline[1:(len(simsline)-1)]
            for item in similarities:
                item = item.strip('\n')
                item = float(item)
                sdtdist.append(1-item)
    return sdtdist    

def getPUD(filepath):
    """
    Calls functions getSDTid and getSDTdist and returns dictionary of pairwise sequence (tuple) keys and pairwise uncorrected distance values.
    
    ::param filepath::     Path to SDT file (input).
    ::param SDTid::       List of sequence ID strings.
    ::param SDTdist::     List of dictionary values (float distances).
    ::param seqpairs::    List of dictionary keys (tuple of sequence IDs).
    ::param seq1Count::   Var to keep count of 1st sequence ID in pair.
    ::param seq2Count::   Var to keep count of 2nd sequence ID in pair.
    ::param eachseqpair:: Tuple of sequence IDs
    ::param PUDS::        Dictionary of pairwise sequence keys and uncorrected distance values, excluding pairwise comparisons of sequence with itself.
    """
    SDTid = getSDTid(filepath)
    SDTdist = getSDTdist(filepath) 
    seqpairs = [] 
    seq1Count = 1
    while seq1Count<len(SDTid):
        seq2Count = 0
        while seq2Count < seq1Count:
            eachseqpair = (SDTid[seq1Count], SDTid[seq2Count])
            seqpairs.append(eachseqpair)
            seq2Count += 1
        seq1Count += 1
    PUDS = dict(itertools.izip(seqpairs, SDTdist))
    return PUDS

def getTPZiddist(filepath2):
    """
    Reads Tree-Puzzle 5.2 maximum likelihood distance file (outdist.txt) and returns list of lists comprising sequence id (first element) and all pairwise distances.
    
    ::param filepath2::      Path to Tree-Puzzle file (input).
    ::param tpzfile::        Input file.
    ::param numseq::         Number of sequences in file.
    ::param line::           Line from input file.
    ::param eachline2::      List of elements from param line.
    ::param item::           Each element in param eachline2.
    ::param eachline2clean:: List of elements from param line, minus empty strings.
    ::param seqid2element:: List of sequence ID string (first element) and pairwise distances for that sequence (subsequent elements).
    ::param seqid2::        List of param seqid2element.    
    """
    with open(filepath2) as tpzfile:
        numseq = tpzfile.readline()
        numseq = numseq.strip('\n')
        numseq = int(numseq)
        seqid2 = []
        seqid2element = []
        eachline2clean = []
        for line in tpzfile:
            eachline2 = line.split(' ')
            for item in eachline2:
                item = item.strip('\n')
                if item != '':
                    eachline2clean.append(item)
        while len(seqid2) < numseq:
            seqid2element = eachline2clean[:numseq+1] 
            seqid2.append(seqid2element)
            del eachline2clean[:numseq+1]
    return seqid2

def getPED(filepath2):
    """
    Calls the function getTPZiddist and returns dictionary of pairwise sequence (tuple) keys and pairwise maximum likelihood distance values. 
    
    ::param filepath2::      Path to Tree-Puzzle file (input).
    ::param TPZiddist::      List of lists of sequence ID string (first element) and distances for that sequence (subsequent elements).
    ::param TPZdistperseq::  List of param distanceList.
    ::param distanceList::   List of distances for one sequence.
    ::param TPZdistances::   List of dictionary values (float distances for all sequences.
    ::param seq1Count::      Var to keep count of 1st sequence ID in pair.
    ::param seq2Count::      Var to keep count of 2nd sequence ID in pair.
    ::param eachseqpair2::   Tuple of sequence string pairs.
    ::param seqpairs2::      List of dictionary keys (param eachseqpair2).
    ::param entries::        List of sequence ID string (first element) and distances for that sequence (subsequent elements).
    ::param countDistances:: Var to keep count of distance values per sequence.
    ::param PEDSinit::       Dictionary of pairwise sequence keys and evolutionary distance values.
    ::param PEDS::           Dictionary of pairwise sequence keys and evolutionary distance values, excluding pairwise comparisons of sequence with itself.
    """
    TPZiddist = getTPZiddist(filepath2)
    TPZdistperseq = [] 
    TPZdistances = []
    seqpairs2 = []
    seq1Count = 0
    while seq1Count < len(TPZiddist):
        seq2Count = 0
        while seq2Count < len(TPZiddist):
            eachseqpair2 = (TPZiddist[seq1Count][0], TPZiddist[seq2Count][0])
            seqpairs2.append(eachseqpair2)
            seq2Count += 1
        seq1Count += 1
    for entries in TPZiddist:
        TPZdistperseq.append(entries[1:])
    for distanceList in TPZdistperseq:
        countDistances = 0
        while countDistances < len(distanceList):
            distanceList[countDistances] = float(distanceList[countDistances])
            countDistances += 1
        TPZdistances = TPZdistances + distanceList
    
    PEDSinit = dict(itertools.izip(seqpairs2, TPZdistances))
    PEDS = {}
    for key,value in PEDSinit.iteritems():
        if key[0] != key[1]:
            PEDS[key] = value
    return PEDS

def PUDxPEDy(filepath, filepath2):
    """
    Calls the functions getPUD and getPED and returns distance values for corresponding keys in both dictionaries, as two lists packed in a tuple.
    
    ::param filepath::  Path to SDT file (input).
    ::param filepath2:: Path to Tree-Puzzle file (input).
    ::param PUDdict::   Dictionary of pairwise sequence keys and uncorrected distance values, excluding pairwise comparisons of sequence with itself.
    ::param PEDdict::   Dictionary of pairwise sequence keys and evolutionary distance values, excluding pairwise comparisons of sequence with itself.
    ::param x::         List of evolutionary distances.
    ::param y::         List of uncorrected distances.
    """
    PUDdict = getPUD(filepath)
    PEDdict = getPED(filepath2)
    x = []
    y = []
    for key in PUDdict.keys():
        y.append(PUDdict[key])
        x.append(PEDdict[key])
    return (x,y)

def reportDistances(filepath, filepath2):
    """
    Prints sequence pairs and distances (uncorrected and evolutionary) to file.
    
    ::param filepath::  Path to SDT file (input).
    ::param filepath2:: Path to Tree-Puzzle file (input).
    ::param PUDdict::   Dictionary of pairwise sequence keys and uncorrected distance values, excluding pairwise comparisons of sequence with itself.
    ::param PEDdict::   Dictionary of pairwise sequence keys and evolutionary distance values, excluding pairwise comparisons of sequence with itself.
    ::param report::    Function output written to file.        
    """
    PUDdict = getPUD(filepath)
    PEDdict = getPED(filepath2)
    with open("PUDPED.txt", "w") as report:
        report.write("sequencePair, uncorrectedDistance, evolutionaryDistance \n")
        for key in PUDdict.keys():
            report.write(str(key) + ", " + str(PUDdict[key]) + ", " + str(PEDdict[key]) + "\n")
    return report

# Output pairwise distances
reportDistances(filepath, filepath2)
print "Pairwise distances printed to PUDPED.txt"

# Output scatterplot of pairwise distances
scatterPoints = PUDxPEDy(filepath, filepath2)
plt.scatter(scatterPoints[0], scatterPoints[1])
plt.xlabel('Pairwise evolutionary distance')
plt.ylabel('Pairwise uncorrected distance')
plt.axis([0, np.ceil(max(scatterPoints[0])), 0, np.ceil(max(scatterPoints[1]))])
plt.plot([0, np.ceil(max(scatterPoints[0]))], [0, np.ceil(max(scatterPoints[1]))], 'k:')
plt.savefig("PUDPEDplot.png")
print "Plot saved as PUDPEDplot.png"