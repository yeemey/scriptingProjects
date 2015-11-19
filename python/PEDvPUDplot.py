# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 12:59:47 2015

Make scatterplot of pairwise uncorrected vs evolutionary distances.

Pairwise uncorrected distances are obtained from SDT 1.0 (filename_mat.txt), while pairwise maximum likelihood distances are obtained from Tree-Puzzle 5.2 (outdist.txt). Both files are read to extract sequence ID and distance values.

"""

import numpy as np
import matplotlib.pyplot as plt
import itertools

SdtFilepath = raw_input("Enter path to SDT matrix file: ")
TpzFilepath = raw_input("Enter path to TreePuzzle distance output file: ")

def getSdtId(SdtFilepath):
    """
    Reads SDT 1.0 matrix output (filename_mat.txt) and returns list of sequence IDs.
    
    ::param SdtFilepath::    Path to SDT file (input).
    ::param SdtFile::        Input file.
    ::param sequenceIds::    List of sequence ID strings.
    ::param eachSequenceId:: Individual sequence ID string.
    ::param line::           Line from input file.
    ::param splitLine::      List of elements from param line.
    """
    with open(SdtFilepath) as SdtFile:
        sequenceIds = []
        eachSequenceId = ""
        for line in SdtFile:
            splitLine = line.split(',')
            eachSequenceId = splitLine[0]
            if len(eachSequenceId) > 10:
                eachSequenceId = eachSequenceId[:10]
                sequenceIds.append(eachSequenceId)
            else:
                sequenceIds.append(eachSequenceId)
    return sequenceIds

def getSdtDist(SdtFilepath):
    """
    Reads SDT 1.0 matrix output (filename_mat.txt) and returns list of pairwise distances as float numbers. 
    
    Order of list elements corresponds to order of distances (left to right) in each matrix row (top to bottom) in the SDT file. Pairwise identities of 0 between the same pairwise sequences are excluded. Pairwise distances are calculated as the complement of identities.
    
    ::param SdtFilepath::  Path to SDT file (input).
    ::param SdtFile::      Input file.
    ::param SdtDistances:: List of pairwise distances (float).
    ::param line::         Line from input file.
    ::param splitLine::    List of elements from param line.
    ::param similarities:: Slice of param simsline to exclude sequence ID & distance of 0 (pairwie comparison with self).
    ::param item::         Individual distance measure.
    """
    with open(SdtFilepath) as SdtFile:
        SdtDistances = []
        for line in SdtFile:
            splitLine = line.split(',')
            similarities = splitLine[1:(len(splitLine)-1)]
            for item in similarities:
                item = item.strip('\n')
                item = float(item)
                SdtDistances.append(1-item)
    return SdtDistances    

def getPud(SdtFilepath):
    """
    Calls functions getSDTid and getSDTdist and returns dictionary of pairwise sequence (tuple) keys and pairwise uncorrected distance values.
    
    ::param SdtFilepath::      Path to SDT file (input).
    ::param SdtIds::           List of sequence ID strings.
    ::param SdtDistances::     List of dictionary values (float distances).
    ::param sequencePairs::    List of dictionary keys (tuple of sequence IDs).
    ::param sequence1Count::   Var to keep count of 1st sequence ID in pair.
    ::param sequence2Count::   Var to keep count of 2nd sequence ID in pair.
    ::param eachSequencePair:: Tuple of sequence IDs
    ::param Puds::             Dictionary of pairwise sequence keys and uncorrected distance values, excluding pairwise comparisons of sequence with itself.
    """
    SdtIds = getSdtId(SdtFilepath)
    SdtDistances = getSdtDist(SdtFilepath) 
    sequencePairs = [] 
    sequence1Count = 1
    while sequence1Count<len(SdtIds):
        sequence2Count = 0
        while sequence2Count < sequence1Count:
            eachSequencePair = (SdtIds[sequence1Count], SdtIds[sequence2Count])
            sequencePairs.append(eachSequencePair)
            sequence2Count += 1
        sequence1Count += 1
    Puds = dict(itertools.izip(sequencePairs, SdtDistances))
    return Puds

def getTpzIdDist(TpzFilepath):
    """
    Reads Tree-Puzzle 5.2 maximum likelihood distance file (outdist.txt) and returns list of lists comprising sequence id (first element) and all pairwise distances.
    
    ::param TpzFilepath::      Path to Tree-Puzzle file (input).
    ::param TpzFile::        Input file.
    ::param numberSequences::         Number of sequences in file.
    ::param line::           Line from input file.
    ::param splitLine::      List of elements from param line.
    ::param item::           Each element in param eachline2.
    ::param splitLinesClean:: List of elements from param line, minus empty strings.
    ::param eachSequenceIdDistance:: List of sequence ID string (first element) and pairwise distances for that sequence (subsequent elements).
    ::param sequenceIdDistances::        List of param eachSequenceId.    
    """
    with open(TpzFilepath) as TpzFile:
        numberSequences = TpzFile.readline()
        numberSequences = numberSequences.strip('\n')
        numberSequences = int(numberSequences)
        sequenceIdDistances = []
        eachSequenceIdDistance = []
        splitLinesClean = []
        for line in TpzFile:
            splitLine = line.split(' ')
            for item in splitLine:
                item = item.strip('\n')
                if item != '':
                    splitLinesClean.append(item)
        while len(sequenceIdDistances) < numberSequences:
            eachSequenceIdDistance = splitLinesClean[:numberSequences+1] 
            sequenceIdDistances.append(eachSequenceIdDistance)
            del splitLinesClean[:numberSequences+1]
    return sequenceIdDistances

def getPed(TpzFilepath):
    """
    Calls the function getTPZiddist and returns dictionary of pairwise sequence (tuple) keys and pairwise maximum likelihood distance values. 
    
    ::param TpzFilepath::      Path to Tree-Puzzle file (input).
    ::param TpzIdDistances::      List of lists of sequence ID string (first element) and distances for that sequence (subsequent elements).
    ::param TpzDistancesPerSequence::  List of param distanceList.
    ::param eachDistanceList::   List of distances for one sequence.
    ::param TPZdistances::   List of dictionary values (float distances for all sequences.
    ::param sequence1Count::      Var to keep count of 1st sequence ID in pair.
    ::param sequence2Count::      Var to keep count of 2nd sequence ID in pair.
    ::param eachSequencePair::   Tuple of sequence string pairs.
    ::param sequencePairs::      List of dictionary keys (param eachseqpair2).
    ::param entries::        List of sequence ID string (first element) and distances for that sequence (subsequent elements).
    ::param countDistances:: Var to keep count of distance values per sequence.
    ::param PEDSinit::       Dictionary of pairwise sequence keys and evolutionary distance values.
    ::param PEDS::           Dictionary of pairwise sequence keys and evolutionary distance values, excluding pairwise comparisons of sequence with itself.
    """
    TpzIdDistances = getTpzIdDist(TpzFilepath)
    TpzDistanceListPerSequence = [] 
    TpzAllDistances = []
    sequencePairs = []
    sequence1Count = 0
    while sequence1Count < len(TpzIdDistances):
        sequence2Count = 0
        while sequence2Count < len(TpzIdDistances):
            eachSequencePair = (TpzIdDistances[sequence1Count][0], TpzIdDistances[sequence2Count][0])
            sequencePairs.append(eachSequencePair)
            sequence2Count += 1
        sequence1Count += 1
    for entries in TpzIdDistances:
        TpzDistanceListPerSequence.append(entries[1:])
    for eachDistanceList in TpzDistanceListPerSequence:
        countDistances = 0
        while countDistances < len(eachDistanceList):
            eachDistanceList[countDistances] = float(eachDistanceList[countDistances])
            countDistances += 1
        TpzAllDistances = TpzAllDistances + eachDistanceList
    
    PedsInit = dict(itertools.izip(sequencePairs, TpzAllDistances))
    Peds = {}
    for key,value in PedsInit.iteritems():
        if key[0] != key[1]:
            Peds[key] = value
    return Peds

def PudXPedY(SdtFilepath, TpzFilepath):
    """
    Calls the functions getPUD and getPED and returns distance values for corresponding keys in both dictionaries, as two lists packed in a tuple.
    
    ::param filepath::  Path to SDT file (input).
    ::param filepath2:: Path to Tree-Puzzle file (input).
    ::param PUDdict::   Dictionary of pairwise sequence keys and uncorrected distance values, excluding pairwise comparisons of sequence with itself.
    ::param PEDdict::   Dictionary of pairwise sequence keys and evolutionary distance values, excluding pairwise comparisons of sequence with itself.
    ::param x::         List of evolutionary distances.
    ::param y::         List of uncorrected distances.
    """
    PudDict = getPud(SdtFilepath)
    PedDict = getPed(TpzFilepath)
    x = []
    y = []
    for key in PudDict.keys():
        y.append(PudDict[key])
        x.append(PedDict[key])
    return (x,y)

def reportDistances(SdtFilepath, TpzFilepath):
    """
    Prints sequence pairs and distances (uncorrected and evolutionary) to file.
    
    ::param filepath::  Path to SDT file (input).
    ::param filepath2:: Path to Tree-Puzzle file (input).
    ::param PUDdict::   Dictionary of pairwise sequence keys and uncorrected distance values, excluding pairwise comparisons of sequence with itself.
    ::param PEDdict::   Dictionary of pairwise sequence keys and evolutionary distance values, excluding pairwise comparisons of sequence with itself.
    ::param report::    Function output written to file.        
    """
    PudDict = getPud(SdtFilepath)
    PedDict = getPed(TpzFilepath)
    with open("PUDPED.txt", "w") as report:
        report.write("sequencePair, uncorrectedDistance, evolutionaryDistance \n")
        for key in PudDict.keys():
            report.write(str(key) + ", " + str(PudDict[key]) + ", " + str(PedDict[key]) + "\n")
    return report

# Output pairwise distances
reportDistances(SdtFilepath, TpzFilepath)
print "Pairwise distances printed to PUDPED.txt"

# Output scatterplot of pairwise distances
scatterPoints = PudXPedY(SdtFilepath, TpzFilepath)
plt.scatter(scatterPoints[0], scatterPoints[1])
plt.xlabel('Pairwise evolutionary distance')
plt.ylabel('Pairwise uncorrected distance')
plt.axis([0, np.ceil(max(scatterPoints[0])), 0, np.ceil(max(scatterPoints[1]))])
plt.plot([0, np.ceil(max(scatterPoints[0]))], [0, np.ceil(max(scatterPoints[1]))], 'k:')
plt.savefig("PUDPEDplot.png")
print "Plot saved as PUDPEDplot.png"