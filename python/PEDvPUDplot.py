# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 12:59:47 2015

Make scatterplot of pairwise uncorrected vs evolutionary distances.

Pairwise uncorrected distances are obtained from SDT 1.0 (filename_mat.txt), while pairwise maximum likelihood distances are obtained from Tree-Puzzle 5.2 (outdist.txt). Examples of both input files are in the folder "PEDvPUDplot_exampleInput". Both files are read to extract sequence ID and distance values. 

"""

import numpy as np
import matplotlib.pyplot as plt
import itertools

def getSdtId(sdtFilepath):
    """
    Reads SDT 1.0 matrix output (filename_mat.txt) and returns list of sequence IDs.
    
    ::param sdtFilepath::   Path to SDT file (input).
    ::returns sequenceIds:: List of sequence ID strings.
    """
    with open(sdtFilepath) as sdtFile:
        sequenceIds = []
        eachSequenceId = ""
        for line in sdtFile:
            splitLine = line.split(',')
            eachSequenceId = splitLine[0]
            if len(eachSequenceId) > 10:
                eachSequenceId = eachSequenceId[:10]
                sequenceIds.append(eachSequenceId)
            else:
                sequenceIds.append(eachSequenceId)
    return sequenceIds

def getSdtDist(sdtFilepath):
    """
    Reads SDT 1.0 matrix output (filename_mat.txt) of pairwise identities, and returns list of pairwise distances as float numbers. 
        
    ::param sdtFilepath::    Path to SDT file (input).
    ::returns sdtDistances:: List of pairwise distances (float).
    """
    with open(sdtFilepath) as sdtFile:
        sdtDistances = []
        for line in sdtFile:
            splitLine = line.split(',')
            similarities = splitLine[1:(len(splitLine)-1)]
            for item in similarities:
                item = item.strip('\n')
                item = float(item)
                sdtDistances.append(1-item)
    return sdtDistances    

def getPud(sdtFilepath):
    """
    Calls functions getSdtId and getSdtDist and returns dictionary of pairwise sequence (tuple) keys and pairwise uncorrected distance values.
    
    ::param sdtFilepath:: Path to SDT file (input).
    ::returns puds::      Dictionary of pairwise sequence keys and uncorrected distance values, excluding pairwise comparisons of sequence with itself.
    """
    sdtIds = getSdtId(sdtFilepath)
    sdtDistances = getSdtDist(sdtFilepath) 
    sequencePairs = [] 
    sequence1Count = 1
    while sequence1Count<len(sdtIds):
        sequence2Count = 0
        while sequence2Count < sequence1Count:
            eachSequencePair = (sdtIds[sequence1Count], sdtIds[sequence2Count])
            sequencePairs.append(eachSequencePair)
            sequence2Count += 1
        sequence1Count += 1
    puds = dict(itertools.izip(sequencePairs, sdtDistances))
    return puds

def getTpzListIdAndDist(tpzFilepath):
    """
    Reads Tree-Puzzle 5.2 maximum likelihood distance file (outdist.txt) and returns list of lists comprising sequence id (first element) and all pairwise distances.
    
    ::param tpzFilepath::           Path to Tree-Puzzle file (input).
    ::returns sequenceIdDistances:: List of lists of sequence IDs and associated pairwise distances.    
    """
    with open(tpzFilepath) as tpzFile:
        numberSequences = tpzFile.readline()
        numberSequences = numberSequences.strip('\n')
        numberSequences = int(numberSequences)
        sequenceIdDistances = []
        eachSequenceIdDistance = []
        splitLinesClean = []
        for line in tpzFile:
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

def getPed(tpzFilepath):
    """
    Calls the function getTpzListIdAndDist and returns dictionary of pairwise sequence (tuple) keys and pairwise maximum likelihood distance values. 
    
    ::param tpzFilepath:: Path to Tree-Puzzle file (input).
    ::returns peds::      Dictionary of pairwise sequence keys and evolutionary distance values, excluding pairwise comparisons of sequence with itself.
    """
    tpzIdDistances = getTpzListIdAndDist(tpzFilepath)
    tpzDistanceListPerSequence = [] 
    tpzAllDistances = []
    sequencePairs = []
    sequence1Count = 0
    while sequence1Count < len(tpzIdDistances):
        sequence2Count = 0
        while sequence2Count < len(tpzIdDistances):
            eachSequencePair = (tpzIdDistances[sequence1Count][0], tpzIdDistances[sequence2Count][0])
            sequencePairs.append(eachSequencePair)
            sequence2Count += 1
        sequence1Count += 1
    for entries in tpzIdDistances:
        tpzDistanceListPerSequence.append(entries[1:])
    for eachDistanceList in tpzDistanceListPerSequence:
        countDistances = 0
        while countDistances < len(eachDistanceList):
            eachDistanceList[countDistances] = float(eachDistanceList[countDistances])
            countDistances += 1
        tpzAllDistances = tpzAllDistances + eachDistanceList
    
    pedsInit = dict(itertools.izip(sequencePairs, tpzAllDistances))
    peds = {}
    for key,value in pedsInit.iteritems():
        if key[0] != key[1]:
            peds[key] = value
    return peds

def pedXPudY(sdtFilepath, tpzFilepath):
    """
    Calls the functions getPud and getPed and returns distance values for corresponding keys in both dictionaries, as two lists packed in a tuple.
    
    ::param sdtFilepath:: Path to SDT file (input).
    ::param tpzFilepath:: Path to Tree-Puzzle file (input).
    ::returns (x,y)::     Tuple of x (list of evolutionary distances) and y (list of uncorrected distances).
    """
    pudDict = getPud(sdtFilepath)
    pedDict = getPed(tpzFilepath)
    x = []
    y = []
    for key in pudDict.keys():
        y.append(pudDict[key])
        x.append(pedDict[key])
    return (x,y)

def reportDistances(sdtFilepath, tpzFilepath):
    """
    Prints sequence pairs and distances (uncorrected and evolutionary) to file.
    
    ::param sdtFilepath:: Path to SDT file (input).
    ::param tpzFilepath:: Path to Tree-Puzzle file (input).
    ::returns report::    Writes output (pairwise sequences and associated distances) to file.        
    """
    pudDict = getPud(sdtFilepath)
    pedDict = getPed(tpzFilepath)
    with open("PUDPED.txt", "w") as report:
        report.write("sequencePair, uncorrectedDistance, evolutionaryDistance \n")
        for key in pudDict.keys():
            report.write(str(key) + ", " + str(pudDict[key]) + ", " + str(pedDict[key]) + "\n")
    return report

"""
This is the script entry point.
"""
sdtFilepath = raw_input("Enter path to SDT matrix file: ")
tpzFilepath = raw_input("Enter path to TreePuzzle distance output file: ")

# Output pairwise distances
reportDistances(sdtFilepath, tpzFilepath)
print "Pairwise distances printed to PUDPED.txt"

# Output scatterplot of pairwise distances
scatterPoints = pedXPudY(sdtFilepath, tpzFilepath)
plt.scatter(scatterPoints[0], scatterPoints[1])
plt.xlabel('Pairwise evolutionary distance')
plt.ylabel('Pairwise uncorrected distance')
plt.axis([0, np.ceil(max(scatterPoints[0])), 0, np.ceil(max(scatterPoints[1]))])
plt.plot([0, np.ceil(max(scatterPoints[0]))], [0, np.ceil(max(scatterPoints[1]))], 'k:')
plt.savefig("PUDPEDplot.png")
print "Plot saved as PUDPEDplot.png"