######################################################################################
## July 18, 2013                                                                    ##
## Extract first, second, third, or first+second codon positions from PHYLIP file.  ##
######################################################################################

filename = raw_input("Enter file name: ")
linesIn =  open(filename, 'r').readlines()
(open(filename,'r')).close()

codonIn = raw_input("Codon positions to extract (1, 2, 3, or 1+2): ")

taxalist = []
seqlist = []
extractedseqs = []

###   Start Function Definitions   ###
def getTaxaSeq(inputLines):
	for eachline in inputLines:
		taxseq = eachline.split()
		if taxseq != []:
			taxalist.append(taxseq[0])
			seqlist.append(taxseq[1])
	return None


def run(codon, reportedTaxa, reportedChar):
	if codon == "1" or codon == "2" or codon == "3" or codon == "1+2":
		#Extract nucleotides from codon positions specified by user
		if reportedChar % 3 != 0:
			print "!!! Reported length of sequences is not a multiple of 3!!!"
			return None
		else:
			if reportedTaxa != len(taxalist):
				print "!!! Reported taxa (" + str(reportedTaxa) + \
				") is NOT the same as actual taxa (" + str(len(taxalist)) + ")!!!"
			else:
				j = 0
				while j < reportedTaxa:
					if reportedChar != len(seqlist[j]):
						print "!!! Reported number of sites (" + str(reportedChar) + \
						") is INCORRECT for sequence " + taxalist[j] + "!!!"
						return None
					else:
						startpos = 0
						newseq = ""
						while startpos < reportedChar:
							if codon == "1":
								newseq += seqlist[j][startpos]
							elif codon == "2":
								newseq += seqlist[j][(startpos + 1)]
							elif codon == "3":
								newseq += seqlist[j][(startpos + 2)]
							elif codon == "1+2":
								newseq += seqlist[j][startpos]
								newseq += seqlist[j][(startpos + 1)]
							startpos += 3
						extractedseqs.append(newseq)
						j += 1
	else:
		print "!!! Invalid codon position selected for extraction!!!"
		print "Valid codon positions are 1, 2, 3, or 1+2."
	return None


def writeResults(codonPos):
	if len(extractedseqs) == 0 or len(extractedseqs) != len(taxalist) or \
	(codonPos != "1" and codonPos != "2" and codonPos != "3" and codonPos != "1+2"):
		print "!!! No sequences extracted!!!"
	else:
		outfile = open('extractedCodon'+codonIn+'_'+filename, 'w')
		outfile.write(str(numtaxa) + "  " + str(len(extractedseqs[0])) + "\n")
		k = 0
		while k < numtaxa:
			outfile.write(taxalist[k] + "  " + extractedseqs[k] + "\n")
			k +=1
		outfile.close()
		print "Sequences extracted to file extractedCodon" + codonIn + "_" + filename
	return None

###   End Function Definitions   ###


getTaxaSeq(linesIn)
# Get number of taxa and characters
numtaxa = int(taxalist[0])
numchar = int(seqlist[0])

del taxalist[0]
del seqlist[0]

run(codonIn, numtaxa, numchar)
writeResults(codonIn)

