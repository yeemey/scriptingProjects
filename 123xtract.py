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


def checkLength(reportedTaxa, reportedChar):
	#Check if reported sequence length is a multiple of 3
		if reportedChar % 3 != 0:
			print "Reported length of sequences is " + str(reportedChar)
			print "Not a multiple of 3!"
			return None
		else:	
			if reportedTaxa == len(taxalist):
				print "Reported taxa (" + str(reportedTaxa) + \
				") is the same as actual taxa (" + str(len(taxalist)) + ")."
				#Check length of each sequence against reported lengths
				i = 0
				while i < len(seqlist):
					if len(seqlist[i]) == reportedChar:
						print "Reported number of sites (" + str(reportedChar) + \
						") is correct for sequence " + taxalist[i]
					else:
						print "!!! Reported number of sites (" + str(reportedChar) + \
						") is INCORRECT for sequence " + taxalist[i] + "!!!"
					i += 1				
			else:
				print "!!! Reported taxa (" + str(reportedTaxa) + \
				") is NOT the same as actual taxa (" + str(len(taxalist)) + ")!!!"
			return None


def run(codon, actualTaxa, actualChar):
	if codon == "1" or codon == "2" or codon == "3" or codon == "1+2":
		#Extract nucleotides from codon positions specified by user
		j = 0
		while j < actualTaxa:
			startpos = 0
			newseq = ""
			while startpos < actualChar:
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
		print "Valid codon positions are 1, 2, 3, or 1+2. Exiting program..."
	return None

###   End Function Definitions   ###


getTaxaSeq(linesIn)
# Get number of taxa and characters
numtaxa = int(taxalist[0])
numchar = int(seqlist[0])
del taxalist[0]
del seqlist[0]

checkLength(numtaxa, numchar)
run(codonIn, numtaxa, numchar)

outfile = open('extractedCodon'+codonIn+'_'+filename, 'w')
outfile.write(str(numtaxa) + "  " + str(len(extractedseqs[0])) + "\n")
k = 0
while k < numtaxa:
	outfile.write(taxalist[k] + "  " + extractedseqs[k] + "\n")
	k +=1
outfile.close()