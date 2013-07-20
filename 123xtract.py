############################################################################
## July 18, 2013                                                          ##
## Extract first+second, or only third codon positions from PHYLIP file.  ##
############################################################################

filename = raw_input("Enter file name: ")
infile = open(filename, 'r')
lines = infile.readlines()
infile.close()

extractpos = raw_input("Codon positions to extract (1, 2, 3, or 1+2): ")

if extractpos is "1" or extractpos is "2" or extractpos is "3" or extractpos is "1+2":
#Unpack all lines
	taxalist = []
	seqlist = []
	for eachline in lines:
		taxseq = eachline.split()
		if taxseq != []:
			taxalist.append(taxseq[0])
			seqlist.append(taxseq[1])

#Unpack number of taxa and characters from the first line, change type to integer
	numtaxa = int(taxalist[0])
	numchar = int(seqlist[0])
	del taxalist[0]
	del seqlist[0]
	#print taxalist[0]

#Check number of sequences
	if numtaxa == len(taxalist):
#Check if reported sequence length is a multiple of 3
		if numchar % 3 == 0:
			print "Reported length of sequences is " + str(numchar)
			print "There should be " + str(numchar/3) + " codons in this alignment."
#Zip taxa and sequence together in a list of tuples
			ziptaxseq = zip(taxalist, seqlist)
#Check length of each sequence against reported lengths
			i = 0
			while i < len(ziptaxseq):
				if len(ziptaxseq[i][1]) == numchar:
					print "Length of " + ziptaxseq[i][0] + " is " + str(len(ziptaxseq[i][1])) + ", as reported."
				else:
					print "Length of sequence " + ziptaxseq[i][0] + " is " + str(len(ziptaxseq[i][1])) + ", not " + str(numchar) + "!!"
					print "Check sequence length in file. Exiting program..."
				i += 1
		else:
			print "Reported length of sequences is " + str(numchar)
			print "Not a multiple of 3! Exiting program..."
		
#Extract nucleotides from codon positions specified by user
		extractedseqs = []
		j = 0
		while j < numtaxa:
			startpos = 0
			newseq = ""
			while startpos < numchar:
				if extractpos == "1":
					newseq += seqlist[j][startpos]
				elif extractpos == "2":
					newseq += seqlist[j][(startpos + 1)]
				elif extractpos == "3":
					newseq += seqlist[j][(startpos + 2)]
				elif extractpos == "1+2":
					newseq += seqlist[j][startpos]
					newseq += seqlist[j][(startpos + 1)]
				startpos += 3
			extractedseqs.append(newseq)
			j += 1

		ziptaxnewseq = zip(taxalist, extractedseqs)
		print ziptaxnewseq
	else:
		print "Reported number of taxa is not the same as number of sequences in file!"
		print "Exiting program..."
else:
	print "Invalid codon position selection for extraction!"
	print "Valid codon positions are 1, 2, 3, or 1+2. Exiting program..."



		

		
		
