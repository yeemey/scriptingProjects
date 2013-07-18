############################################################################
## July 18, 2013                                                          ##
## Extract first+second, or only third codon positions from PHYLIP file.  ##
############################################################################


filename = raw_input("Enter file name: ")
infile = open(filename, 'r')
lines = infile.readlines()
infile.close()

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

#Check if reported sequence length is a multiple of 3
if numchar % 3 != 0:
	print "Reported length of sequences is " + str(numchar)
	print "Not a multiple of 3!"
else:
	print "Reported length of sequences is " + str(numchar)
	print "There should be " + str(numchar/3) + " codons in this alignment."

ziptaxseq = zip(taxalist, seqlist)

i = 1
while i < len(ziptaxseq):
	if len(ziptaxseq[i][1]) == numchar:
		print "Length of " + ziptaxseq[i][0] + " is " + str(len(ziptaxseq[i][1])) + ", as reported."
	else:
		print "Length of sequence " + ziptaxseq[i][0] + " is " + str(len(ziptaxseq[i][1])) + ", not " + str(numchar) + "!!"
		print "Check sequence length in file. Exiting program..."
		break
	i += 1


