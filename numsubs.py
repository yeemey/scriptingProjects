# script to count number of changes
# input seq_pars.log file
# first get outgroup name (str)
# search for <nt outgroup> and replace with <*****nt outgroup>
# save changes in input file
# now count number of "nt ==> nt" and "nt --> nt"
# print output number of nt changes to numsubs.txt
import sys
print sys.argv
import re

filename = raw_input("Enter log file name: ")
outgroup = raw_input("Enter outgroup: ")


infile = open(filename, 'r')
lines = infile.readlines()
infile.close()

#Only consider the change list at the bottom of the log file
chglist = []
chgliststart = lines.index("Character change lists:\n")
while chgliststart < len(lines):
	chglist.append(lines[chgliststart])
	chgliststart += 1

nuc = [' A ', ' C ', ' G ', ' T ']

for eachline in chglist:
	for eachnuc in nuc:
#Identify all changes to outgroup and mark with asterisks
		matchout = re.search(eachnuc + outgroup, eachline)
		if matchout:
			index = chglist.index(eachline)
			chglist[index] = eachline.replace(eachnuc, " ****" + eachnuc)
			print chglist[index]

'''
***************Mark's Version***************
'''
def printResults(type, dict, order):
	print "\nPrinting " + type + " counts..."
	for k in order:
		print type + "["+ k + "]: " + str(dict[k]) 
	

# Create
A = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
C = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
G = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
T = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
e = re.compile("([ATCG]\s*==>\s*[ATCG]|[ATCG]\s*-->\s*[ATCG])")
for eachline in chglist:
	matchin = e.findall(eachline)
	if matchin:
		matchStr = ''.join(str(x) for x in matchin)
		print "Reading entry: " + matchStr
		key = matchStr[len(matchStr) - 1].strip()
		print "Which bin are we counting in? " + key
		if matchStr[0] == "A":
			A[key] += 1
		elif matchStr[0] == "C":
			C[key] += 1
		elif matchStr[0] == "G":
			G[key] += 1
		elif matchStr[0] == "T":
			T[key] += 1
order = ["A", "C", "G", "T"]
printResults("A", A, order)
printResults("C", C, order)
printResults("G", G, order)
printResults("T", T, order)

# print "\nPrinting A counts..."
# for k in A.keys():
# 	print "A["+ k + "]: " + str(A[k])
# 	
# print "\nPrinting C counts..."
# for k in C.keys():
# 	print "C["+ k + "]: " + str(C[k])
# 
# print "\nPrinting G counts..."	
# for k in G.keys():
# 	print "G["+ k + "]: " + str(G[k])
# 	
# print "\nPrinting T counts..."
# for k in T.keys():
# 	print "T["+ k + "]: " + str(T[k])
	

'''
END Mark's Version
'''

#Count number of ingroup changes

# Achg = {'AtoC': 0, 'AtoG': 0, 'AtoT': 0}
# Cchg = {'CtoA': 0, 'CtoG': 0, 'CtoT': 0}
# Gchg = {'GtoA': 0, 'GtoC': 0, 'GtoT': 0}
# Tchg = {'TtoA': 0, 'TtoC': 0, 'TtoG': 0}

# for eachline in chglist:
# 	for eachnuc in nuc:
# 		i=0
# 		while i < len(nuc):
# 			matchin = re.search(eachnuc + ("==>" or "-->") + nuc[i], eachline)
# 			if matchin:
# 				if eachnuc == "A ":
# 					if nuc[i] == "C ":
# 						Achg['AtoC'] += 1
# 					elif nuc[i] == "G ":
# 						Achg['AtoG'] += 1
# 					elif nuc[i] == "T ":
# 						Achg['AtoT'] += 1
# 				elif eachnuc == "C ":
# 					if nuc[i] == "A ":
# 						Cchg['CtoA'] += 1
# 					elif nuc[i] == "G ":
# 						Cchg['CtoG'] += 1
# 					elif nuc[i] == "T ":
# 						Cchg['CtoT'] += 1
# 				elif eachnuc == "G ":
# 					if nuc[i] == "A ":
# 						Gchg['GtoA'] += 1
# 					elif nuc[i] == "C ":
# 						Gchg['GtoC'] += 1
# 					elif nuc[i] == "T ":
# 						Gchg['GtoT'] += 1
# 				elif eachnuc == "T ":
# 					if nuc[i] == "A ":
# 						Tchg['TtoA'] += 1
# 					elif nuc[i] == "C ":
# 						Tchg['TtoC'] += 1
# 					elif nuc[i] == "C ":
# 						Tchg['TtoG'] += 1
# 			
# 			i += 1