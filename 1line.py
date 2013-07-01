#convert block of text into one line.

import sys

filename = raw_input("Enter file name: ")

infile = open(filename, 'r')
lines = infile.readlines()
infile.close()

for eachline in lines:
	index = lines.index(eachline)
	lines[index] = eachline.replace("\n", "")

outfile = open('1line.txt', 'w')
outfile.write(''.join(lines))
outfile.close()
