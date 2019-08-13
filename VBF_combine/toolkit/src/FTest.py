import ROOT
from ROOT import *
import sys

#name of general function class
fcname = "Pol"
#how many orders there are
norders = 5
#what order to start at, eg 2 for Pol, 1 for expPol
storder = 2

#beginning of names of input files
filenames = ["ndof", "chi2", "rech", "prob"]
#make separate table for each category
for cat in range(9):
	outfile = open("FTable_" + str(cat) + ".txt", "w")
	outfile.write("FTest Table for Category " + str(cat) + "\n")
	outfile.write("%s order \tndof \t \t \tchi2 \t \t \tchi2/ndof \t \tprob \t \t \tFprob \n"%fcname)
	infiles = []
	#now print table row for each order
	for name in filenames:
		infiles.append( open(name + "_" + str(cat) + ".txt", "r"))
	#2d array to hold all stats for this category.
	stats = []
	stats.append([0.0 for i in range(norders)]) #ndof
	stats.append([0.0 for i in range(norders)]) #chi2
	stats.append([0.0 for i in range(norders)]) #reduced chi2
	stats.append([0.0 for i in range(norders)]) #prob
	stats.append([0.0 for i in range(norders)]) #fval
	#read from each input file
	for i in range(4):
		j = 0
		for line in infiles[i]:
			if j > norders: sys.exit("Error! More than %d entries in %s_%d.txt" %(norders,filenames[i],cat))
			stats[i][j] = float(line)
			j += 1
	for i in range(norders):
		if i == 0:
			stats[4][i] = float("nan")
		else:
			stats[4][i] = TMath.Prob( (stats[1][i-1] - stats[1][i]), 1 )
		outfile.write(str(i+storder) + "\t \t")
		for j in range(5):
			outfile.write("%f \t \t" %stats[j][i])
		outfile.write("\n")
	outfile.close()
	for fil in infiles:
		fil.close()
