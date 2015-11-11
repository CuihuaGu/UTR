# 2015-10-20
# Cuihua Gu
# Remove apa locates before last exon
# INPUT
	# 1. apaFile="/data/guc1/UTR/Output/1_ApaAnnotate/AssignedAPA_refseq.txt"
# OUTPUT
	# 2. filtered apaFile
# USAGE

import sys
#---------------------USER INPUT
apaFile = sys.argv[1]
outFile = sys.argv[2]

#---------------------FUNCTION
def filterApa(apaFile,outFile):
	f = open(apaFile,"r")
	fout = open(outFile,"w")
	for line in f.readlines():
		line2 = line.rstrip("\n").split("\t")
		seq = line2[6].split(",")[1]
		if seq == "":
			fout.write(line)
		elif int(seq) == 1:
			fout.write(line)
		else:
			pass
	f.close()
	fout.close()

#---------------------MAIN
filterApa(apaFile,outFile)


