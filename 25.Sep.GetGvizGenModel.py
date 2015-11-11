# 2015-9-5
# Cuihua Gu
# get refSeq gene annotation data frame from refseq annotation file
# INPUT
	# 1. refSeq annotation file
	# 2. output prefix
# OUTPUT
	# 2. table file
# USAGE
#chromosome\tstart\tend\twidth\tstrand\tfeature\tgene\texon\ttranscript\tsymbol
#chr7 26591441 26591829   389      + lincRNA ENSG00000233760	ENSE00001693369 ENST00000420912 AC004947.2
import sys

inputFile = sys.argv[1] # /data/guc1/UTR/Data/Annotation/RefSeqGenesAnno.txt
outpref = sys.argv[2]

def getModel(line):
	res = []
	chrom = line[2]
	symbol = line[12]
	gene = line[12]
	transcript = line[1]
	strand = line[3]
	exonStart = line[9].rstrip(',').split(',')
	exonEnd = line[10].rstrip(',').split(',')
	for i in range(0,len(exonStart)):
		l = exonStart[i]
		r = exonEnd[i]
		res.append([chrom,str(l),str(r),str(int(r)-int(l)),strand,"exon",gene,''.join(["exon",str(i+1)]),transcript,symbol])
	return res


f = open(inputFile,"r")
count=1
Result = []
for line in f.readlines():
	if count==1:
		pass
	else:
		line = line.rstrip("\n").split("\t")
		Result= Result + getModel(line)
	count = count + 1
	if count==1000:
		print "1000 lines done"
	if count==10000:
		print "10000 lines done"
	if count==20000:
		print "20000 lines done"

f.close()

f = open(outpref,"w")
f.write("chromosome\tstart\tend\twidth\tstrand\tfeature\tgene\texon\ttranscript\tsymbol\n")
for elem in Result:
	f.write('\t'.join(elem))
	f.write('\n')
f.close()