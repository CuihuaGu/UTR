# Cuihua Gu 
# 2015-11-3
# Process LUSC..25.0.000001.pval.txt and LUSC..25.0.000001.drec.txt files, find the shortened sequence

import sys
import numpy
# -----------------VARIABLES
pvalFile=sys.argv[1] # LUSC..25.0.000001.pval.txt
drecFile=sys.argv[2] # LUSC..25.0.000001.drec.txt
pvalCutoff=float(sys.argv[3]) # 0.00001
outpref=sys.argv[4] # LUSC..25.0.000001
print sys.argv
faRef="/export/data/guc1/HCC/RNAseq/CHIFIND/CHIFIND/lib/hg19_bt_index/hg19.fa"

# -----------------FUNCTION
def getSigApaPair(pvalFile,drecFile,pvalCutoff,outpref):
	fp = open(pvalFile,"r")
	fd = open(drecFile,"r")
	fpd = {}
	fdd = {}
	count = 0
	for line in fp.readlines():
		if count ==0:
			namep = line.rstrip("\n").split('\t')
			count = count +1
		else:
			line = line.rstrip("\n").split("\t")
			fpd[line[0]]=line[1:]
	count = 0
	for line in fd.readlines():
		if count ==0:
			named = line.rstrip("\n").split('\t')
			count = count +1
		else:
			line = line.rstrip("\n").split("\t")
			fdd[line[0]]=line[1:]
	fp.close()
	fd.close()
	if len(fpd)!=len(fdd):
		print "error: pval file is not in the same dimention with drec file."
	if namep[:-2] != named:
		print "error: pval name doesn't equal with drec file."
		print namep[:-2]
		print named
	# find number
	#fout_l = open('_'.join([outpref,str(minDynamicNum),str(pvalCutoff),"lengthen.txt"]),"w")
	#fout_s = open('_'.join([outpref,str(minDynamicNum),str(pvalCutoff),"shorten.txt"]),"w")
	fout = open('_'.join([outpref,str(pvalCutoff),"summary.txt"]),"w")
	for key in fpd:
		#print key
		short_sig = 0
		short_notsig = 0
		long_sig = 0
		long_notsig = 0
		short_p = []
		long_p = []
		tmp1 = [float(ele) for ele in fpd[key]]
		for i in range(0,len(fpd[key])-2):
			if fdd[key][i]=="shorten":
				short_p.append(float(fpd[key][i]))
				if tmp1[i]<=pvalCutoff:
					short_sig = short_sig + 1
				else:
					short_notsig = short_notsig + 1
			elif fdd[key][i]=="lengthen":
				long_p.append(float(fpd[key][i]))
				if tmp1[i]<=pvalCutoff:
					long_sig = long_sig + 1
				else:
					long_notsig = long_notsig + 1
		if len(short_p)==0:
			short_p_min='NA'
			short_p_id = 'NA'
			short_p_mean = 'NA'
		else:
			short_p_min = min(short_p)
			short_p_id = namep[tmp1.index(min(short_p))]
			short_p_mean = numpy.mean(short_p)
		if len(long_p)==0:
			long_p_min='NA'
			long_p_id ='NA'
			long_p_mean='NA'
		else:
			long_p_min = min(long_p)
			long_p_id = namep[tmp1.index(min(long_p))]
			long_p_mean = numpy.mean(long_p)
		#print min(short_p)
		#print fpd[key]
		#print [key,str(short_sig),str(short_notsig),str(long_sig),str(long_notsig),str(short_p_mean),str(short_p_min),str(short_p_id),str(long_p_mean),str(long_p_min),str(long_p_id)]
		fout.write('\t'.join([key,str(short_sig),str(short_notsig),str(long_sig),str(long_notsig),str(short_p_mean),str(short_p_min),str(short_p_id),str(long_p_mean),str(long_p_min),str(long_p_id)])+'\n')
	fout.close()
print "running main"
getSigApaPair(pvalFile,drecFile,pvalCutoff,outpref)



