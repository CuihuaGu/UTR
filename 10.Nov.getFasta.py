import os 
import sys
pvalCutoff="1e-06" # 0.00001
outpref=sys.argv[1] # LUSC..25.0.000001
fc = 2
def getSigApaSequence(outpref,pvalCutoff,fc):
	apaPair = "/data/guc1/UTR/Output/2_ApaPair/PairwiseApa_rmOtherThanLastExon_less20kb.bed"
	shortFile = '_'.join([outpref,str(pvalCutoff),str(fc)+"fold_short.txt"])
	longFile = '_'.join([outpref,str(pvalCutoff),str(fc)+"fold_long.txt"])
	shortBed = '_'.join([outpref,str(pvalCutoff),str(fc)+"fold_short.bed"])
	longBed = '_'.join([outpref,str(pvalCutoff),str(fc)+"fold_long.bed"])
	shortFa = '_'.join([outpref,str(pvalCutoff),str(fc)+"fold_short.fa"])
	longFa = '_'.join([outpref,str(pvalCutoff),str(fc)+"fold_long.fa"])
	f = open(apaPair,"r")
	d = {}
	s = []
	l = []
	for line in f.readlines():
		line = line.rstrip("\n").split("\t")
		d[line[3]]=line
	f.close()
	f = open(shortFile,"r")
	for line in f.readlines():
		line = line.rstrip("\n").split("\t")
		s.append(line[0]+"_dis")
	f.close()
	f = open(longFile,"r")
	for line in f.readlines():
		line = line.rstrip("\n").split("\t")
		l.append(line[0]+"_dis")
	f.close()
	f = open(shortBed,"w")
	for ele in s:
		f.write('\t'.join(d[ele])+'\n')
	f.close()
	f = open(longBed,"w")
	for ele in l:
		f.write('\t'.join(d[ele])+'\n')
	f.close()
	cmd1 = "bedtools getfasta -fi /export/data/guc1/HCC/RNAseq/CHIFIND/CHIFIND/lib/hg19_bt_index/hg19.fa -name -bed "+shortBed+" -s -fo "+shortFa
	cmd2 = "bedtools getfasta -fi /export/data/guc1/HCC/RNAseq/CHIFIND/CHIFIND/lib/hg19_bt_index/hg19.fa -name -bed "+longBed+" -s -fo "+longFa
	os.system(cmd1)
	os.system(cmd2)

getSigApaSequence(outpref,pvalCutoff,fc)
	
