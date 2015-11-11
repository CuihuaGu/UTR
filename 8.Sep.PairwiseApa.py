# 2015-9-8
# Cuihua Gu
# Pair assigned APA peaks to proximal and distal pairs
# INPUT
	# 1. /data/guc1/UTR/Output/1_ApaAnnotate/AssignedAPA_refseq.txt
	# 2. output name
# OUTPUT
	# 1. Result
	# 2. Log file
# USAGE
	# python script.py AssignedAPA_refseq.txt PairwiseApa

import sys
import time
import datetime
import copy
import re

def pairTwoPeaks(li1raw,li2raw):
	region = []
	if int(li1raw[1]) < int(li2raw[1]):
		li1 = li1raw
		li2 = li2raw
	else:
		li2 = li1raw
		li1 = li2raw
	li1_1 = int(li1[6])
	li1_2 = int(li1[7])
	li2_1 = int(li2[6])
	li2_2 = int(li2[7])
	pk1 = int(li1[1])
	pk2 = int(li2[1])
	if li1[5] == li2[5]:
		r_ = determRegionOverl(li1_1,li1_2,li2_1,li2_2)
		if r_ == "tandem":
			region = determInSameRegion(li1_1,li1_2,li2_1,li2_2,pk1,pk2,li1[3])
		elif r_ == "mutually_exclusive":
			region = determInDiffRegion(li1_1,li1_2,li2_1,li2_2,pk1,pk2,li1[3])
	elif sorted(list([li1[5],li2[5]]))==['ExtR,', 'UR,']:
		region = determInSameRegion(li1_1,li1_2,li2_1,li2_2,pk1,pk2,li1[3])
	else:
		region = determInDiffRegion(li1_1,li1_2,li2_1,li2_2,pk1,pk2,li1[3])
	return region, determRegionOverl(li1_1,li1_2,li2_1,li2_2)

# a1 < a2, b1 < b2, p1 < p2
def determInSameRegion(a1,a2,b1,b2,p1,p2,strand):
	if strand == "+":
		start1 = min(a1,b1)
		end1 = p1
		start2 = p1
		end2 = p2
	else:
		start1 = p1
		end1 = p2
		start2 = p2
		end2 = max(a2,b2)
	return [start1,end1,start2,end2]

def determInDiffRegion(a1,a2,b1,b2,p1,p2,strand):
	if strand == "+":
		start1 = a1
		end1 = p1
		start2 = b1
		end2 = p2
	else:
		start1 = p1
		end1 = a2
		start2 = p2
		end2 = b2
	return [start1,end1,start2,end2]

def determRegionOverl(a1,a2,b1,b2):
	if a1<=b1 and a2>=b1:
		r = "tandem"
	elif a2>=b2 and a1<=b2:
		r = "tandem"
	elif a2<b1 and a1<b1:
		r = "mutually_exclusive"
	elif a2>b2 and a1>b2:
		r = "mutually_exclusive"
	else:
		print "error"
	return r

def apaToBed(apaFile):
        f = open(apaFile,"r")
        bed = []
        for line in f.readlines():
                line = line.rstrip("\n").split("\t")
                name_p = line[0]+"_"+line[4]+"_pro"
                name_d = line[0]+"_"+line[4]+"_dis"
                chrom = line[1]
                strand = line[2]
                startp = line[5]
                endp = line[6]
                startd = line[7]
                endd = line[8]
                typ = line[9]
                if strand=="+":
                        tmp_p = [chrom,startp,endp,name_p,0,strand]
                        tmp_d = [chrom,startd,endd,name_d,0,strand]
                        bed.append(tmp_p)
                        bed.append(tmp_d)
                else:
                        tmp_p = [chrom,startd,endd,name_p,0,strand]
                        tmp_d = [chrom,startp,endp,name_d,0,strand]
                        bed.append(tmp_p)
                        bed.append(tmp_d)
        f.close()
        return bed

def writeBed(bed,outPref):
	f = open(outPref,'w')
	for i in range(0,len(bed)):
		tmp = [str(elem) for elem in bed[i]]
		f.write("\t".join(tmp)+"\n")
	f.close

# --------------------1. Read AssignedAPA_refseq.txt

logFile = open(sys.argv[2]+".log","w")
ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
logFile.write(ts+'--------------------Start\n')
f = open(sys.argv[1],"r")
AssdApa = {}
for line in f.readlines():
	line = line.rstrip("\n").split("\t")
	if AssdApa.has_key(line[0]):
		AssdApa[line[0]].append(line[1:9])
	else:
		AssdApa[line[0]] = []
		AssdApa[line[0]].append(line[1:9])

# --------------------2. Remove sigle APA entry
Result = []
for key in AssdApa.keys():
	if len(AssdApa[key]) == 1:
		pass
	else:
		chroms = AssdApa[key][0][2]
		stran = AssdApa[key][0][3]
		typ = AssdApa[key][0][4]
		l = len(AssdApa[key])
		for i in range(0,l-1):
			for j in range(i+1,l):
				pks = ",".join([AssdApa[key][i][0],AssdApa[key][j][0]])
				strin = re.compile('chr.{1,2}_peaks')
				pks = strin.sub('',pks)
				resraw = [key, chroms, stran, typ, pks]
				#print AssdApa[key][i],AssdApa[key][i]
				res,r_r = pairTwoPeaks(AssdApa[key][i],AssdApa[key][j])
				Result.append(resraw+res+[r_r])

f = open(sys.argv[2]+".txt","w")
for elem in Result:
	elem=[str(tmp) for tmp in elem]
	f.write('\t'.join(elem))
	f.write('\n')
f.close()

bed = apaToBed(sys.argv[2]+".txt")
writeBed(bed,sys.argv[2]+".bed")

ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
logFile.write(ts+'--------------------End\n')
logFile.close()
