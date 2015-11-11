# 2015-9-23
# Cuihua Gu
# Assign APA peaks to refSeq transcripts.
# INPUT
	# 1. /data/guc1/UTR/Data/Annotation/RefSeqGenesAnno.txt
	# 2.1. /data/guc1/UTR/Data/ApaRaw/unified-atlas.minus_pks.bed
	# 2.2. /data/guc1/UTR/Data/ApaRaw/unified-atlas.plus_pks.bed
# OUTPUT
	# 1. Result
	# 2. Log file
# USAGE
	# python script.py 

import sys
import time
import datetime
import copy

# --------------------1. Read annotation
def buildModel(annoFile, ext_len):
	ext_len = int(ext_len)
	f = open(annoFile,"r")
	Anno={}
	for line in f.readlines():
		line = line.rstrip("\n").split("\t")
		if Anno.has_key(line[2]):
			Anno[line[2]]['-'.join([line[1],line[4]])] = [line[2],line[3],line[4],line[5],line[8],line[9].rstrip(',').split(','),line[10].rstrip(',').split(','),line[12],line[6],line[7]]
			# chr, strand, txstart, txend, exon_number, exon_start, exon_end, gene_name
		else:
			Anno[line[2]] = {}
			Anno[line[2]]['-'.join([line[1],line[4]])] = [line[2],line[3],line[4],line[5],line[8],line[9].rstrip(',').split(','),line[10].rstrip(',').split(','),line[12],line[6],line[7]]
			# chr, strand, txstart,txend, exon_number, exon_start, exon_end, gene_name, cds_end
	f.close()
	# ----------1.1 Convert annotation to 5kb extended model
	for key in Anno.keys():
		plus = 0
		minus = 0
		plus_ext = 0
		plus_gen = 0
		minus_ext = 0
		minus_gen = 0
		txStart = []
		txEnd = []
		txStart2 = []
		txEnd2 = []
		exonStart = [Anno[key][key2][5] for key2 in Anno[key].keys()]
		exonEnd = [Anno[key][key2][6] for key2 in Anno[key].keys()]
		for k in range(0,len(exonEnd)):
			for k2 in range(0,len(exonEnd[k])):
				txStart.append(exonStart[k][k2])
				txEnd.append(exonEnd[k][k2])
		l = len(txStart)
		for key2 in Anno[key].keys():
			#print key2
			txStart2 = []
			txEnd2 = []
			Strand = copy.copy(Anno[key][key2][1])
			End = copy.copy(Anno[key][key2][6])
			Start = copy.copy(Anno[key][key2][5])
			CdsEnd = copy.copy(Anno[key][key2][9])
			CdsStart = copy.copy(Anno[key][key2][8])
			Anno2 = copy.copy(Anno)
			#print txStart
			if Strand == '+':
				plus = plus + 1
				ext = int(End[-1]) + ext_len
				for k3 in range(0,l):
					#print k3
					if int(txStart[k3]) > int(End[-1]):
						#print txStart[k3], End[-1]
						txStart2.append(txStart[k3])
				#print txStart2
				l2 = len(txStart2)
				#print l2
				temp = [End[-1] for i in range(0,l2)]
				if l2==0:
					nextGene=int(End[-1]) + ext_len + 1
				else:
					nextGene = int(End[-1]) + min([abs(int(txStart2[j]) - int(temp[j])) for j in range(0,l2)])
				if ext > nextGene:
					plus_gen = plus_gen + 1
					if CdsEnd == CdsStart:
						Anno2[key][key2][5].append(Start[-1])
						Anno2[key][key2][6].append(End[-1])
						Anno2[key][key2][5].append(End[-1])
						Anno2[key][key2][6].append(nextGene)
						del Anno2[key][key2][5][-3]
						del Anno2[key][key2][6][-3]
					else:
						Anno2[key][key2][5].append(Start[-1])
						Anno2[key][key2][6].append(CdsEnd)
						Anno2[key][key2][5].append(CdsEnd)
						Anno2[key][key2][6].append(End[-1])
						Anno2[key][key2][5].append(End[-1])
						Anno2[key][key2][6].append(nextGene)
						del Anno2[key][key2][5][-4]
						del Anno2[key][key2][6][-4]
					#print 'plusgen'
				else:
					plus_ext = plus_ext + 1
					if CdsEnd == CdsStart:
						Anno2[key][key2][5].append(Start[-1])
						Anno2[key][key2][6].append(End[-1])
						Anno2[key][key2][5].append(End[-1])
						Anno2[key][key2][6].append(ext)
						del Anno2[key][key2][5][-3]
						del Anno2[key][key2][6][-3]
					else:
						Anno2[key][key2][5].append(Start[-1])
						Anno2[key][key2][6].append(CdsEnd)
						Anno2[key][key2][5].append(CdsEnd)
						Anno2[key][key2][6].append(End[-1])
						Anno2[key][key2][5].append(End[-1])
						Anno2[key][key2][6].append(ext)
						del Anno2[key][key2][5][-4]
						del Anno2[key][key2][6][-4]
				Anno2[key][key2][5][-1] = Anno2[key][key2][5][-2]
			elif Strand == '-':
				minus = minus + 1 
				ext = int(Start[0]) - ext_len
				for k4 in range(0,l):
					if int(txEnd[k4]) < int(Start[0]):
						txEnd2.append(txEnd[k4])
				l2 = len(txEnd2)
				temp2 = [Start[0] for i in range(0,l2)]
				if l2==0:
					nextGene=0
				else:
					nextGene = int(Start[0]) - min([abs(int(txEnd2[j]) - int(temp2[j])) for j in range(0,l2)])
				if ext > nextGene:
					minus_ext = minus_ext + 1
					if CdsEnd == CdsStart:
						Anno2[key][key2][5].insert(0,ext)
						Anno2[key][key2][6].insert(0,Start[0])
						Anno2[key][key2][5].insert(1,Start[0])
						Anno2[key][key2][6].insert(2,End[0])
						del Anno2[key][key2][5][2]
						del Anno2[key][key2][6][2]
					else:
						Anno2[key][key2][5].insert(0,ext)
						Anno2[key][key2][6].insert(0,Start[0])
						Anno2[key][key2][5].insert(1,Start[0])
						Anno2[key][key2][6].insert(1,CdsStart)
						Anno2[key][key2][5].insert(2,CdsStart)
						Anno2[key][key2][6].insert(2,End[0])
						del Anno[key][key2][5][3]
						del Anno[key][key2][6][3]
					#print 'mlusext'
				else:
					minus_gen = minus_gen + 1
					if CdsEnd == CdsStart:
						Anno2[key][key2][5].insert(0,nextGene)
						Anno2[key][key2][6].insert(0,Start[0])
						Anno2[key][key2][5].insert(1,Start[0])
						Anno2[key][key2][6].insert(2,End[0])
						del Anno2[key][key2][5][2]
						del Anno2[key][key2][6][2]
					else:
						Anno2[key][key2][5].insert(0,nextGene)
						Anno2[key][key2][6].insert(0,Start[0])
						Anno2[key][key2][5].insert(1,Start[0])
						Anno2[key][key2][6].insert(1,CdsStart)
						Anno2[key][key2][5].insert(2,CdsStart)
						Anno2[key][key2][6].insert(2,End[0])
						del Anno2[key][key2][5][3]
						del Anno2[key][key2][6][3]
				Anno2[key][key2][6][0] = Anno2[key][key2][6][1]
	buildModel_log =[str(plus), str(minus), str(plus_gen), str(plus_ext), str(minus_gen), str(minus_ext)]
	return Anno2, buildModel_log

# --------------------2. Read APA peaks
def readPeaks(pksFile):
	Apa = []
	f = open(pksFile,"r")
	for line in f.readlines():
		line = line.rstrip("\n").split("\t")
		Apa.append(line)
	f.close()
	return Apa

# --------------------3. Fit Apa into Annotation model
def judgeIfIn(li_a,li_b,numb,strand):
	numb = int(numb)
	if strand == '+':
		for i in range(0,len(li_a)):
			if numb > int(li_a[i]) and numb < int(li_b[i]):
				target = i+1
				break
			else:target = "Na"
	else:
		for i in list(reversed(range(0,len(li_a)))):
			if numb > int(li_a[i]) and numb < int(li_b[i]):
				target = len(li_a) - i
				break
			else:target = "Na"
	return target, [li_a[i],li_b[i]]

def assignApaToRefseq(Anno, Apa):
	Result = []
	for i in range(0,len(Apa)):
		for key in Anno.keys():
			if key == Apa[i][0]:
				for key2 in Anno[key].keys():
					Res = [Apa[i][3],Apa[i][0],Apa[i][5],int(Apa[i][6])]
					if Anno[key][key2][1] == Apa[i][5]:
						target, boundary_r = judgeIfIn(Anno[key][key2][5],Anno[key][key2][6],Apa[i][6],Apa[i][5])
						# print boundary_r
						if target == "Na":
							pass
						else:
							Res = Res + [key2,int(Anno[key][key2][4]),int(boundary_r[0]),int(boundary_r[1]),Anno[key][key2][7],len(Anno[key][key2][5]),target]
							Result.append(Res)
							Res = []
	return Result					

# --------------------4. Convert Data structure
def convertAssdAPA(AssdApa):
	for i in range(0,len(AssdApa)):
		exonNum = AssdApa[i][5]
		extExonNum = AssdApa[i][9]
		assdNum = AssdApa[i][10]
		if extExonNum - exonNum == 2:
			AssdApa[i].append("Coding")
			if assdNum == exonNum + 2:
				AssdApa[i].append("ExtR")
			elif assdNum == exonNum + 1:
				AssdApa[i].append("UR")
			elif assdNum <= exonNum:
				AssdApa[i].append("ExonR")
				AssdApa[i].append(str(exonNum - assdNum + 1))
			else:
				AssdApa[i].append("Unknown")
		elif extExonNum - exonNum == 1:
			AssdApa[i].append("Noncoding")
			if assdNum == exonNum + 1:
				AssdApa[i].append("ExtR")
			elif assdNum <= exonNum:
				AssdApa[i].append("ExonR")
				AssdApa[i].append(str(exonNum - assdNum + 1))
			else:
				AssdApa[i].append("Unknown")
		else:
			AssdApa[i].append("Unknown")
			AssdApa[i].append("Unknown")
	return AssdApa

def mergingCov(AssdApa):
	D = {}
	for i in range(0,len(AssdApa)):
		key1 = AssdApa[i][8]
		key2 = AssdApa[i][0]
		exonStart = AssdApa[i][6]
		exonEnd = AssdApa[i][7]
		Cod = AssdApa[i][11]
		Loci = AssdApa[i][12]
		LociNum = AssdApa[i][-1]
		if D.has_key(key1):
			if D[key1].has_key(key2):
				D[key1][key2].append([AssdApa[i][1],AssdApa[i][2],AssdApa[i][3],exonStart,exonEnd,Cod,Loci,LociNum])
			else:
				D[key1][key2] = [[AssdApa[i][1],AssdApa[i][2],AssdApa[i][3],exonStart,exonEnd,Cod,Loci,LociNum]]
		else:
			D[key1] = {}
			D[key1][key2] = [[AssdApa[i][1],AssdApa[i][2],AssdApa[i][3],exonStart,exonEnd,Cod,Loci,LociNum]]
	return D

# --------------------5. Merging transcript to gene level
def removeNoncoding(D):
	for key1 in D.keys():
		for key2 in D[key1].keys():
			Cod = []
			Loci = []
			temp = []
			temp2 = []
			# Removing Noncoding
			for i in range(0,len(D[key1][key2])):
				Cod.append(D[key1][key2][i][5])
			Cod = sorted(list(set(Cod)))
			if len(Cod) == 1 and Cod[0] == 'Noncoding':
				pass
			elif len(Cod) == 1 and Cod[0] == 'Coding':
				pass
			elif len(Cod) == 2:
				for i2 in range(0,len(D[key1][key2])):
					if D[key1][key2][i2][5] == 'Coding':
						temp.append(D[key1][key2][i2])
				D[key1][key2] = temp
			# removing Exon/Extend is UTR exists
			for i3 in range(0,len(D[key1][key2])):
				Loci.append(D[key1][key2][i3][6])
			Loci = sorted(list(set(Loci)))
			if 'UR' in Loci:
				for i4 in range(0,len(D[key1][key2])):
					if D[key1][key2][i4][6] == 'UR':
						temp.append(D[key1][key2][i4])
				D[key1][key2] = temp
			elif Loci == ['ExonR']:
				pass
			elif Loci == ['ExtR']:
				pass
			elif 'UR' not in Loci and 'ExonR' in Loci and 'ExtR' in Loci:
				print 'ExonR and ExtR exists at the same time!!'
			else:
				print "Unknwon situation!!"
	return D

def mergingToGene(D):
	### merging to gene level
	Result = []
	for key1 in D.keys():
		for key2 in D[key1].keys():
			geneName = key1
			peakName = key2
			chrom = D[key1][key2][0][0]
			strand = D[key1][key2][0][1]
			peak = D[key1][key2][0][2]
			Cod = D[key1][key2][0][5]
			Typ = D[key1][key2][0][6]
			Start_ = []
			End_ = []
			Typ_num = []
			for i in range(0,len(D[key1][key2])):
				Start_.append(D[key1][key2][i][3])
				End_.append(D[key1][key2][i][4])
				if Typ == 'ExonR':
					Typ_num.append(D[key1][key2][i][-1])
				else:
					pass
			Res = [geneName,peakName,peak,chrom,strand,Cod,','.join([Typ,','.join(sorted(list(set(Typ_num))))]),min(Start_),max(End_)]
			Result.append(Res)
	return Result

# -------------------- MAIN
# Writing log function
logFile = open(sys.argv[3]+".log","w")
# running start time
ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
logFile.write(' '.join(sys.argv[0:])+'\n')

logFile.write(ts+'----------------------------------------Start\n')

ext_len = sys.argv[4]
Anno, buildModel_log = buildModel(sys.argv[1], ext_len)  

ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')

logFile.write('\t'.join(["plus", "minus", "plus_gen", "plus_ext", "minus_gen", "minus_ext"])+'\n')
logFile.write('\t'.join(buildModel_log)+'\n')
logFile.write(ts+'--------------------buildModel done\n')

Apa = readPeaks(sys.argv[2])
AssdApa = assignApaToRefseq(Anno, Apa)

AssdApa2 = convertAssdAPA(AssdApa)
D = mergingCov(AssdApa2)
D = removeNoncoding(D)
D = mergingToGene(D)

ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
logFile.write(ts+'--------------------assign APA done\n')



ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
logFile.write(ts+'--------------------mergeing APA done\n')

logFile.write(ts+'------------------------------------------End\n')
# -------------------- Write output files
f = open(sys.argv[3]+".txt","w")
#f.write("geneName\tpeakName\tpeak\tchrom\tstrand\tcodingType\tApaRegion\tApaRegionStart\tApaRegionEnd\n")
for elem in D:
	elem=[str(tmp) for tmp in elem]
	f.write('\t'.join(elem))
	f.write('\n')
f.close()
logFile.close()




