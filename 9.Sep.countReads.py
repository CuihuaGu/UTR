# 2015-9-5
# Cuihua Gu
# Assign APA peaks to refSeq transcripts.
# INPUT
	# 1. pair information file
	# 2. TCGA cancer data directory
	# 3. TCGA normal data directory
	# 4. APA information bed file
	# 5. Output directory and prefix
# OUTPUT

# USAGE
	# python script.py 

import sys
import copy
import time
import datetime
from subprocess import call

#---------------------USER INPUT
pairFile = sys.argv[1]
dataDirCan = sys.argv[2]
dataDirNor = sys.argv[3]
bedFile = sys.argv[4]
outPref = sys.argv[5]
outID = int(sys.argv[6])

def getTCGAData(pairFile,dataDirCan,dataDirNor):
	f = open(pairFile)
	files = []
	for line in f.readlines():
		line = line.rstrip("\n").split("\t")
		pID = line[0]
		canD = line[1]
		canF = line[2]
		norD = line[3]
		norF = line[4]
		canFile = dataDirCan+"/"+canD+"/"+canF
		norFile = dataDirNor+"/"+norD+"/"+norF
		files.append([pID,canFile,norFile])
	f.close()
	return files


#---------------------MAIN
# 1. convert apa to bed
# 2. read TCGA data dir
tcgaFile = getTCGAData(pairFile,dataDirCan,dataDirNor)

# run intersectBed
Nm = []
logFile = open(outPref+".log","w")
ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
logFile.write(ts+'--------------------Start\n')

for i in range(0,len(tcgaFile)):
	canNm = "T"+str(outID+i)
	norNm = "N"+str(i+outID)
	pID = tcgaFile[i][0]
	canFile = tcgaFile[i][1]
	norFile = tcgaFile[i][2]
	#cmd1 = ["intersectBed","-abam",canFile,"-b",bedFile,"-bed","-wo","-s",">",outPref+canNm+".bed"]
	#cmd2 = ["intersectBed","-abam",norFile,"-b",bedFile,"-bed","-wo","-s",">",outPref+norNm+".bed"]
	#call(cmd1)
	#call(cmd2)
	cmd1 = ["sh","9.Sep.IntersectBed.sh",canFile,bedFile,outPref+"."+canNm]
	ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
	logFile.write(canNm+'----------------------'+ts+'\n')
	logFile.write(" ".join(cmd1)+'\n')
	call(cmd1)
	cmd2 = ["sh","9.Sep.IntersectBed.sh",norFile,bedFile,outPref+"."+norNm]
	ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
	logFile.write(norNm+'----------------------'+ts+'\n')
	logFile.write(" ".join(cmd2)+'\n')
	call(cmd2)
	Nm.append([canNm+"/"+norNm,pID])

for elem in Nm:
	logFile.write(": ".join(elem)+"\n")

ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
logFile.write(ts+'--------------------End\n')

logFile.close()
































