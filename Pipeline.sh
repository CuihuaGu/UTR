# 2015-9-24
# Cuihua Gu
# Run the pipeline
# STEPS
	# 1. 9.Sep.countReads.py 1 2 apaPair 3 4
	# 2. 15.Sep.MergeItstBed.sh 4 5 6
	# 3. 21.Sep.CalcDiff.R table 7 outpref
	# 4. 21.Sep.plotGenome.R pvalFile drecFile 1 outpref 2 apaPair 9
# INPUT
	# 1. /data/guc1/UTR/Data/TCGA/PairInfor/KIRP_Tumor_Normal_dir.txt
	# 2. /data/guc1/UTR/Data/TCGA/RNAseq/KIRP/SolidTumor
	# 3. /data/guc1/UTR/Data/TCGA/RNAseq/KIRP/Normal
	# 4. outDir /data/guc1/UTR/Output/3_countReads/KIRP/KIRP
	# 5. outpref
	# 6. IDStart 
	# 7. IDEnd
	# 8. chi-square p-value cutoff
	# 9. Number of patient cutoff
	# 10. Extend
	# 11. minimun dynamic number of sample
# OUTPUT
	# 2. 
# USAGE
	# sh run_Pipes.sh /data/guc1/UTR/Data/TCGA/PairInfor/KIRP_Tumor_Normal_dir.txt /data/guc1/UTR/Data/TCGA/RNAseq/KIRP/SolidTumor /data/guc1/UTR/Data/TCGA/RNAseq/KIRP/Normal /data/guc1/UTR/Output/KIRP 1 2 0.01 2 1000 2
# -----------------INPUT VARIABLES
Pair=$1
canDir=$2
norDir=$3
#apaPair=/data/guc1/UTR/Output/2_ApaPair/26.Sep.Check.bed
#apaPair=/data/guc1/UTR/Output/2_ApaPair/PairwiseApa.bed
apaPair=/data/guc1/UTR/Output/2_ApaPair/PairwiseApa_rmOtherThanLastExon_less20kb.bed
outDir=$4 # /data/guc1/UTR/Output/KIRP
outPref=$5 
idStart=$6
idEnd=$7
pCutoff=$8
numCutoff=$9
extend=$10
minDynamicNum=$11
fcCutoff=$12
filterFile=$13
nm=`basename $filterFile|sed 's/\.txt//g'`
# -----------------DEFINE VARIABLES
mkdir -p ${outDir}/ReadCounts
mkdir -p ${outDir}/CalcDiff
mkdir -p ${outDir}/Visualize
logFile=${outDir}/${outPref}.`date|awk -F ' ' '{printf($2"."$3"."$6"."$4)}'`.log
readDir=${outDir}/ReadCounts/${outPref}
diffDir=${outDir}/CalcDiff/${outPref}.${nm}
picDir=${outDir}/Visualize/${outPref}.${nm}
sed -n "${idStart},${idEnd}p" $Pair > ${outDir}/${outPref}.tmp
pairFile=${outDir}/${outPref}.tmp
echo $diffDir
# -----------------Run pipes
echo ----------------------------Start >> $logFile
echo `date` >> $logFile

# 1. 9.Sep.countReads.py 1 2 apaPair 3 4
# 2. 15.Sep.MergeItstBed.sh 4 5 6
echo -------------Counting reads >> $logFile
echo `date` >> $logFile
#python /data/guc1/UTR/Code/9.Sep.countReads.py ${pairFile} ${canDir} ${norDir} ${apaPair} ${readDir} ${idStart}
#sh /data/guc1/UTR/Code/15.Sep.MergeItstBed.sh ${outDir}/ReadCounts/${outPref} ${idStart} ${idEnd}

# 3. 21.Sep.CalcDiff.R table 8 9 outpref
echo -------------Calculating dynamic APA >> $logFile
echo `date` >> $logFile
Rscript /data/guc1/UTR/Code/21.Sep.CalcDiff.R ${outDir}/ReadCounts/${outPref}.TdTpNdNp.txt ${pCutoff} ${numCutoff} ${diffDir} ${fcCutoff} ${nm} ${filterFile}

# 4. 21.Sep.plotGenome.R pvalFile drecFile 1 outpref 2 apaPair 9
echo -------------Genome visualization >> $logFile
echo `date` >> $logFile
#Rscript /data/guc1/UTR/Code/27.Oct.plotGenome.R ${diffDir}.${numCutoff}.${pCutoff}.pval.txt ${diffDir}.${numCutoff}.${pCutoff}.drec.txt ${pairFile} ${picDir} ${canDir} ${norDir} ${extend} ${minDynamicNum}


echo ----------------------------End >> $logFile
echo `date` >> $logFile



