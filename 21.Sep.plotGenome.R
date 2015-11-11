# 2015-9-21
# Cuihua Gu
# Genome visualization
# INPUT
	# 1. Pval file
	# 2. Drec file
	# 3. TCGA paired sample file
	# 4. output prefix
	# 5. cancer directory
	# 6. normal directory
	# 7. extend
	# 8. minDynamicNum
# OUTPUT
	# 2. pdf
# USAGE

args <- commandArgs(TRUE)
pvalFile = args[1] # /data/guc1/UTR/Code/tmp.pval.txt
drecFile = args[2] # /data/guc1/UTR/Code/tmp.drec.txt
pairFile = args[3] #"/data/guc1/UTR/Data/TCGA/PairInfor/KIRP_Tumor_Normal_dir.txt"
outpref = args[4] # ""
canDir = args[5] # "/data/guc1/UTR/Data/TCGA/RNAseq/KIRP/SolidTumor/"
norDir = args[6] # "/data/guc1/UTR/Data/TCGA/RNAseq/KIRP/Normal/"
extend = args[7] # 0
minDynamicNum = args[8]

findArgument <- function(pair,pval,canDir,norDir,extend,outpref){
	canFile = c()
	norFile = c()
	apa = c()
	out = c()
	exten = c()
	for(i in 1:dim(pval)[1]){
		pval_li = pval[i,]
		ind= as.integer(which(pval_li == min(pval_li)))
		canFile1 = paste(canDir,"/",pair[ind,2],"/",pair[ind,3],sep="")
		norFile1 = paste(norDir,"/",pair[ind,4],"/",pair[ind,5],sep="")
		canFile = c(canFile,canFile1)
		norFile = c(norFile,norFile1)
		apa = c(apa,rownames(pval)[i])
		out = c(out,paste(outpref,"_",rownames(pval)[i],".pdf",sep=""))
		exten = c(exten,extend)
	}
	res = cbind("Rscript","/data/guc1/UTR/Code/15.Sep.GenomeVisulization.R",canFile,norFile,apa,out,exten)
	return(res)
}

callCmd <- function(LI){
	for(i in 1:length(LI)){
		system(LI[i])
	}
}

pval = read.table(pvalFile,header=T,row.names=1,sep="\t")
drec = read.table(drecFile,header=T,row.names=1,sep="\t")
pair = read.table(pairFile,header=F)
s = sum(rowSums(drec=="shorten")>=minDynamicNum)
l = sum(rowSums(drec=="lengthen")>=minDynamicNum)
pval_s = pval[rowSums(drec=="shorten")>=minDynamicNum,]
pval_l = pval[rowSums(drec=="lengthen")>=minDynamicNum,]

short_arg = findArgument(pair,pval_s,canDir,norDir,extend,paste(outpref,"_short",sep=""))
long_arg = findArgument(pair,pval_l,canDir,norDir,extend,paste(outpref,"_long",sep=""))
short_arg = apply(short_arg,1,paste,collapse=" ")
long_arg = apply(long_arg,1,paste,collapse=" ")

callCmd(short_arg)
callCmd(long_arg)

