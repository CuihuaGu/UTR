# 2015-11-04
# Cuihua Gu
# Genome visualization one by one
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
	# 2. scripts to run
# USAGE

args <- commandArgs(TRUE)
pvalFile = args[1] # /data/guc1/UTR/Code/tmp.pval.txt
drecFile = args[2] # /data/guc1/UTR/Code/tmp.drec.txt
pairFile = args[3] #"/data/guc1/UTR/Data/TCGA/PairInfor/KIRP_Tumor_Normal_dir.txt"
outpref = args[4] # ""
canDir = args[5] # "/data/guc1/UTR/Data/TCGA/RNAseq/KIRP/SolidTumor/"
norDir = args[6] # "/data/guc1/UTR/Data/TCGA/RNAseq/KIRP/Normal/"
extend = args[7] # 0
minDynamicNum = as.integer(args[8])
shortenFile = args[9] # LUSC_1e-06.2fold_short.txt
lengthenFile = args[10] # LUSC_1e-06.2fold_long.txt
pCutoff = as.numeric(args[11])

findArgument <- function(pair,pval,canDir,norDir,extend,outpref,drec,d,pCutoff){
	canFile = c()
	norFile = c()
	apa = c()
	out = c()
	exten = c()
	ind = c()
	for(i in 1:dim(pval)[1]){
		drec_li = drec[i,]
		pval_li = pval[i,]
		if(min(pval_li[,drec_li==d])<=pCutoff){
		ind_li = as.integer(which(pval_li == min(pval_li[,drec_li==d])))[1]
		ind_li = as.integer(strsplit(colnames(pval)[ind_li],"p")[[1]][2])
		canFile1 = paste(canDir,"/",pair[ind_li,2],"/",pair[ind_li,3],sep="")
		norFile1 = paste(norDir,"/",pair[ind_li,4],"/",pair[ind_li,5],sep="")
		canFile = c(canFile,canFile1)
		norFile = c(norFile,norFile1)
		apa = c(apa,rownames(pval)[i])
		out = c(out,paste(outpref,"_",rownames(pval)[i],"_",min(pval_li[,drec_li==d]),sep=""))
		exten = c(exten,extend)
		ind = c(ind,ind_li)
		}
	}
	res = cbind("Rscript","/data/guc1/UTR/Code/15.Sep.GenomeVisulization.R",canFile,norFile,apa,out,exten,ind)
	return(res)
}

callCmd <- function(LI){
	for(i in 1:length(LI)){
		system(LI[i])
	}
}

pval = read.table(pvalFile,header=T,row.names=1,sep="\t")
drec = read.table(drecFile,header=T,row.names=1,sep="\t")
Lon = as.character(read.table(lengthenFile)[,1])
Shor = as.character(read.table(shortenFile)[,1])
pair = read.table(pairFile,header=F)
s = length(Shor)
l = length(Lon)
pval_s = pval[Shor,]
drec_s = drec[Shor,]
pval_l = pval[Lon,]
drec_l = drec[Lon,]

short_arg = findArgument(pair,pval_s,canDir,norDir,extend,paste(outpref,"_short",sep=""),drec_s,"shorten",pCutoff)
long_arg = findArgument(pair,pval_l,canDir,norDir,extend,paste(outpref,"_long",sep=""),drec_l,"lengthen",pCutoff)
short_arg = apply(short_arg,1,paste,collapse=" ")
long_arg = apply(long_arg,1,paste,collapse=" ")
write.table(short_arg,paste(outpref,"_short_command.txt",sep=""),quote=F,row.names=F,col.names=F)
write.table(long_arg,paste(outpref,"_long_command.txt",sep=""),quote=F,row.names=F,col.names=F)
print(s)
print(l)
#callCmd(short_arg)
#callCmd(long_arg)
