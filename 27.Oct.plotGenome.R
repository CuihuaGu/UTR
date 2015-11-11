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
suppressMessages(library(Gviz))
suppressMessages(library(GenomicFeatures))

#-------------------------INPUT
args <- commandArgs(TRUE)
pvalFile = args[1] # /data/guc1/UTR/Code/tmp.pval.txt
drecFile = args[2] # /data/guc1/UTR/Code/tmp.drec.txt
pairFile = args[3] #"/data/guc1/UTR/Data/TCGA/PairInfor/KIRP_Tumor_Normal_dir.txt"
outpref = args[4] # ""
canDir = args[5] # "/data/guc1/UTR/Data/TCGA/RNAseq/KIRP/SolidTumor/"
norDir = args[6] # "/data/guc1/UTR/Data/TCGA/RNAseq/KIRP/Normal/"
extend = args[7] # 0
minDynamicNum = args[8]

#-------------------------FUNCTION
ApaDb <- function(apaList){
	d = read.table("/data/guc1/UTR/Output/2_ApaPair/PairwiseApa_GViz.bed",sep="\t",row.names=4)
	d2 = d[as.character(apaList),]
	return(d2)
}

ApaHighLiPro <- function(apaList){
	d = read.table("/data/guc1/UTR/Output/2_ApaPair/PairwiseApa.bed",sep="\t",row.names=4)
	d2 = d[paste(as.character(apaList),"pro",sep="_"),]
	return(d2)
}

ApaHighLiDis <- function(apaList){
	d = read.table("/data/guc1/UTR/Output/2_ApaPair/PairwiseApa.bed",sep="\t",row.names=4)
	d2 = d[paste(as.character(apaList),"dis",sep="_"),]
	return(d2)
}

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
	res = cbind(out,apa,exten,canFile,norFile)
	return(res)
}

plotGviz <- function(outpref,apaList,extend,cancerFiles,normalFiles,txdb,atlas){
	cat(format(Sys.time(), "%Y-%m-%d %I:%M:%S%p"))
	cat("\n")
	ApaList = ApaDb(apaList)
	pro = ApaHighLiPro(apaList)
	dis = ApaHighLiDis(apaList)
	extend = as.integer(extend)
	pdf(paste(outpref,".pdf",sep=""),width=30,height=20)
	chrom=as.character(ApaList[1,1])
	afrom=as.integer(ApaList[1,2])-extend
	ato=as.integer(ApaList[1,3])+extend
	hlStart=c(as.integer(pro[1,2]),as.integer(dis[1,2]))
	hlEnd=c(as.integer(pro[1,3]),as.integer(dis[1,3]))
	cat(paste("-----------",rownames(ApaList[1,]),sep=""))
	cat("\n")
	canTrack=AlignmentsTrack(cancerFiles,fill="hotpink2",name="cancer")
	norTrack=AlignmentsTrack(normalFiles,fill="deepskyblue2",name="normal")
	cat("    1. cancer and normal track loaded\n")
	txTr <- GeneRegionTrack(txdb,genome='hg19',chromosome=chrom)
	atlasTrack <- AlignmentsTrack(atlas,fill="lightseagreen",name="APA")
	cat("    2. txtr and atlas loaded\n")
	ht <- HighlightTrack(c(canTrack,norTrack,atlasTrack,txTr), start = hlStart,end=hlEnd,chromosome = chrom)
	cat("    3. start to plotting\n")
	print(ht)
	plotTracks(ht,chromosome=chrom,from=afrom,to=ato,type="coverage",main=rownames(ApaList[1,]))
	dev.off()
}

#-------------------------MAIN
txdb <- loadDb("~/lib/hg19db_knownGene.sqlite")
atlas = "/home/guc1/public_html/unified-atlas.bam"
pval = read.table(pvalFile,header=T,row.names=1,sep="\t")
drec = read.table(drecFile,header=T,row.names=1,sep="\t")
pair = read.table(pairFile,header=F)
s = sum(rowSums(drec=="shorten")==minDynamicNum)
l = sum(rowSums(drec=="lengthen")==minDynamicNum)
pval_s = pval[rowSums(drec=="shorten")==minDynamicNum,]
pval_l = pval[rowSums(drec=="lengthen")==minDynamicNum,]
print(paste("short",s,sep=":"))
print(paste("long",l,sep=":"))

short_arg = findArgument(pair,pval_s,canDir,norDir,extend,paste(outpref,"_short",sep=""))
long_arg = findArgument(pair,pval_l,canDir,norDir,extend,paste(outpref,"_long",sep=""))
for(i in 1:dim(short_arg)[1]){
	plotGviz(short_arg[i,1],short_arg[i,2],short_arg[i,3],short_arg[i,4],short_arg[i,5],txdb,atlas)
	print(short_arg[i,2])
}
for(i in 1:dim(long_arg)[1]){
	plotGviz(long_arg[i,1],long_arg[i,2],long_arg[i,3],long_arg[i,4],long_arg[i,5],txdb,atlas)
	print(long_arg[i,2])
}

