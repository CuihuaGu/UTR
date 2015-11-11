# 2015-9-5
# Cuihua Gu
# Script to visualize genome data
# INPUT
	# 1. Cancer bam files directory and name list
	# 2. Normal bam files directory and name list
	# 3. Name of apa which needs visualization
	# 4. output prefix
# OUTPUT
	# 2. Pdf file
# USAGE
library(Gviz)
library(GenomicFeatures)

args <- commandArgs(TRUE)
apaList = read.table(args[1])
outpref = args[2]
extend = as.integer(args[3])

txdb <- loadDb("~/lib/hg19db_knownGene.sqlite")
#genMod = read.table("/data/guc1/UTR/Data/Annotation/RefSeqGenesGvizModel.txt",header=T)
atlas = "/home/guc1/public_html/unified-atlas.bam"
#apa = read.table("/data/guc1/UTR/Output/2_ApaPair/PairwiseApa.bed",row.names=4,sep="\t")

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

ApaList = ApaDb(apaList[,1])
pro = ApaHighLiPro(apaList[,1])
dis = ApaHighLiDis(apaList[,1])

pdf(paste(outpref,".pdf",sep=""),width=20,height=8)
for(i in 1:dim(ApaList)[1]){
	print(ApaList[i,])
	#canTrack = c()
	#norTrack = c()
	chrom=as.character(ApaList[i,1])
	afrom=as.integer(ApaList[i,2])-extend
	ato=as.integer(ApaList[i,3])+extend
	hlStart=c(as.integer(pro[i,2]),as.integer(dis[i,2]))
	hlEnd=c(as.integer(pro[i,3]),as.integer(dis[i,3]))
	#for(j in 1:length(cancerFiles)){
	#	canTrack=c(canTrack,AlignmentsTrack(cancerFiles[j],genome="hg19",chromosome=chrom,type='coverage',col="hotpink2",fill="hotpink2",name=paste("Cancer",j,sep="-"),background.title ="white", fontcolor.title="grey31", col.axis="grey31"))
	#	norTrack=c(norTrack,AlignmentsTrack(normalFiles[j],genome="hg19",chromosome=chrom,type='coverage',col="deepskyblue2",fill="deepskyblue2",name=paste("Normal",j,sep="-"),background.title ="white", fontcolor.title="grey31", col.axis="grey31"))
	#}
	atlasTrack <- AlignmentsTrack(atlas,genome="hg19",chromosome=chrom,type='coverage',col="lightseagreen",fill="lightseagreen",name="APA",background.title ="white", fontcolor.title="grey31", col.axis="grey31")
	#print("atlasTrack")
	#grtrack <- GeneRegionTrack(genMod,genome='hg19',chromosome=chrom,name="refSeq",col="black",fill="black",background.title ="white", fontcolor.title="grey31", col.axis="grey31")
	txTr <-GeneRegionTrack(txdb,genome="hg19",chromosome=chrom,showId=TRUE,geneSymbol=TRUE,name="UCSC",col="black",fill="black",background.title ="white", fontcolor.title="grey31", col.axis="grey31")
	#print("txTr")
	hlTrack <- HighlightTrack(list(atlasTrack,txTr), genome="hg19",start = hlStart,end=hlEnd,chromosome = chrom)
	#print("ht")
	#print(ht)
	plotTracks(hlTrack,chromosome=chrom,genome="hg19",from=afrom,to=ato,type="coverage",main=rownames(ApaList[i,]))

}
dev.off()

