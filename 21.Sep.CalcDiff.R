# 2015-9-21
# Cuihua Gu
# Calculate chi-square test using R
# INPUT
	# 1. table file (KIRP.TdTpNdNp.txt)
	# 2. chi-square p-value cutoff
	# 3. Number of patient cutoff
	# 4. output prefix
	# 5. fcCutoff
# OUTPUT
	# 2. a table contains chi-square p-value
# USAGE

args <- commandArgs(TRUE)
InputFile = args[1]
CutOff = args[2] # pvalue cutoff
NumCutoff = args[3] # number cutoff
outpref = args[4]
fcCutoff=args[5]
filtername=args[6]
print(length(args))
if(length(args)==7){
	FilterFile=args[7]
	Filter_="True"
}
if(length(args)==5){
	Filter_="False"
}
print(args)
print(Filter_)
filterApa <- function(data,Filter_,FilterFile){
	if(Filter_=="False"){
		return(data)
	}
	else{
		f=as.character(read.table(FilterFile,sep="\t")[,2])
		return(data[f,])
	}
}

chisqPval <- function(li){
	li = as.numeric(strsplit(as.character(li),split=",")[[1]])
	if(sum(li) == 0){
		pval = 1
	}
	else{
		dim(li) = c(2,2)
		li = t(li)
		pval = chisq.test(li)$p.value
	}
	return(pval)
}

adjP <- function(pval){
	pval2 = p.adjust(pval,"BH")
	dim(pval2) = dim(pval)
	colnames(pval2) = colnames(pval)
	rownames(pval2) = rownames(pval)
	return(pval2)
}

chisqPval_multiple <- function(LI){
	t = c()
	n = c()
	for(i in 1:length(LI)){
		temp = as.numeric(strsplit(as.character(LI[i]),split=",")[[1]])
		t = c(t,c(temp[c(1,2)]))
		n = c(n,c(temp[c(3,4)]))
	}
	pval = chisq.test(rbind(t,n))$p.value
	return(pval)
}

changDirec <- function(li){
	li = as.numeric(strsplit(as.character(li),split=",")[[1]])
	canRatio = (li[1]+1)/(li[2]+1)
	norRatio = (li[3]+1)/(li[4]+1)
	if(sum(li) == 0){
		drec = "even"
	}
	else if(canRatio > norRatio){
		drec = "lengthen"
	}
	else if(canRatio < norRatio){
		drec = "shorten"
	}
	else{
		drec = "unknown"
	}
	return(drec)
}

norByLen <- function(li,ref){
	li = as.numeric(strsplit(as.character(li),split=",")[[1]])
	pro = as.integer(ref[1,1])
	dis = as.integer(ref[1,2])
	tmp = rep(0,nchar(round(mean(c(pro,dis))))-1)
	tmp2 = as.numeric(paste("1",gsub(", ","",toString(tmp)),sep=""))
	li[1] = li[1]/dis * tmp2 
	li[2] = li[2]/pro * tmp2 
	li[3] = li[3]/dis * tmp2 
	li[4] = li[4]/pro * tmp2 
	return(gsub(" ","",toString(li)))
}

getSigf <- function(pval,CutOff){
	CutOff = as.numeric(CutOff)
	SigNum = rowSums(pval<CutOff)
	pvalMean = rowMeans(pval)
	return(cbind(pval,pvalMean,SigNum))
}

getCanRatio <- function(li){
	li = as.numeric(strsplit(as.character(li),split=",")[[1]])
	canRatio = (li[1]+1)/(li[2]+1)
	return(canRatio)
}

getNorRatio <- function(li){
	li = as.numeric(strsplit(as.character(li),split=",")[[1]])
	norRatio = (li[3]+1)/(li[4]+1)
	return(norRatio)
}

plotRatio <- function(d,outpref,pval,NumCutoff,CutOff,fcCutoff){
	short=0
	lengthen=0
	fcCutoff=as.numeric(fcCutoff)
	pdf(paste(outpref,NumCutoff,CutOff,fcCutoff,"pdf",sep="."))
	canR = rowMeans(apply(d,c(1,2),getCanRatio))
	norR = rowMeans(apply(d,c(1,2),getNorRatio))
	plot(log2(canR),log2(norR),pch=20,cex=0.6,col="darkgrey",xlab="Cancer:log2(Distal/Proximal)",ylab="Normal:log2(Distal/Proximal)")
	ind = which(pval[,"SigNum"]>=as.numeric(NumCutoff))
	points(log2(canR[ind]),log2(norR[ind]),pch=20,cex=0.6,col="black")
	lines(c(-10,10),c(-10,10),col="red",lty=2)
	for(i in ind){
		if(canR[i]/norR[i]<(1/fcCutoff))
		{
		short=short+1
		points(log2(canR[i]),log2(norR[i]),pch=20,col="red",cex=0.8)
		}
		else if(canR[i]/norR[i]>fcCutoff)
		{
		lengthen=lengthen+1
		points(log2(canR[i]),log2(norR[i]),pch=20,col="blue",cex=0.8)
		}
	}
	legend("topleft",c(paste("Shorter 3'UTR in cancer: ",short),paste("Longer 3'UTR in cancer: ",lengthen)),pch=20,col=c("red","blue"),cex=1,bty="n")
	legend("bottomright",c(paste("p-value < ",CutOff),paste("Fold change > ",fcCutoff,sep="")),cex=1,bty="n")
	dev.off()
}

# MAIN FUCNTION
ref = read.table("/data/guc1/UTR/Output/2_ApaPair/PairwiseApa_Length.bed",sep="\t",row.names=1,header=T)
d = read.table(InputFile, header=T, row.names=1, sep='\t')
d=filterApa(d,Filter_,FilterFile)
d2 = as.matrix(d)
print("1")
# 1. normalize by length
for(i in 1:dim(d)[1]){
	ref_sub = ref[rownames(d[i,]),]
	for(j in 1:dim(d)[2]){
		d2[i,j] = norByLen(d[i,j],ref_sub)
	}
}
d = d2
print("2")
# 2.1. Calculate p-value using chi-square single
pval = apply(d, c(1,2), chisqPval)
pval[is.nan(pval)] <- 1
pval = adjP(pval)
pval = getSigf(pval,CutOff)
#pval2 = pval[which(pval[,"SigNum"]>=as.numeric(NumCutoff)),]
write.table(pval,paste(outpref,filtername,"pval.txt",sep="."),row.names=T,col.names=T,quote=F,sep="\t")

# 2.2 Calculate p-value using chi-square multiple


# 3.1. Decide change direction using single
drec = apply(d, c(1,2), changDirec)
#drec = drec[which(pval[,"SigNum"]>=as.numeric(NumCutoff)),]
write.table(drec,paste(outpref,filtername,"drec.txt",sep="."),row.names=T,col.names=T,quote=F,sep="\t")

# 3.2. Decide change direction using multiple

# 4. Plot 
plotRatio(d,outpref,pval,NumCutoff,CutOff,fcCutoff)

