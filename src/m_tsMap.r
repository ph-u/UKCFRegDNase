#!/bin/env Rscript
# author: ph-u
# script: m_tsMap.r (cp from v_tsMap.r)
# desc: time-series simulations on data
# in: Rscript m_tsMap.r
# out: ?
# arg: 0
# date: 20230224

##### import #####
source(paste0("../../03_cf425_mucus/src/plotSRC.r")); source("../../00_biLVC/src/src.r"); library(deSolve)
ptIN = "../p_20230224/data/"; ptOT = "../p_20230224/out/"
load(paste0(ptIN,"../m_ts.rda"))
nRep = 7; simO = 500; gpMc = c("beta agonist","deoxyribonuclease","mucolytic agent","others")
f = list.files(ptIN,"filter")
yR = c(); for(i in 1:length(rEc)){yR = c(yR,rEc[[i]][,1])};rm(i);yR = range(yR)

##### plot by group #####
for(i in 1:length(rEc)){
	cat("Plotting",names(rEc)[i],": data, ")
## data points
	nAm = ifelse(strsplit(names(rEc)[i],"")[[1]]==0,"",gpMc)
	pdf(paste0(ptOT,"m_ts_",names(rEc)[i],".pdf"), width=14, height=14)
	par(mar=c(5,5,1,0)+.1, mfrow = c(2,1), cex.axis=1.4, xpd=F)
	matplot(rEc[[i]][,1],rEc[[i]][,-1]*100, "p", pch=3+2:ncol(rEc[[i]]), xaxt="n",xlab="", ylab="Prevalence (%)", cex=5, col=cBp, cex.axis=2.1, cex.lab=2.1, ylim=c(0,100), xlim=yR)
	abline(h=0, col="#000000ff")
	axis(1,at=seq(yR[1],yR[2]), padj=.7, labels=paste0(seq(yR[1],yR[2]),"\n(",sAm[,which(colnames(sAm)==names(rEc)[i])],")\n[",rowSums(sAm[,-1]),"]"))
	mtext(paste0("Year (Group Sample Size) [Total Sample Size]\n",paste0(capFirst(nAm)[nAm!=""], collapse=" + ")), side=1, padj=2.8, cex=2.1)
	text(yR[1]-.05, 93, labels="n =", cex=2)

## simulations
	cat("simulations, ")
	f0 = f[grep(names(rEc)[i],f)]
	for(i0 in 1:length(f0)){
		stEd = strsplit(strsplit(f0[i0],"-")[[1]][1],"_")[[1]][2]
		stEd = 2000+c(as.numeric(substr(stEd,1,2)),as.numeric(substr(stEd,3,4)))
		x = rEc[[names(rEc)[which(names(rEc)==strsplit(f0[i0],"_")[[1]][1])]]]
		x0 = as.numeric(x[which(x[,1]==stEd[1]),-1])*100
		a = read.csv(paste0(ptIN,f0[i0]), header=T)
		if(nrow(a)>0){
			s = read.csv(paste0(ptIN,sub("filter","seed",f0[i0])), header=T)
			for(i1 in 1:nrow(a)){
				set.seed(s[a[i1,1],2])
				a0 = solveLV(x0, as.numeric(a[i1,-1]),stEd,"g")
				matplot(a0[,1], a0[,-1], type="l", add=T, col=cBl, lty=1+(1:(ncol(a0)-1))%%5)
		}}else{s=i1=a0=0}; text(.1+stEd[1],ifelse(stEd[1]%%2==0, 97,95), labels=nrow(a), cex=2)
	};rm(i0,stEd,x,x0,i1,a0,s,a)

## legend
	cat("legend\n"); plot.new(); x = colnames(rEc[[i]])[-1]
	lPt = legPlot(x, nDim=3)
	legend(x=c(-.12,1.04),y=c(.85,.45), legend=capFirst(x), title=paste("Microbial categories:",nRep,"replicates,",simO,"simulation each; N =",nRep*simO), border=NA, xpd=T, ncol=lPt[[2]], pch=3+2:ncol(rEc[[i]]), col=cBp, cex=2.1)
	invisible(dev.off())
};rm(i)
