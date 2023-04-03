#!/bin/env Rscript
# author: ph-u
# script: dataSimplify.r
# desc: remap cf425 data to microbiota, functional drugs category, mutation
# in: Rscript dataSimplify.r
# out: ?
# arg: 0
# date: 20230222

##### import #####
cat("data reading, environment preparation:",date(),"\n")
ptIN="../../01_cf425/data/"; ptOT="../raw/"
load(paste0(ptIN,"cf425MedMic.rda"));
muT = read.csv(paste0(ptIN,"mutPWCF.csv"), header=T, stringsAsFactors=F);

##### remap medication history 20230210 #####
cat("remap medication:",date(),"\n")
dRug = read.csv(paste0(ptIN,"F508dd_drugs.csv"), header=T)
d0 = unique(unlist(strsplit(unique(dRug$src), split = ";")), na.rm=T)
d0 = d0[order(d0)];d0 = d0[!is.na(d0)] # 210 medication groups
mEd0 = as.data.frame(matrix(0,nr=nrow(mEdic),nc=length(d0)+2))
colnames(mEd0) = c(colnames(mEdic)[1:2],d0)
for(i in which(!is.na(dRug[,2]))){
        d1 = unlist(strsplit(dRug[i,2],";"))
        mEd0[,d1] = mEd0[,d1] + mEdic[,dRug[i,1]]};rm(i,d1)
mEd0 = as.data.frame(mEd0>0)
mEd0[,1:2] = mEdic[,1:2]

## plot medication time-series
d1 = as.data.frame(matrix(0,nr=length(unique(mEd0$year)),nc=ncol(mEd0)))
colnames(d1) = c(colnames(mEd0)[-1],"total"); d1[,1] = unique(mEd0$year)
for(i in 1:nrow(d1)){
        x = mEd0[which(mEd0$year==d1$year[i]),-c(1:2)]
        d1[i,-1] = c(colSums(x),nrow(x))
};rm(i,x)
d2 = d1[,-c(1,ncol(d1))]
for(i in 1:nrow(d2)){d2[i,] = d2[i,]/d1$total[i]};rm(i)
d3 = t(d2);colnames(d3) = d1$year
d4 = which(apply(d3,1,function(s){any(s>.3)})>0)
pdf(paste0(ptOT,"ts_med.pdf"), width=21);par(mar=c(5,5,1,0)+.1)
matplot(1:nrow(d3),d3*100, type="l", ylim = c(-40,100), ylab = "Presence (%)", xlab = "Medication (alphabetical order)", cex.axis = 2, cex.lab=2)
text(d4,-23,names(d4), srt=90)
abline(h=30);text(28,55,"Threshold = 30%", cex=2)
legend("topleft",legend=colnames(d3),lty=1:ncol(d3),lwd=7, hori=T, cex=1.5, col=1:ncol(d3))
invisible(dev.off())


##### remap microbiota #####
cat("remap microbiota:",date(),"\n")
m0 = strsplit(colnames(mIcro)[-c(1:2)], " ")
m1 = c();for(i in 1:length(m0)){m1 = c(m1,m0[[i]][1])};m0 = unique(m1);rm(i)

## genus (184) concat
m2 = as.data.frame(matrix(0,nr=nrow(mIcro),nc=length(m0)+2))
colnames(m2) = c(colnames(mIcro)[1:2],m0); m2[,1:2] = mIcro[,1:2]
for(i in 1:length(m0)){ if(length(grep(m0[i],m1))>1){
		m2[,i+2] = (rowSums(mIcro[,2+which(m1==m0[i])])>0)
	}else{
		m2[,i+2] = mIcro[,2+which(m1==m0[i])]
}};rm(i)

## time-series genera abundances
m3 = as.data.frame(matrix(0,nr=length(unique(m2$year)),nc=ncol(m2)))
colnames(m3) = c(colnames(m2)[-1],"total"); m3[,1] = unique(m2$year)
for(i in 1:nrow(m3)){
	x = m2[which(m2$year==m3$year[i]),-c(1:2)]
	m3[i,-1] = c(colSums(x),nrow(x))
};rm(i,x)

## plot genus time-series
m4 = m3[,-c(1,ncol(m3))]
for(i in 1:nrow(m4)){m4[i,] = m4[i,]/m3$total[i]};rm(i)
m4 = t(m4);colnames(m4) = m3$year
m5 = which(apply(m4,1,function(s){any(s>.05)})>0)
pdf(paste0(ptOT,"ts_genus.pdf"), width=21);par(mar=c(5,5,1,0)+.1)
matplot(1:nrow(m4),m4*100, type="l", ylim = c(-30,100), ylab = "Presence (%)", xlab = "Genus (alphabetical order)", cex.axis = 2, cex.lab=2)
text(m5,-18,names(m5), srt=90)
abline(h=5);text(60,10,"Threshold = 5%", cex=2)
legend("topleft",legend=colnames(m4),lty=1:ncol(m4),lwd=7, hori=T, cex=1.5, col=1:ncol(m4))
invisible(dev.off())

## sort genus
m6 = as.data.frame(matrix(0,nr=nrow(m2),nc=length(m5)+3))
colnames(m6) = c(colnames(m2)[1:2],names(m5),"others")
m6[,1:2] = m2[,1:2]
for(i in 3:ncol(m2)){if(colnames(m2)[i] %in% names(m5)){
	m6[,colnames(m2)[i]] = m2[,i]
}else{
	m6$others = m6$others + m2[,i]
}};rm(i)
m6$Mycobacteria = ((m6$Mycobacterium + m6$Mycobacteroides)>0)
m6$Mycobacterium = m6$Mycobacteroides = NULL
m6$others = (m6$others>0)

##### export #####
gEnus = m2; mIc0 = m6[,c(1:5,12,6:8,10,9,11)]
write.csv(d1,paste0(ptOT,"ts_med.csv"), quote=F, row.names=F)
write.csv(m3,paste0(ptOT,"ts_genus.csv"), quote=F, row.names=F)
save(mIc0,mEd0,gEnus,muT,file=paste0(ptOT,"cf425_20230222.rda"))
