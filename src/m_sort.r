#!/bin/env Rscript
# author: ph-u
# script: m_sort.r
# desc: remap cf425 data to microbiota, functional drugs category, mutation
# in: Rscript m_sort.r
# out: ?
# arg: 0
# date: 20230224

##### import #####
ptIN = "../../02_cf425_vitamin/"
load(paste0(ptIN,"raw/cf425_20230222.rda"))
source(paste0(ptIN,"src/plotSRC.r"))
ptOT = "../p_20230224/"
cBp = palette.colors(palette = "Okabe-Ito", alpha=1, recycle = T)[c(2,9,3,4,5,6,7,1,8)]

##### segregate #####
# DF508 homozygotes
# CFTRm[0], antimicrobials[1], pancrelipase[1], others[1] -> group vitamins
ddF508 = muT$regid_anon[which(muT[,3]==1 & rowSums(muT[,-c(1:2)])==1)] # 6419/12651 (all with Micro & Med data)
fIlter = which(mEd0$`CFTR modulators`==0 & mEd0$antimicrobials>0 & mEd0$pancrelipase>0 & mEd0$regid_anon %in% ddF508) # 43322
mIc1 = mIc0[fIlter,]
mEd1 = mEd0[fIlter,-which(colnames(mEd0) %in% c("CFTR modulators","antimicrobials","pancrelipase"))]

##### med groups: beta agonist, deoxyribonuclease, mucolytic agent, others #####
gP = c("beta agonist","deoxyribonuclease","mucolytic agent","others")
mEd2 = as.data.frame(matrix(0,nr=nrow(mEd1),nc=length(gP)+2))
colnames(mEd2) = c(colnames(mEd1)[1:2],gP); mEd2[,1:2] = mEd1[,1:2]
gP0 = c();for(i in 1:(length(gP)-1)){
	gP0 = c(gP0,grep(gP[i],colnames(mEd1)))
	mEd2[,2+i] = mEd1[,gP0[length(gP0)]]
};rm(i)
mEd2$others = rowSums(mEd1[,-c(1,2,gP0)])
mEd2[,-c(1:2)] = mEd2[,-c(1:2)]>0

yEs = which(rowSums(mEd2[,-c(1:2)])>0)
mIc1 = mIc1[yEs,]; mEd2 = mEd2[yEs,] # 41188 (exact group 0111)

##### prescriptions timeline #####
mEdts = as.data.frame(matrix(0,nr=length(unique(mEd2$year)),nc=ncol(mEd2)))
colnames(mEdts) = c(colnames(mEd2)[-1],"total");mEdts$year = unique(mEd2$year)
for(i in 1:nrow(mEdts)){
        x = mEd2[which(mEd2$year==mEdts$year[i]),-c(1:2)]
        mEdts[i,-1] = c(colSums(x),nrow(x))
};rm(i,x)
mEdts = mEdts[order(mEdts$year),]

## plot time-series
pdf(paste0(ptOT,"ts_med.pdf"), width=14, height=21)
medtsPlot(mEdts)
invisible(dev.off())

##### extract useable year fraction #####
yR = mEdts$year[apply(mEdts,1,function(s){all(s>0)})]
mEd3 = mEd2[which(mEd2$year %in% yR),]
mIc2 = mIc1[which(mIc1$year %in% yR),]

mEd3$group = paste0(as.numeric(mEd3[,3]),as.numeric(mEd3[,4]),as.numeric(mEd3[,5]),as.numeric(mEd3[,6]))
rEc = vector(mode="list",length=length(unique(mEd3$group)))
sAm = as.data.frame(matrix(0,nr=length(yR),nc=length(rEc)+1))
names(rEc) = colnames(sAm)[-1] = unique(mEd3$group)
colnames(sAm)[1] = "year"; sAm$year = yR
for(i in 1:length(rEc)){
        rEc[[i]] = as.data.frame(matrix(0,nr=length(yR),nc=ncol(mIc2)-1))
        colnames(rEc[[i]]) = colnames(mIc2)[-1]; rEc[[i]]$year = yR
        for(i0 in 1:length(yR)){
                x = mIc2[which(mEd3$year==yR[i0] & mEd3$group==names(rEc)[i]),]
                rEc[[i]][i0,-1] = c(colSums(x[,-c(1:2)])/nrow(x))
                sAm[which(sAm$year==yR[i0]),which(colnames(sAm)==names(rEc)[i])] = nrow(x)
}};rm(i,i0,x)

## eliminate entries/data under threshold (each year sample size >= 10)
i = which(sAm<10);i0 = i%%nrow(sAm); i0[i0==0] = nrow(sAm)
if(length(i)>0){ for(i1 in 1:length(i)){
        rEc[[ceiling(i[i1]/nrow(sAm))-1]][i0[i1],] = NA
};rm(i,i0,i1)}
i0 = c();for(i in 1:length(rEc)){
        rEc[[i]] = rEc[[i]][which(!is.na(rEc[[i]]$year)),]
        if(nrow(rEc[[i]])<3){i0 = c(i0,i)}
};for(i in rev(i0)){rEc[[i]] = NULL};rm(i,i0)

##### export #####
write.csv(sAm,paste0(ptOT,"m_sample.csv"), quote=F, row.names=F)
for(i in 1:length(rEc)){ for(i0 in 1:(nrow(rEc[[i]])-2)){
        write.csv(rEc[[i]][i0:(i0+2),],paste0(ptOT,"ts/",names(rEc)[i],"_",substr(rEc[[i]]$year[i0],3,4),substr(rEc[[i]]$year[i0+2],3,4),".csv"), quote=F, row.names=F)
}};rm(i,i0)
save(rEc,sAm,file=paste0(ptOT,"m_ts.rda"))
