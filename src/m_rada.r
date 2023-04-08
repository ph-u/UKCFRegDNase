#!/bin/env Rscript
# author: ph-u
# script: m_rada.r (cp from F508ddRadar.r)
# desc: radar chart
# in: Rscript m_rada.r
# out: ?
# arg: 0
# date: 20230403

##### env #####
library(fmsb);source(paste0("../../03_cf425_mucus/src/plotSRC.r"))
nRep = 7;simO = 500
ptOT = "../p_20230224/out/"
sEco = read.csv(paste0(ptOT,"m_taxa_PC.csv"), header=T, colClasses=c(rep("character",3),rep("numeric",10)))[,-1]
colnames(sEco) = c("Med_Gp","Species",substr(colnames(sEco)[-c(1:2)],4,nchar(colnames(sEco)[-c(1:2)])-2))

##### med group names #####
medGp = unique(sEco[,1])
spNam = unique(sEco[,2])
totSimu = nRep*simO*length(spNam)
lgele = c("beta agonists","deoxyribonuclease","mucolytic agents","supplementaries")
lG0 = c();lG1 = strsplit(medGp,"");for(i in 1:length(medGp)){
        lG0[i] = paste(capFirst(lgele[which(lG1[[i]]==1)]),collapse=" + ")};rm(i,lG1)

##### Overall plot #####
sOver = as.data.frame(matrix(0,nr=length(medGp)+2, nc=ncol(sEco)-2))
colnames(sOver) = colnames(sEco)[-c(1,2)]
row.names(sOver)[1:2] = c("max","min");sOver[1,] = max(sEco[,-c(1:2)])

pdf(paste0(ptOT,"radar_overalL.pdf"), width=14, height=14)
plot.new()
lPt = legPlot(medGp, nDim=length(lG0)) # capFirst(lG0)
legend("top", inset=c(0,0), legend = lG0, title=paste("Medication Group - 100% =",totSimu,"simulations"), border=NA, xpd=T, cex=2, ncol=lPt[[2]], pch = rep(19,length(lG0)), col = cBp)
invisible(dev.off())

pdf(paste0(ptOT,"radar_map.pdf"), width=14, height=28)
par(mar=c(1,1,1,1)+.1, mfrow=c(5,2), cex=1.8, xpd=T)
for(i in 1:length(medGp)){
	s0 = sEco[which(sEco$Med_Gp==medGp[i]),-c(1:2)]
	sOver[i+2,] = colSums(s0)/nrow(s0)
	radarchart(sOver[c(1,2,i+2),], pcol=cBp[i], pfcol=cBl[i], plty=1, cglcol = "#000000ff", cglwd=1.4, cglty=1, axislabcol="grey", vlcex=.8, axistype=1, maxmin=T, caxislabels=paste0(seq(0,round(sOver[1,1]*100),round(sOver[1,1]*100)/4),"%"), title=lG0[i])
};rm(i)
invisible(dev.off())

##### med group focus plot #####
pdf(paste0(ptOT,"radar_taxaL.pdf"), width=14, height=14)
plot.new()
lPt = legPlot(spNam, nDim=length(spNam)) # capFirst(lG0)
legend("top", inset=c(0,0), legend = spNam, title=paste("Taxonomic Category - 100% =",totSimu,"simulations"), border=NA, xpd=T, cex=2, ncol=lPt[[2]], pch = rep(19,length(lG0)), col = cBp)
invisible(dev.off())

sOver = as.data.frame(matrix(0,nr=length(spNam)+2, nc=ncol(sEco)-2))
colnames(sOver) = colnames(sEco)[-c(1,2)]
row.names(sOver)[1:2] = c("max","min");sOver[1,] = max(sEco[,-c(1:2)])

pdf(paste0(ptOT,"radar_taxaM.pdf"), width=14, height=28)
par(mar=c(1,1,1,1)+.1, mfrow=c(5,2), cex=1.8, xpd=T)
for(i in 1:length(spNam)){
        s0 = sEco[which(sEco$Med_Gp=="0100" & sEco$Species==spNam[i]),-c(1:2)]
        sOver[i+2,] = colSums(s0)/nrow(s0)
        radarchart(sOver[c(1,2,i+2),], pcol=cBp[i], pfcol=cBl[i], plty=1, cglcol = "#000000ff", cglwd=1.4, cglty=1, axislabcol="grey", vlcex=.8, axistype=1, maxmin=T, caxislabels=paste0(seq(0,round(sOver[1,1]*100),round(sOver[1,1]*100)/4),"%"), title=spNam[i])
};rm(i)
invisible(dev.off())

##### taxa clade plot #####
pdf(paste0(ptOT,"radar_comp.pdf"), width=14, height=10)
par(mar=c(1,1,1,1)+.1, mfrow=c(1,2), cex=1.8, xpd=T)

s1 = list("Fungi"=c(1,2,8),"Gram_Negatives"=c(3,5,7))
for(i in 1:length(s1)){
	sOver = sEco[which(sEco$Med_Gp=="0100" & sEco$Species %in% spNam[s1[[i]]]),]
	row.names(sOver) = sOver$Species; sOver = sOver[,-c(1:2)]
	radarchart(rbind(rep(max(sOver),ncol(sOver)),rep(0,ncol(sOver)),sOver), pcol=cBp[s1[[i]]], pfcol=cBl[s1[[i]]], plty=1, cglcol = "#000000ff", cglwd=1.4, cglty=1, axislabcol="grey", vlcex=.8, axistype=1, maxmin=T, caxislabels=paste0(seq(0,round(max(sOver)*100),round(max(sOver)*100)/4),"%"), title=gsub("_"," ",names(s1)[i]))
	legend("topright", inset=c(0,0), legend=spNam[s1[[i]]], border=NA, xpd=T, cex=1, pch=rep(19,nrow(sOver)), col=cBp[s1[[i]]])
};rm(i)

invisible(dev.off())
