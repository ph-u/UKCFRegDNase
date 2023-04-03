#!/bin/env Rscript
# author: ph-u
# script: m_clus.r (cp from v_clus.r)
# desc: multivariate analysis (PCA, PLS-DA, hclust)
# in: Rscript m_clus.r
# out: ?
# arg: 0
# date: 20230224

##### import #####
source(paste0("../../03_cf425_mucus/src/plotSRC.r"))
library(ggbiplot); library(mixOmics)
ptIN = "../p_20230224/data/"; ptOT = "../p_20230224/out/"
nRep = 7; simO = 500; load(paste0(ptIN,"../m_ts.rda"))
f = list.files(ptIN,"eco");f = f[-grep("0110",f)]
nAm = matrix(unlist(strsplit(matrix(unlist(strsplit(f,"-")),nr=2)[1,],"_")),nr=2)

## coagulate pairwise ecology
for(i in 1:length(f)){
	x = read.csv(paste0(ptIN,f[i]), header=T)
	if(nrow(x)>0){
		x$group = nAm[1,i];x$end = as.numeric(substr(nAm[2,i],3,4))
		if(i==1){f0 = x}else{f0 = rbind(f0,x)}
}};rm(i,x);f0$end = 2000+f0$end
f0$pair = paste0(f0$category1," - ",f0$category2)
eCo = unique(f0$c1_is); eCo = data.frame(src=eCo,role=eCo[c(1,4,7,2,5,8,3,6,9)],type=paste0("<-",c("Mutualism","Commensalism","Predatory/Parasitism","Commensalism","Neutral/no interaction","Amensalism","Predatory/Parasitism","Amensalism","Competition"),"->"))

## prescriptions
gpMc = c("beta agonist","deoxyribonuclease","mucolytic agent","others")
gP = data.frame(src=unique(f0$group),nAm=NA)
for(i in 1:nrow(gP)){
	x = ifelse(strsplit(gP$src[i],"")[[1]]==0,"",gpMc)
	gP[i,2] = paste0(gsub("[.]"," ",capFirst(x)[x!=""]), collapse=" + ")
};rm(i,x)

##### pairwise ecology #####
pE = f0[,c("group","end","pair","c1_is","count")]
uqPE = unique(pE[,c("group","pair")]); row.names(uqPE) = 1:nrow(uqPE)
uqPE = cbind(uqPE,matrix(0,nr=nrow(uqPE),nc=length(unique(eCo$type))))
colnames(uqPE)[-c(1:2)] = unique(eCo$type)
for(i in 1:nrow(uqPE)){
	x = f0[which(f0$group==uqPE$group[i] & f0$pair==uqPE$pair[i]),]
	for(i0 in 1:nrow(eCo)){
		uqPE[i,which(colnames(uqPE)==eCo$type[i0])] = sum(x$count[which(x$c1_is==eCo$src[i0])])/(nRep*simO*length(unique(x$end)))
}};rm(i,x,i0)
uqPE[,"<-Unpredictable->"] = 1-rowSums(uqPE[,-c(1:2)])
for(i in 1:nrow(gP)){uqPE$group[which(uqPE$group==gP$src[i])] = gP$nAm[i]};rm(i)

## medication PCA
p = prcomp(uqPE[,-c(1:2)], scale.=T)
pCa0 = ggbiplot(p, var.scale=1, groups=uqPE$group, ellipse=T, ellipse.prob=.95, labels.size=4, varname.size = 4, varname.adjust = c(3,1.5,1.9,1,1.4,3,1.8))+ #labels=row.names(uqPE)
	scale_color_manual(values=setNames(cBp[1:length(unique(uqPE$group))], unique(uqPE$group)))+
	guides(color=guide_legend(title="Prescription (Antimicrobials + Pancrelipase + ...)")) + theme_bw()+
	scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10))+
	coord_cartesian(xlim=c(-4.5,4.5))+
	theme(legend.position = 'bottom', legend.direction = "vertical",
        	panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
	        legend.title = element_text(size=18),
	        axis.text=element_text(size=16),axis.title=element_text(size=14),
	        plot.margin=margin(t=0,r=0,b=0,l=1))
ggsave(paste0(ptOT,"m_pair_PCA.pdf"), plot=pCa0, width=6, height=7)

## medication PLS-DA
p0 = plsda(uqPE[,-c(1:2)], as.factor(uqPE$group))
pdf(paste0(ptOT,"m_pair_PLSDA.pdf"), width=14, height=14)
plotIndiv(p0, ellipse=T, ellipse.level=.95, col=cBp[c(9,10,8,7,6,3,5,4,2,1)], size.xlabel=rel(1.4), size.ylabel=rel(4), size.axis=rel(2), legend=F, style="graphics", cex=3, pch=20)
invisible(dev.off())

##### medication cluster #####
sP = unique(c(f0$category1,f0$category2))
eRol = data.frame(group=rep(gP$src,each=length(sP)),taxa=sP)
eRol = cbind(eRol,matrix(0,nr=nrow(eRol),nc=nrow(eCo)+1))
colnames(eRol)[-c(1:2)] = c(eCo$src,"unpredictable")

for(i in 1:nrow(eRol)){
	x = f0[intersect(grep(eRol$taxa[i],f0$pair),which(f0$group==eRol$group[i])),]
	x0 = which(x$category1!=x$category2 & x$category2==eRol$taxa[i])
	if(length(x0)>0){ for(i0 in x0){x[i0,1:2] = x[i0,2:1];x$c1_is[i0] = eCo$role[which(eCo$src==x$c1_is[i0])]}}
	for(i0 in 3:(ncol(eRol)-1)){
		eRol[i,i0] = sum(x$count[which(x$c1_is==colnames(eRol)[i0])])
}};rm(i,i0,x,x0)
eRol[,-c(1:2,ncol(eRol))] = eRol[,-c(1:2,ncol(eRol))]/(nRep*simO*length(sP)*(nrow(sAm)-2))
eRol$unpredictable = 1-rowSums(eRol[,-c(1:2,ncol(eRol))])

colnames(eRol)[-c(1:2)] = paste0("<-",capFirst(c(gsub("_"," ",sub("_of","",sub("_c2","",eCo$src))),"unpredictable")),"->")

## taxa PCA
pTx = prcomp(eRol[,-c(1:2)], scale.=T)
pCaTx = ggbiplot(pTx, var.scale=1, groups=eRol$group, ellipse=T, ellipse.prob=.95, labels.size=4, varname.size = 3, varname.adjust = 2+c(0,-.7,0,1,-1,-.5,0,1,.5,0))+ # labels=row.names(eRol)
	scale_color_manual(values=setNames(cBp[1:length(unique(eRol$group))], unique(eRol$group)))+
	guides(color=guide_legend(title="Prescription (Antimicrobials + Pancrelipase + ...)")) + theme_bw()+
	scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10))+
	coord_cartesian(xlim=c(-4.5,4.5))+
	theme(legend.position = 'bottom', legend.direction = "vertical",
        	panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
	        legend.title = element_text(size=18),
	        axis.text=element_text(size=16),axis.title=element_text(size=14),
	        plot.margin=margin(t=0,r=0,b=0,l=1))
ggsave(paste0(ptOT,"m_taxa_PCA.pdf"), plot=pCaTx, width=6, height=7)

## taxa PLS-DA
p0Tx = plsda(eRol[,-c(1:2)], as.factor(eRol$group))
pdf(paste0(ptOT,"m_taxa_PLSDA.pdf"), width=14, height=14)
plotIndiv(p0Tx, ellipse=T, ellipse.level=.95, col=cBp[1:length(unique(eRol$group))], size.xlabel=rel(1.4), size.ylabel=rel(4), size.axis=rel(2), legend=F, style="graphics", cex=3, pch=20)
invisible(dev.off())

## hclust
row.names(eRol) = paste0("(",row.names(eRol),") ",eRol$taxa)
hc = hclust(dist(eRol[,-c(1:2)], method = "euclidean"), method = "complete");
pdf(paste0(ptOT,"m_taxa_clus.pdf"), width=14, height=14);
plot(hc, lwd=3);
invisible(dev.off());

## taxa plot
x = capFirst(gsub("_"," ",sub("_of","",sub("_c2","",eCo$src))))
lPt = legPlot(x, nDim=nrow(eCo))
for(i in 1:length(sP)){
	tX = eRol[grep(sP[i],eRol$taxa),]
	tX0 = as.matrix(t(tX[,-c(1:2)]))
	colnames(tX0) = gP$src#gsub(" + ","\n",gP$nAm)
	pdf(paste0(ptOT,"v_taxa_",sP[i],".pdf"), width=14, height=14)
	par(mar=c(5,4.5,2,1)+.1, mfrow=c(2,1), cex.axis=1.4, xpd=T)
	barplot(tX0[-nrow(tX0),]*100, ylab="Ecological role (%)", col=cBp, ylim=c(0,100), cex.lab=2, cex.axis=2)
	plot.new()
	legend("top", inset=c(0,0), legend = x, title=paste("Ecological Relationship - 100% =",nRep*simO*(nrow(sAm)-2),"simulations"), border=NA, xpd=T, cex=2, ncol=lPt[[2]], pch = rep(19,nrow(eCo)), col = cBp)
	invisible(dev.off())
};rm(i,x,lPt,tX,tX0)

##### export #####
write.csv(f0,paste0(ptOT,"m_pair_eco.csv"),quote=F,row.names=F)
write.csv(uqPE,paste0(ptOT,"m_pair_PC.csv"),quote=F,row.names=T)
write.csv(eRol,paste0(ptOT,"m_taxa_PC.csv"),quote=F,row.names=T)
save(f0,eCo,uqPE,file=paste0(ptOT,"../m_pairEco.rda"))
