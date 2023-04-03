#!/bin/env Rscript
# author: ph-u
# script: plotSRC.r
# desc: plotting functions
# in: source("plotSRC.r")
# out: //
# arg: 0
# date: 20230223

cBp = c(palette.colors(palette = "Okabe-Ito", alpha=1, recycle = T)[c(2,9,3,4,5,6,7,1,8)],palette.colors(palette = "alphabet", alpha=1, recycle = T))
cBl = c(palette.colors(palette = "Okabe-Ito", alpha=.01, recycle = T)[c(2,9,3,4,5,6,7,1,8)],palette.colors(palette = "alphabet", alpha=.01, recycle = T))

##### f: capitalise first letter #####
capFirst = function(x){return(paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x))))}

##### f: plot legend #####
legPlot = function(x,nDim=3){
        nDim = c(nDim, ceiling(length(x)/nDim))
        legMx = matrix(1:prod(nDim), nrow=nDim[1], ncol=nDim[2], byrow=F)
        legBd = rep("#000000ff", prod(nDim))
        legBd[legMx>length(x)] = legMx[legMx>length(x)] = NA
        return(list(legMx,nDim[2]))
}

##### f: medication plot #####
medtsPlot = function(x){
	par(mar=c(10,9,3,3)+.1, mfrow=c(2,1), xpd=F)
	matplot(x[,1],x[,-1],type="b", pch=20+(1:(ncol(x)-1)), lty=1, col=cBp, lwd=12, cex=7, xaxt="n", yaxt="n", xlab="", ylab="", cex.axis=3.5)
	axis(1, at=x[,1], labels=x[,1], lwd=10, cex.axis=4.2, padj=1.2, tck=-.02)
	xLab=x[seq(1,nrow(x), by=3),1]
	axis(1, at=xLab, labels=xLab, lwd=10, cex.axis=4.2, padj=1.2, tck=-.05)
	yLab=seq(0,ceiling(max(x[,-1])/1000)*1000, by=1000)
	axis(2, at=yLab, labels=yLab, lwd=10, cex.axis=4.2, padj=-.7)
	mtext("Year", side = 1, cex = 4.2, padj=3)
	mtext("Samples", side = 2, cex = 4.2, padj=-2.5)
	abline(h=0, col="#000000ff", lwd=5)
	plot.new()
	lPt = legPlot(colnames(x)[-1], ncol(x)-1)
	legend("top", legend=sub("[.]"," ",capFirst(colnames(x)[-1]))[lPt[[1]]], border=NA, ncol=lPt[[2]], lty=c(rep(1,length(1:(ncol(x)-1))),NA), title="Annual Review Medication Records", col=cBp[1:(ncol(x)-1)], lwd=5, cex=3.5, pch=c(1:(ncol(x)-1),NA)+20)
}
