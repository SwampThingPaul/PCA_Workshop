## 
## PCA Example
##
## Code was compiled by Paul Julian
## contact info: pauljulianphd@gmail.com

# Libraries ---------------------------------------------------------------
#devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper)

library(reshape)
library(vegan)
library(REdaS)

#devtools::install_github("karthik/wesanderson")
library(wesanderson)

# Import Data -------------------------------------------------------------
dat<-read.csv("./data/lake_data.csv") 

# Data massage and clean-up -----------------------------------------------
# Cross tabulate the data based on parameter name
dat.xtab <- cast(dat,Station.ID+LAKE+Date.EST~param,value="HalfMDL",mean)

# Cleaning up/calculating parameters
dat.xtab$TN <- with(dat.xtab,TN_Combine(NOx,TKN,TN))
dat.xtab$DIN <- with(dat.xtab, NOx+NH4)

# More cleaning of the dataset 
vars <- c("Alk","Cl","Chla","DO","pH","SRP","TP","TN","DIN")
dat.xtab <- dat.xtab[,c("Station.ID","LAKE","Date.EST",vars)]

head(dat.xtab[,c("Station.ID","LAKE",vars)],4L)

#How many rows of data do we have?
nrow(dat.xtab)

#How many rows contain NAs
nrow(na.omit(dat.xtab))

dat.xtab <- na.omit(dat.xtab)


# scatterplots ------------------------------------------------------------
par(family="serif",mar=c(1,1.5,0.1,0.1),oma=c(3,3.5,0.75,0.5));
layout(matrix(1:72,8,8))

params=c(names(dat.xtab)[4:12])
axis.lab=c("Alk\n(mg L\u207B\u00B9)","Cl\n(mg L\u207B\u00B9)", "Chl-a\n(\u03BCg L\u207B\u00B9)","DO\n(mg L\u207B\u00B9)","pH\n(unitless)","SRP\n(mg L\u207B\u00B9)","TP\n(mg L\u207B\u00B9)", "TN\n(mg L\u207B\u00B9)","DIN\n(mg L\u207B\u00B9)")

for(j in 1:8){
  if(j!=1){for(k in 1:(j-1)){plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)}}
  
  params2=params[-1:-j]
  axis.lab2=axis.lab[-1:-j]
  xlim.val=c(0,175);by.x=50;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
  lim.min=c(0,0,0,0,5,0,0,0,0)
  lim.max=c(175,80,125,13,10,0.15,1.5,10,0.8);by.val=c(75,40,75,6,2.5,0.06,0.75,5,0.4)
  for(i in 1:length(params2)){
    xlim.val=c(lim.min[j],lim.max[j]);by.x=by.val[j];xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
    ylim.val=c(lim.min[-1:-j][i],lim.max[-1:-j][i]);by.y=by.val[-1:-j][i];ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y)
    plot(dat.xtab[,params[j]],dat.xtab[,params2[i]],xlim=xlim.val,ylim=ylim.val,axes=F,type="n",ylab=NA,xlab=NA)
    abline(h=ymaj,v=xmaj,lty=3,col="grey")
    points(dat.xtab[,params[j]],dat.xtab[,params2[i]],pch=21,bg=adjustcolor("dodgerblue1",0.25),col=adjustcolor("grey",0.5),lwd=0.01,cex=0.8)
    if(i==length(params2)){axis_fun(1,xmaj,xmin,format(xmaj),cex=0.8,line=-0.75)}else{axis_fun(1,xmaj,xmin,NA,cex=0.8,line=-0.5)}
    if(j==1){axis_fun(2,ymaj,ymin,format(ymaj),cex=0.8)}else{axis_fun(2,ymaj,ymin,NA,cex=0.8)}
    box(lwd=1)
    if(j==1){mtext(side=2,line=2.25,cex=0.8,axis.lab2[i])}
  }
  mtext(side=1,line=2.5,cex=0.8,axis.lab[j])
}


# PCA Assumption checks ---------------------------------------------------
# KMOS
KMOS(dat.xtab[,vars])

# Bartlett's Test Of Sphericity
bart_spher(dat.xtab[,vars])


# PCA ---------------------------------------------------------------------
dat.xtab.pca <- rda(dat.xtab[,vars],scale = T)

#dat.xtab.pca <- princomp(dat.xtab[,vars])

summary(dat.xtab.pca)$cont

#Extract eigenvalues (see definition above)
eig <- dat.xtab.pca$CA$eig

# Percent of variance explained by each compoinent
variance <- eig*100/sum(eig)

# The cumulative variance of each component (should sum to 1)
cumvar <- cumsum(variance)

# Combine all the data into one data.frame
eig.pca <- data.frame(eig = eig, variance = variance,cumvariance = cumvar)

# Not explored in presentation but discusses in blog post
# https://swampthingecology.org/blog/pca-basics-in-rstats/

#library(eigenprcomp)
#boot_pr_comp(as.matrix(dat.xtab[,vars]))

# Scree Plot --------------------------------------------------------------
layout(matrix(1:2,1,2))
par(family="serif",mar=c(1,2,0.1,1),oma=c(3,1,0.75,0.5));

ylim.val=c(0,5);by.y=1;ymaj=seq(ylim.val[1],100,by.y);ymin=seq(ylim.val[1],100,by.y/2)
x=barplot(eig.pca$eig,ylim=ylim.val,col="grey",yaxt="n")
abline(h=ymaj,lty=3,col="grey")
x=barplot(eig.pca$eig,ylim=ylim.val,col="grey",yaxt="n",add=T)
abline(h=1,lty=2,col="red",lwd=2)
axis_fun(1,line=-0.7,x,x,seq(1,length(x),1),0.7)
axis_fun(2,ymaj,ymin,ymaj,0.75);box(lwd=1)
mtext(side=1,line=1.5,"Principal Components")
mtext(side=2,line=1.5,"Eigenvalue")

ylim.val=c(0,110);by.y=25;ymaj=seq(ylim.val[1],100,by.y);ymin=seq(ylim.val[1],100,by.y/2);#set y limit and delineates the major and minor ticks
x=barplot(eig.pca$variance,ylim=ylim.val,col="white",border=0,yaxt="n")# inital plot to get the measurements
abline(h=ymaj,lty=3,col="grey")#makes vertical lines from y axis
x=barplot(eig.pca$variance,ylim=ylim.val,col="grey",yaxt="n",add=T)# the real plot that matters
lines(x,eig.pca$cumvariance,col="indianred1",lwd=2)# adds the cumulative variance for each factor
points(x,eig.pca$cumvariance,pch=21,bg="indianred1",cex=1.25)
abline(h=80,lty=2,col="red",lwd=2)
axis_fun(1,line=-0.7,x,x,seq(1,length(x),1),0.7)
axis_fun(2,ymaj,ymin,ymaj,0.75);box(lwd=1)
mtext(side=1,line=1.5,"Principal Components")
mtext(side=2,line=1.75,"Percentage of Variances")
legend.text=c("Absolute","Cumulative");#helper vaiable for legend
pt.col=c("grey","indianred1")#helper vaiable for legend
legend("topleft",legend=legend.text,pch=c(22,21),pt.bg=pt.col,col=c("black",pt.col[2]),lty=c(0,1),lwd=1.5,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend("topleft",legend=legend.text,pch=c(22,21),pt.bg=pt.col,col="black",lty=0,lwd=0.5,pt.cex=1.55,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)



# Score Extraction --------------------------------------------------------
scrs=scores(dat.xtab.pca,display=c("sites","species"),choices=c(1,2,3));


# Biplot ------------------------------------------------------------------
# version 1
par(family="serif",mar=c(1,3,0.1,0.5),oma=c(3,1.5,0.75,0.5));
layout(matrix(1:2,1,2))

xlim.val=c(-0.5,1);by.x=0.5;xmaj=c(0,seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2);
ylim.val=c(-0.5,1);by.y=0.5;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2);
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA);
abline(h=0,v=0,lty=3,col="grey");
points(scrs$sites[,c(1,2)],pch=21,bg="grey",cex=1,lwd=0.5); #plots the points
#arrows(0,0,scrs$species[,1],scrs$species[,2],length = 0.05, angle = 15, code = 2,col="indianred1",lwd=1.5);# makes the arrows
#with(scrs,text(species[,1]-0.1,species[,2],labels=rownames(species),cex=0.75));#adds labels to the arrows; 
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1); #adds x axis ticks
axis_fun(2,ymaj,ymin,format(ymaj),1); #adds y axis ticks
mtext(side=1,line=1.8,paste0("PCA 1 (",round(eig.pca$variance[1],1),"%)"));#adds x axis label with percent variance
mtext(side=2,line=2.25,paste0("PCA 2 (",round(eig.pca$variance[2],1),"%)"));#adds y axis label with percent variance


xlim.val=c(-1,1.5);by.x=0.5;xmaj=c(0,seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2);
ylim.val=c(-1.5,1.5);by.y=0.5;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2);
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA);
abline(h=0,v=0,lty=3,col="grey");
points(scrs$sites[,c(1,3)],pch=21,bg="grey",cex=1,lwd=0.5); 
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1); 
axis_fun(2,ymaj,ymin,format(ymaj),1); 
mtext(side=1,line=1.8,paste0("PCA 1 (",round(eig.pca$variance[1],1),"%)"));
mtext(side=2,line=2.25,paste0("PCA 3 (",round(eig.pca$variance[3],1),"%)"));

#version 2
par(family="serif",mar=c(1,3,0.1,0.5),oma=c(2.5,1.5,0.75,0.5));
layout(matrix(1:2,1,2))

labs=c("Alk","Cl","Chl-a","DO","pH","SRP","TP","TN","DIN")
xlim.val=c(-0.5,1);by.x=0.5;xmaj=c(0,seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2);
ylim.val=c(-0.5,1);by.y=0.5;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2);
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA);
abline(h=0,v=0,lty=3,col="grey");
points(scrs$sites[,c(1,2)],pch=21,bg=adjustcolor("grey",0.25),col=adjustcolor("grey",0.75),cex=1,lwd=0.01); #plots the points
arrows(0,0,scrs$species[,1]/3,scrs$species[,2]/3,length = 0.05, angle = 15, code = 2,col="indianred1",lwd=1.5);# makes the arrows
with(scrs,text((species[,1]+0.15)/3,species[,2]/3,labels=labs,cex=0.75,font=3));#adds labels to the arrows; 
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1); #adds x axis ticks
axis_fun(2,ymaj,ymin,format(ymaj),1); #adds y axis ticks
mtext(side=1,line=1.8,paste0("PCA 1 (",round(eig.pca$variance[1],1),"%)"));#adds x axis label with percent variance
mtext(side=2,line=2.25,paste0("PCA 2 (",round(eig.pca$variance[2],1),"%)"));#adds y axis label with percent variance


xlim.val=c(-1,1.5);by.x=0.5;xmaj=c(0,seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2);
ylim.val=c(-1.5,1.5);by.y=0.5;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2);
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA);
abline(h=0,v=0,lty=3,col="grey");
points(scrs$sites[,c(1,3)],pch=21,bg=adjustcolor("grey",0.25),col=adjustcolor("grey",0.75),cex=1,lwd=0.01); #plots the points
arrows(0,0,scrs$species[,1]/2,scrs$species[,3]/2,length = 0.05, angle = 15, code = 2,col="indianred1",lwd=1.5);# makes the arrows
with(scrs,text((species[,1]+0.15)/2,species[,3]/2,labels=labs,cex=0.75,font=3));#adds labels to the arrows; 
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1); 
axis_fun(2,ymaj,ymin,format(ymaj),1); 
mtext(side=1,line=1.8,paste0("PCA 1 (",round(eig.pca$variance[1],1),"%)"));
mtext(side=2,line=2.25,paste0("PCA 3 (",round(eig.pca$variance[3],1),"%)"));


##combine scores with orginal data
dat.xtab=cbind(dat.xtab,scrs$sites)

head(dat.xtab,3L)

#version 3
par(family="serif",mar=c(1.5,3,0.1,0.5),oma=c(0,1.5,0.75,0.5));
layout(matrix(c(1:2,3,3),2,2,byrow=T),heights=c(1,0.3))

#length(unique(dat.xtab$LAKE))
cols=wesanderson::wes_palette("Zissou1",6,"continuous")

labs=c("Alk","Cl","Chl-a","DO","pH","SRP","TP","TN","DIN")
xlim.val=c(-0.5,1);by.x=0.5;xmaj=c(0,seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2);
ylim.val=c(-0.5,1);by.y=0.5;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2);
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA);
abline(h=0,v=0,lty=3,col="grey");
for(i in 1:6){
  with(subset(dat.xtab,LAKE==unique(dat.xtab$LAKE)[i]),points(PC1,PC2,pch=21,bg=adjustcolor(cols[i],0.25),col=adjustcolor(cols[i],0.5),lwd=0.1,cex=1.25))
}
arrows(0,0,scrs$species[,1]/3,scrs$species[,2]/3,length = 0.05, angle = 15, code = 2,col="indianred1",lwd=1.5);# makes the arrows
with(scrs,text((species[,1]+0.15)/3,species[,2]/3,labels=labs,cex=0.75,font=3));#adds labels to the arrows; 
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1); #adds x axis ticks
axis_fun(2,ymaj,ymin,format(ymaj),1); #adds y axis ticks
mtext(side=1,line=1.8,paste0("PCA 1 (",round(eig.pca$variance[1],1),"%)"));#adds x axis label with percent variance
mtext(side=2,line=2.25,paste0("PCA 2 (",round(eig.pca$variance[2],1),"%)"));#adds y axis label with percent variance


xlim.val=c(-1,1.5);by.x=0.5;xmaj=c(0,seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2);
ylim.val=c(-1.5,1.5);by.y=0.5;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2);
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA);
abline(h=0,v=0,lty=3,col="grey");
for(i in 1:6){
  with(subset(dat.xtab,LAKE==unique(dat.xtab$LAKE)[i]),points(PC1,PC3,pch=21,bg=adjustcolor(cols[i],0.25),col=adjustcolor(cols[i],0.5),lwd=0.1,cex=1.25))
}
arrows(0,0,scrs$species[,1]/2,scrs$species[,3]/2,length = 0.05, angle = 15, code = 2,col="indianred1",lwd=1.5);# makes the arrows
with(scrs,text((species[,1]+0.15)/2,species[,3]/2,labels=labs,cex=0.75,font=3));#adds labels to the arrows; 
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1); 
axis_fun(2,ymaj,ymin,format(ymaj),1); 
mtext(side=1,line=1.8,paste0("PCA 1 (",round(eig.pca$variance[1],1),"%)"));
mtext(side=2,line=2.25,paste0("PCA 3 (",round(eig.pca$variance[3],1),"%)"));

plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)
legend(0.5,0,legend=unique(dat.xtab$LAKE),
       pch=c(21),
       col=adjustcolor(cols,0.5),
       pt.bg=adjustcolor(cols,0.25),
       lwd=c(0.01),lty=c(NA),pt.cex=1.5,ncol=3,cex=0.75,bty="n",y.intersp=1.75,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.1)


# END