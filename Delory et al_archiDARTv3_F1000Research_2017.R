#########################################################################################################################################
#R code used in Delory et al (2017) archiDART v3.0: a new data analysis pipeline allowing the topological analysis of plant root systems
#########################################################################################################################################

##########################################################################################
#To run the code properly, the latest version of archiDART must be installed (version 3.0)
#The following code can be used to install the version of archiDART used for this paper:
#devtools::install_github("archidart/archidart")
#Set your working directory
#Set the paths to RSML files
#and run the code! :-)
##########################################################################################

#Load R packages
library("archiDART")
library("vegan")
library("RColorBrewer")
library("ggplot2")
library("plotly")
library("FactoMineR")
library("grDevices")
library("ellipse")

#!!! Set working directory !!!

pathwd<-"SET_YOUR_WORKING_DIRECTORY_HERE"
setwd(pathwd)

#!!! Set paths to RSML files !!!

path1<-"SET_PATH_TO_FOLDER_CONTAINING_RSML_FOR_FIGURE_3"
path2<-"SET_PATH_TO_FOLDER_CONTAINING_100_RSML_LIBRARY_HERE"
path3<-"SET_PATH_TO_FOLDER_CONTAINING_4_RSML_FOR_FIGURE_5"

##########################################################
#Create plotly animation for geodesic distance (figure 3)
##########################################################

table<-data.frame()

for (geo in 0:44){
  
  table1<-rsmlToTable(inputrsml=path1, unitlength="cm", rsml.date="age", rsml.connect=T, fitter=FALSE)
  
  for (i in 1:nrow(table1)){
    
    if (table1$geodesic[i]>geo & (table1$geodesic[i]-table1$length[i])<geo){
      
      dist<-geo-(table1$geodesic[i]-table1$length[i])
      slope<-(table1$y2[i]-table1$y1[i])/(table1$x2[i]-table1$x1[i])
      invslope<-1/slope
      intercept<-table1$y2[i]-slope*table1$x2[i]
      ycoord<-table1$y1[i]+sqrt(dist^2/(1+invslope^2))
      xcoord<-(ycoord-intercept)/slope
      
      table1[nrow(table1)+1,]<-table1[i,]
      
      table1$y2[i]<-ycoord
      table1$y1[nrow(table1)]<-ycoord
      table1$x2[i]<-xcoord
      table1$x1[nrow(table1)]<-xcoord
      table1$length[i]<-dist
      table1$length[nrow(table1)]<-table1$length[nrow(table1)]-dist
      table1$geodesic[i]<-geo}}
  
  table1$geo<-45-geo
  table1$color[table1$geodesic>geo]<-"yes"
  table1$color[table1$geodesic<=geo]<-"no"
  
  table<-rbind(table, table1)}

minx<-min(table$x2)
miny<-min(table$y1)

table$y1<-table$y1-miny
table$y2<-table$y2-miny
table$x1<-table$x1-minx
table$x2<-table$x2-minx

ph<-perhomology(rsmlToTable(inputrsml=path1, unitlength="cm", rsml.date="age", rsml.connect=T, fitter=FALSE))
table1<-as.data.frame(ph$`ModelRootSystem_1`)
tableph<-data.frame()

for (geo in 0:44){
  
  table1<-as.data.frame(ph$`ModelRootSystem_1`)
  table1$endgeo<-rep(geo, nrow(table1))
  table1$begingeo<-rep(geo, nrow(table1))
  table1$geo<-45-geo
  table1$colorbegin<-rep("yes", nrow(table1))
  table1$colorend<-rep("no", nrow(table1))
  table1$hzero<-c(1:9)
  
  for (i in 1:nrow(table1)){
    
    if (geo>table1$birth[i]){
      table1$endgeo[i]<-table1$birth[i]
      table1$begingeo[i]<-table1$birth[i]}
    if (geo<table1$death[i]){
      table1$endgeo[i]<-table1$death[i]
      table1$begingeo[i]<-table1$death[i]}}
  
  tableph<-rbind(tableph, table1)}

tableph1<-tableph
n<-nrow(tableph1)*2
tableph<-data.frame(x1=rep(NA, n), y1=rep(NA, n), x2=rep(NA, n), y2=rep(NA, n), geo=rep(NA, n), color=rep(NA, n))

for (i in 1:nrow(tableph1)){
  
  tableph$x1[2*i-1]<-tableph1$birth[i]
  tableph$y1[2*i-1]<-tableph1$hzero[i]
  tableph$x2[2*i-1]<-tableph1$endgeo[i]
  tableph$y2[2*i-1]<-tableph1$hzero[i]
  tableph$geo[2*i-1]<-tableph1$geo[i]
  tableph$color[2*i-1]<-tableph1$colorbegin[i]
  
  tableph$x1[2*i]<-tableph1$begingeo[i]
  tableph$y1[2*i]<-tableph1$hzero[i]
  tableph$x2[2*i]<-tableph1$death[i]
  tableph$y2[2*i]<-tableph1$hzero[i]
  tableph$geo[2*i]<-tableph1$geo[i]
  tableph$color[2*i]<-tableph1$colorend[i]}

#Combine plots

col1<-"green4"

a <- list(
  text = "Model root system",
  font = list(size=18, color="black"),
  xref = "paper",
  yref = "paper",
  yanchor = "bottom",
  xanchor = "center",
  align = "center",
  x = 0.5,
  y = 1,
  showarrow = FALSE)

b <- list(
  text = "Persistence barcode",
  font = list(size=18, color="black"),
  xref = "paper",
  yref = "paper",
  yanchor = "bottom",
  xanchor = "center",
  align = "center",
  x = 0.5,
  y = 0.9,
  showarrow = FALSE)

gg1<-plot_ly(data=table, color=~color, x= ~x1, y=~-y1, xend=~x2, yend=~-y2, colors=c("gray70", col1)) %>%
  add_segments(line=list(width=4), frame=~geo) %>% 
  layout(annotations=a, xaxis=list(title="X (cm)", zeroline=FALSE, showline=TRUE, scaleanchor="y", gridcolor = toRGB("gray80"), range=c(-5,30), domain=c(0.1,0.9), dtick=5), yaxis=list(title="Y (cm)", zeroline=FALSE, showline=TRUE, gridcolor = toRGB("gray80"), range=c(-46,0)))

gg2<-plot_ly(data=tableph, color=~color, x= ~x1, y=~y1, xend=~x2, yend=~y2, colors=c("gray70", col1)) %>%
  add_segments(line=list(width=5), frame=~geo) %>% 
  layout(annotations=b, xaxis=list(title="Geodesic distance (cm)", zeroline=FALSE, showline=TRUE, gridcolor = toRGB("gray80"), range=c(46,0), domain=c(0.08,0.92), dtick=5), yaxis=list(title="H0", zeroline=FALSE, showline=TRUE, gridcolor = toRGB("gray80"), domain=c(0.1,0.9), range=c(9.5,0.5), dtick=1))

p<-subplot(gg1, gg2, nrows=1, titleX=TRUE, titleY=TRUE) %>%
  animation_opts(300, transition=0, easing="linear", redraw=FALSE) %>%
  hide_legend() %>%
  animation_slider(currentvalue = list(prefix = "Geodesic distance (cm): ", font = list(color="green", size=16)))

for (i in 1:45){p$x$layout$sliders[[1]]$steps[[i]]$label[1]<-45-i}

config(p, showLink=TRUE)

htmlwidgets::saveWidget(as_widget(p), "Figure3.html")

############################################
#Persistent homology on root system library
############################################

table<-rsmlToTable(inputrsml=path2, unitlength="cm", rsml.date="age", rsml.connect=T, fitter=FALSE, show.progress=TRUE)

archi<-architect(inputrsml=table, fitter=TRUE)

ph1<-perhomology(table, show.progress=TRUE)

dist1<-bottleneckdist(ph1, show.progress=TRUE)

########################################################
#Multidimensional scaling on bottleneck distance matrix
########################################################

mds3<-metaMDS(as.dist(dist1), trymax=1000, autotransform=FALSE, wascores = FALSE, noshare=FALSE, expand=FALSE)
type<-c(rep("taproot", 50),rep("fibrous", 50))
data3<-data.frame(mds3$points, Type=type)
data3$file<-rownames(data3)

stressplot(mds3)

########################################################################################################
#Principal component analysis on selected root systems using aggregated metrics calculated by architect
########################################################################################################

#Need to choose the root systems at the last observation date

files<-unique(archi$FileName)

datapca<-archi[1:100,]

for (i in 1:100){
  
  sub<-archi[archi$FileName==files[i],]
  datapca[i,]<-sub[nrow(sub),]}

rownames(datapca)<-c(1:100)

#PCA with Fitter indices

datapca1<-datapca[,c(3,5,7:13,18:25,28,32,36)]
rownames(datapca1)<-datapca$FileName
datapca1$Type<-c(rep("taproot", 50), rep("fibrous", 50))

pca<-PCA(datapca1, scale.unit=TRUE, ncp=3, quali.sup=21, graph=TRUE)
contrib<-pca$var$contrib
correlation<-pca$var$cor
ind<-data.frame(pca$ind$coord, Type=datapca1$Type)
colnames(ind)<-c("PC1", "PC2", "PC3", "Type")
ind$file<-rownames(ind)
ind<-ind[order(ind$file),]

#Select root systems (based on PCA and NMDS)

distpca<-as.matrix(dist(ind[,1:2])) #Distance matrix for PCA (traits)
distmds<-as.matrix(dist(data3[,1:2])) #Distance matrix for MDS (topology)
distindex<-(distpca-distmds)/(distpca+distmds) #Compute an index to detect root systems that cluster on the PCA map, but not on MDS map

mini<-min(as.dist(distindex[1:50, 1:50]))
maxi<-max(as.dist(distindex[1:50, 1:50]))
which(distindex==mini, arr.ind=T) #Select 2 dicot root systems: 108-10-25 and 102-8-25

mini<-min(as.dist(distindex[51:100, 51:100]))
maxi<-max(as.dist(distindex[51:100, 51:100]))
which(distindex==mini, arr.ind=T) #Select 2 root systems: 35-7-25 and 192-2-25

#################
#Create figure 5
#################

table<-rsmlToTable(inputrsml=path3, unitlength="cm", rsml.connect=TRUE, rsml.date="age", show.progress=TRUE)
table1<-table

table1$system<-table1$file
table1$system[table1$file=="monocot-sim-35-7-25"]<-"Fibrous 1"
table1$system[table1$file=="monocot-sim-192-2-25"]<-"Fibrous 2"
table1$system[table1$file=="dicot-sim-108-10-25"]<-"Taproot 1"
table1$system[table1$file=="dicot-sim-102-8-25"]<-"Taproot 2"

tiff(filename="Figure5.tif", res=1000, height=9, width=12, units="cm", compression="none", pointsize=6)

p<-ggplot(table1) + 
  geom_segment(aes(x = x1, y = -y1, xend = x2, yend = -y2, colour = geodesic), size = 0.5, alpha = 1) + 
  coord_fixed() + 
  theme_bw() + 
  facet_wrap(~system, ncol=4) +
  xlab("X (cm)") +
  ylab("Y (cm)") +
  scale_y_continuous(breaks=c(-50, -100, -150, -200, -250, -300, -350, -400)) +
  scale_colour_gradient(name = "Geodesic\ndistance (cm)", breaks=c(100, 200, 300)) +
  theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"), panel.spacing = unit(1.3, "lines"), legend.position=c(0.11,0.24), text=element_text(size=8), legend.text=element_text(size=7), axis.text=element_text(size=9), strip.text=element_text(size=9), legend.key.size=unit(1.4, "line"))
p

dev.off()

#################
#Create figure 4
#################

data3$pch<-c(rep(1, 50), rep(16,50))
ind$pch<-c(rep(1, 50), rep(16,50))
colors<-c("green4", "green4", "orangered", "orangered")

tiff(filename="Figure4.tif", res=1000, height=16, width=16, units="cm", compression="none", pointsize=9)
layout(matrix(1:4, ncol=2, nrow=2))

plot(ind$PC1, ind$PC2, xlab="PC1 (32.42%)", ylab="PC2 (24.58%)", type="n", las=1, main="PCA (individual map)", ylim=c(-7,9), xlim=c(-7,10))
abline(h=0, v=0, lty=2)
for (i in 1:100) {
  if (ind$file[i]!="monocot-sim-35-7-25_1" & ind$file[i]!="monocot-sim-192-2-25_1" & ind$file[i]!="dicot-sim-108-10-25_1" & ind$file[i]!="dicot-sim-102-8-25_1") {points(ind$PC1[i], ind$PC2[i], col="gray60", pch=ind$pch[i], cex=1)}}
points(ind$PC1[ind$file=="monocot-sim-35-7-25_1"], ind$PC2[ind$file=="monocot-sim-35-7-25_1"], col=colors[1], pch=3, cex=2.5, lwd=2)
points(ind$PC1[ind$file=="monocot-sim-192-2-25_1"], ind$PC2[ind$file=="monocot-sim-192-2-25_1"], col=colors[2], pch=4, cex=2.5, lwd=2)
points(ind$PC1[ind$file=="dicot-sim-108-10-25_1"], ind$PC2[ind$file=="dicot-sim-108-10-25_1"], col=colors[3], pch=3, cex=2.5, lwd=2)
points(ind$PC1[ind$file=="dicot-sim-102-8-25_1"], ind$PC2[ind$file=="dicot-sim-102-8-25_1"], col=colors[4], pch=4, cex=2.5, lwd=2)
ind1<-ind[ind$Type=="fibrous",]
ind2<-ind[ind$Type=="taproot",]
par(new=TRUE)
plot(ellipse(x=cor(ind1$PC1, ind1$PC2), scale=c(sd(ind1$PC1), sd(ind1$PC2)), centre=c(mean(ind1$PC1), mean(ind1$PC2)), level=0.95), type="l", col="forestgreen", ylim=c(-7,9), xlim=c(-7,10), axes=FALSE, lwd=1, xlab="", ylab="")
par(new=TRUE)
plot(ellipse(x=cor(ind2$PC1, ind2$PC2), scale=c(sd(ind2$PC1), sd(ind2$PC2)), centre=c(mean(ind2$PC1), mean(ind2$PC2)), level=0.95), type="l", col="orangered", ylim=c(-7,9), xlim=c(-7,10), axes=FALSE, lwd=1, xlab="", ylab="")
legend("topleft", legend=c("Fibrous 1", "Fibrous 2", "Taproot 1", "Taproot 2"), col=colors, pch=c(3,4,3,4), bty="n", pt.cex=1.3, cex=1)
legend("bottomleft", legend=c("Fibrous", "Taproot"), col=c("gray60", "gray60"), pch=c(16,1), bty="n", pt.cex=1.3, cex=1)
title("A", cex.main=1.7, font=2, adj=0)

plot(data3$MDS1, data3$MDS2, xlab="NMDS1", ylab="NMDS2", type="n", las=1, main="NMDS", ylim=c(-40,70), xlim=c(-150,100))
abline(h=0, v=0, lty=2)
for (i in 1:100) {
  if (data3$file[i]!="monocot-sim-35-7-25_1" & data3$file[i]!="monocot-sim-192-2-25_1" & data3$file[i]!="dicot-sim-108-10-25_1" & data3$file[i]!="dicot-sim-102-8-25_1") {points(data3$MDS1[i], data3$MDS2[i], col="gray60", pch=data3$pch[i], cex=1)}}
points(data3$MDS1[data3$file=="monocot-sim-35-7-25_1"], data3$MDS2[data3$file=="monocot-sim-35-7-25_1"], col=colors[1], pch=3, cex=2.5, lwd=2)
points(data3$MDS1[data3$file=="monocot-sim-192-2-25_1"], data3$MDS2[data3$file=="monocot-sim-192-2-25_1"], col=colors[2], pch=4, cex=2.5, lwd=2)
points(data3$MDS1[data3$file=="dicot-sim-108-10-25_1"], data3$MDS2[data3$file=="dicot-sim-108-10-25_1"], col=colors[3], pch=3, cex=2.5, lwd=2)
points(data3$MDS1[data3$file=="dicot-sim-102-8-25_1"], data3$MDS2[data3$file=="dicot-sim-102-8-25_1"], col=colors[4], pch=4, cex=2.5, lwd=2)
text(-160, 67, paste("Stress: ", round(mds3$stress, 4), sep=""), cex=1, pos=4)
data31<-data3[data3$Type=="fibrous",]
data32<-data3[data3$Type=="taproot",]
par(new=TRUE)
plot(ellipse(x=cor(data31$MDS1, data31$MDS2), scale=c(sd(data31$MDS1), sd(data31$MDS2)), centre=c(mean(data31$MDS1), mean(data31$MDS2)), level=0.95), type="l", col="forestgreen", ylim=c(-40,70), xlim=c(-150,100), axes=FALSE, lwd=1, xlab="", ylab="")
par(new=TRUE)
plot(ellipse(x=cor(data32$MDS1, data32$MDS2), scale=c(sd(data32$MDS1), sd(data32$MDS2)), centre=c(mean(data32$MDS1), mean(data32$MDS2)), level=0.95), type="l", col="orangered", ylim=c(-40,70), xlim=c(-150,100), axes=FALSE, lwd=1, xlab="", ylab="")
title("C", cex.main=1.7, font=2, adj=0)

plot(0,0,type="n", xlim=c(-1.2,1.2), ylim=c(-1,1), xlab="PC1 (32.42%)", ylab="PC2 (24.58%)", las=1, main="PCA (variable map)")
abline(v=0, h=0, lty=2)
for (i in 1:20){arrows(x0=0, y0=0, x1=correlation[i,1], y1=correlation[i,2], angle=30, length=0.07, col="gray60")}
var<-rownames(correlation)
pos<-c(4,2,1,4,1,1,2,4,2,4,3,2,3,1,3,1,4,3,2,3)
text(correlation[,1], correlation[,2], labels=var, pos=pos)
title("B", cex.main=1.7, font=2, adj=0)

plot(ph1$`dicot-sim-108-10-25_1`, xlab="Geodesic distance (cm)", xlim=c(250,0), ylim=c(89,0), las=1, main="Persistence barcodes", col=adjustcolor("orangered", alpha=0.6), lwd=2)
par(new=TRUE)
plot(ph1$`monocot-sim-35-7-25_1`, xlab="", xlim=c(250,0), ylim=c(89,0), las=1, main="", col=adjustcolor("forestgreen", alpha=0.5), axes=FALSE, lwd=2)
axis(2, at=c(10,30,50,70,90), las=1)

legend("bottomleft", legend=c("Fibrous 1", "Taproot 1"), col=c("forestgreen", "orangered"), lty=1, lwd=1, bty="n", cex=1)
title("D", cex.main=1.7, font=2, adj=0)

dev.off()