#Limpieza de variables

rm (list = ls())

#Librerias

library(limma)
library(edgeR)
library(statmod)
library(Glimma)
library("org.Sc.sgd.db")
library(gplots)

#Funciones

make.dir <- function(fp) {
  if(!file.exists(fp)) {  # If the folder does not exist, create a new one
    make.dir(dirname(fp))
    dir.create(fp)
  } else {   # If it existed, delete and replace with a new one  
    unlink(fp, recursive = TRUE)
    dir.create(fp)
  }
} 


#Script para el genoma

setwd("~/Documentos/TFM/Analisis_R_stranded")
files <- list.files(path = "./Genome/counts_genome")
filesT <- list.files(path = "./Genome/counts_trinity_genome")

datac<-NA
datat<-NA
for (i in 1:length(files)){
  tmp<-read.csv(paste("Genome/counts_genome/",files[i],sep=""),sep="\t",header=F)
  tmp2<-read.csv(paste("Genome/counts_trinity_genome/",filesT[i],sep=""),sep="\t",header=F)
  print(i)
  scnames<-tmp[,1]
  tnames<-tmp2[,1]
  datac<-cbind(datac,tmp[,2])
  datat<-cbind(datat,tmp2[,2])
  colnames(datac)[dim(datac)[2]]<-sub("_counts","",files[i])
  colnames(datat)[dim(datat)[2]]<-sub("_counts","",files[i])
}
row.names(datac)<-as.character(scnames)
row.names(datat)<-as.character(tnames)
datac <- datac[,-1]
datat <- datat[,-1]
datac<-rbind(datac,datat)
colnames(datac) <- gsub(".counts.table","",colnames(datac))
colnames(datac) <- gsub("-1-","-",colnames(datac))

#check distribution of reads per sample
barplot(apply(datac,2,sum,na.rm=T),las=2)

#experiment design
design<- as.matrix(read.csv2("layout.csv", sep ="\t"))
rownames(design) <- design[,1]
design <- design[,-1]
design <- as.data.frame(design)

#sort countdata colnames = design rownames
oo<-order(colnames(datac))
datac<-datac[,oo]
oo<-order(rownames(design))
design <- design[oo,]
cbind(rownames(design),colnames(datac))

#Elimina las reads no alineadas
malos <- grepl ("^_", rownames (datac))
datac[malos,]
datac <- datac[!malos,]

###Differential expression with limma

#define read count matrix
d0 <- DGEList(datac)
#remove genes with low expression (esto lo pongo aqqu? por defecto, 
#pero quiz? la mejor forma de eliminar genes con baja expresi?n sea
#quitar quellos genes con 0 lecturas en alguna de las tres cepas)
#Viendo el grafico voom al final el criterio 
#keep <- rowSums(cpm(d0)>2)>=3
#keep <- rowSums(d0$counts[,1:9])>0 & rowSums(d0$counts[,10:18])>0 & rowSums(d0$counts[,19:27])>0
keep <- rowSums(d0$counts[,1:9])>3 & rowSums(d0$counts[,10:18])>3 & rowSums(d0$counts[,19:27])>3
#normalization step from edgeR to scale the library sizes
d <- calcNormFactors(d0[keep,])
#experiment design for limma, defining group for variable interaction
#analysis
group <- factor(paste0(design$Strain, design$Temperature))
design2 <- model.matrix(~0+ group, data=design)
#define comparisons to make (aqu? solo pongo 2 ejemplos, habr?a que
#completar resto de temperaturas y cepas)
cont.matrix <- makeContrasts(
  "ETOHRED12CvsCEMPK12C"= groupETHANOLRED12C - groupCENPK12C,
  "ETOHRED30CvsCEMPK30C"=  groupETHANOLRED30C - groupCENPK30C,
  "ETOHRED39CvsCEMPK39C"=  groupETHANOLRED39C - groupCENPK39C,
  "ETOHRED12CvsADY512C"= groupETHANOLRED12C - groupADY512C,
  "ETOHRED30CvsADY530C"= groupETHANOLRED30C - groupADY530C,
  "ETOHRED39CvsADY539C"= groupETHANOLRED39C - groupADY539C,
  "ADY512CvsCEMPK12C"= groupADY512C - groupCENPK12C,
  "ADY530CvsCEMPK30C"= groupADY530C - groupCENPK30C,
  "ADY539CvsCEMPK39C"= groupADY539C - groupCENPK39C,
  levels = as.data.frame(design2))
#transform reads into logcounts per million, estimate mean-variance to compute
#appropiate weights at observation level
make.dir("voom_genome")
jpg<-"voom_genome/Voom_genome.jpg"
jpeg(jpg, width = 1500, height = 1500)
v <- voom(d, design2, plot = T, normalize = "quantile")
dev.off()

#linear model for every gene and smoothing
fit <- lmFit(v, design2)
fit_c <- contrasts.fit(fit,cont.matrix)
fit_c <- eBayes(fit_c)

make.dir("plotSA_genome")
jpg<-"plotSA_genome/plotSA_genome"
jpeg(jpg, width = 1500, height = 1500)
plotSA(fit_c, main="Final model: Mean-variance trend")
dev.off()

#obtain DE analysis
results <- decideTests(fit_c, adjust.method = "BH")

make.dir("vennDiagram_genome")
jpg<-"vennDiagram_genome/ETHvsCEM.jpg"
jpeg(jpg, width = 1500, height = 1500)
vennDiagram(results[,c("ETOHRED12CvsCEMPK12C", "ETOHRED30CvsCEMPK30C", "ETOHRED39CvsCEMPK39C")], cex = 2.5, circle.col=c("turquoise", "salmon", "blue"))
dev.off()
jpg<-"vennDiagram_genome/ETHvsADY.jpg"
jpeg(jpg, width = 1500, height = 1500)
vennDiagram(results[,c("ETOHRED12CvsADY512C", "ETOHRED30CvsADY530C", "ETOHRED39CvsADY539C")], cex = 2.5, circle.col=c("turquoise", "salmon", "blue"))
dev.off()
jpg<-"vennDiagram_genome/ADYvsCEM.jpg"
jpeg(jpg, width = 1500, height = 1500)
vennDiagram(results[,c("ADY512CvsCEMPK12C", "ADY530CvsCEMPK30C", "ADY539CvsCEMPK39C")], cex = 3, circle.col=c("turquoise", "salmon", "blue"))
dev.off()

summary(results)

#PCA clustering samples
make.dir("plotMDS_genome")
jpg<-"plotMDS_genome/plotMDS.jpg"
jpeg(jpg, width = 1500, height = 1500)
par(mar=c(6,7,3,1)+.1)
plotMDS(d, cex=4, cex.lab=3.2, cex.axis=2)
dev.off()

#See differentially expressed genes in a table with relevant info for every
#contrast

contrastes<-c("ETOHRED12CvsCEMPK12C", "ETOHRED30CvsCEMPK30C", "ETOHRED39CvsCEMPK39C", "ETOHRED12CvsADY512C", "ETOHRED30CvsADY530C",
              "ETOHRED39CvsADY539C", "ADY512CvsCEMPK12C", "ADY530CvsCEMPK30C", "ADY539CvsCEMPK39C")

for (i in 1:length(contrastes)) {
  tabla<-topTable(fit_c, coef=i, number = length(fit_c$F.p.value), sort.by = "p", p = 0.05)
  assign(contrastes[i], tabla)
}


#See differentially expressed genes in every contrast
make.dir("plotMD_genome")
for (i in 1:length(contrastes)) {
  jpg<-paste0("plotMD_genome/plotMD_",contrastes[i],".jpg")
  jpeg(jpg, width = 1500, height = 1500)
  par(mar=c(7,8,5,1)+.1)
  plotMD(fit_c,coef=i,status=results[,contrastes[i]], values = c(-1, 1), hl.cex=1, cex.lab=4, cex.main=4.5, legend=FALSE, cex.axis=2)
  dev.off()
}

#Creamos una tabla con anotaciones. Despues generamos un plotMD interactivo en html.
make.dir("glMDPlot_genome")
for (i in 1:length(contrastes)) {
  geneid<-rownames(results[,contrastes[i]])
  genes<-select(org.Sc.sgd.db, keys=geneid, columns=c("GENENAME","DESCRIPTION"), 
                keytype="ORF")
  folderg<-paste0("glMDPlot_genome/","glMDPlot_", contrastes[i])
  glMDPlot(fit_c, coef=i, status=results[,contrastes[i]], main = contrastes[i], anno = genes, groups=group, launch=FALSE, folder= folderg)
}

#Es una prueba y no si servira de algo. Se generan heatmaps.
lcpm<-cpm(d)
make.dir("heatmaps_genome")
for (j in 1:length(contrastes)) {
  top<-get(contrastes[j])[1:100,]
  i <- which(genes$ORF %in% rownames(top))
  mycol <- colorpanel(1000,"blue","white","red")
  jpg<-paste0("heatmaps_genome/",contrastes[j],".jpg")
  print(jpg)
  jpeg(jpg, width = 1500, height = 1500)
  heatmap.2(lcpm[i,], scale="row",
            labRow=genes$ORF[i], labCol=group, 
            col=mycol, trace="none", density.info="none", 
            margin=c(10,8), lhei=c(2,10), dendrogram="column")
  dev.off()
}

#See distribution of differential expression, and differentially expressed 
#gene names for every contrast
make.dir("volcanoplot_genome")
for (j in 1:length(contrastes)) {
  jpg<-paste0("volcanoplot_genome/",contrastes[j],".jpg")
  print(jpg)
  jpeg(jpg, width = 1500, height = 1500)
  par(mar=c(6,7,5,1)+.1)
  volcanoplot(fit_c,coef=j,highlight=150,names=labels(fit_c$Amean), main= contrastes[j], cex.lab=3, cex.main=3.7, cex.axis=2)
  dev.off()
}


############################################################################################################


#Limpieza de variables

rm (list = ls())

#Funciones

make.dir <- function(fp) {
  if(!file.exists(fp)) {  # If the folder does not exist, create a new one
    make.dir(dirname(fp))
    dir.create(fp)
  } else {   # If it existed, delete and replace with a new one  
    unlink(fp, recursive = TRUE)
    dir.create(fp)
  }
} 

#Script pangenoma

files <- list.files(path = "./Pangenome/counts_pangenome")
filesT <- list.files(path = "./Pangenome/counts_trinity_pangenome")

datac<-NA
datat<-NA
for (i in 1:length(files)){
  tmp<-read.csv(paste("Pangenome/counts_pangenome/",files[i],sep=""),sep="\t",header=F)
  tmp2<-read.csv(paste("Pangenome/counts_trinity_pangenome/",filesT[i],sep=""),sep="\t",header=F)
  print(i)
  scnames<-tmp[,1]
  tnames<-tmp2[,1]
  datac<-cbind(datac,tmp[,2])
  datat<-cbind(datat,tmp2[,2])
  colnames(datac)[dim(datac)[2]]<-sub("_counts","",files[i])
  colnames(datat)[dim(datat)[2]]<-sub("_counts","",files[i])
}
row.names(datac)<-as.character(scnames)
row.names(datat)<-as.character(tnames)
datac <- datac[,-1]
datat <- datat[,-1]
datac<-rbind(datac,datat)
colnames(datac) <- gsub(".counts.table","",colnames(datac))
colnames(datac) <- gsub("-1-","-",colnames(datac))

#check distribution of reads per sample
barplot(apply(datac,2,sum,na.rm=T),las=2)

#experiment design
design<- as.matrix(read.csv2("layout.csv", sep ="\t"))
rownames(design) <- design[,1]
design <- design[,-1]
design <- as.data.frame(design)

#sort countdata colnames = design rownames
oo<-order(colnames(datac))
datac<-datac[,oo]
oo<-order(rownames(design))
design <- design[oo,]
cbind(rownames(design),colnames(datac))

#Elimina las reads no alineadas
malos <- grepl ("^_", rownames (datac))
datac[malos,]
datac <- datac[!malos,]

###Differential expression with limma

#define read count matrix
d0 <- DGEList(datac)
#remove genes with low expression (esto lo pongo aqqu? por defecto, 
#pero quiz? la mejor forma de eliminar genes con baja expresi?n sea
#quitar quellos genes con 0 lecturas en alguna de las tres cepas)
#Viendo el grafico voom al final el criterio 
#keep <- rowSums(cpm(d0)>2)>=3
#keep <- rowSums(d0$counts[,1:9])>0 & rowSums(d0$counts[,10:18])>0 & rowSums(d0$counts[,19:27])>0
keep <- rowSums(d0$counts[,1:9])>3 & rowSums(d0$counts[,10:18])>3 & rowSums(d0$counts[,19:27])>3
#normalization step from edgeR to scale the library sizes
d <- calcNormFactors(d0[keep,])
#experiment design for limma, defining group for variable interaction
#analysis
group <- factor(paste0(design$Strain, design$Temperature))
design2 <- model.matrix(~0+ group, data=design)
#define comparisons to make (aqu? solo pongo 2 ejemplos, habr?a que
#completar resto de temperaturas y cepas)
cont.matrix <- makeContrasts(
  "ETOHRED12CvsCEMPK12C"= groupETHANOLRED12C - groupCENPK12C,
  "ETOHRED30CvsCEMPK30C"=  groupETHANOLRED30C - groupCENPK30C,
  "ETOHRED39CvsCEMPK39C"=  groupETHANOLRED39C - groupCENPK39C,
  "ETOHRED12CvsADY512C"= groupETHANOLRED12C - groupADY512C,
  "ETOHRED30CvsADY530C"= groupETHANOLRED30C - groupADY530C,
  "ETOHRED39CvsADY539C"= groupETHANOLRED39C - groupADY539C,
  "ADY512CvsCEMPK12C"= groupADY512C - groupCENPK12C,
  "ADY530CvsCEMPK30C"= groupADY530C - groupCENPK30C,
  "ADY539CvsCEMPK39C"= groupADY539C - groupCENPK39C,
  levels = as.data.frame(design2))
#transform reads into logcounts per million, estimate mean-variance to compute
#appropiate weights at observation level
make.dir("voom_pagenome")
jpg<-"voom_pagenome/Voom_pangenome.jpg"
jpeg(jpg, width = 1500, height = 1500)
v <- voom(d, design2, plot = T, normalize = "quantile")
dev.off()

#linear model for every gene and smoothing
fit <- lmFit(v, design2)
fit_c <- contrasts.fit(fit,cont.matrix)
fit_c <- eBayes(fit_c)

make.dir("plotSA_pangenome")
jpg<-"plotSA_pangenome/plotSA_pangenome"
jpeg(jpg, width = 1500, height = 1500)
plotSA(fit_c, main="Final model: Mean-variance trend")
dev.off()

#obtain DE analysis
results <- decideTests(fit_c, adjust.method = "BH")

make.dir("vennDiagram_pangenome")
jpg<-"vennDiagram_pangenome/ETHvsCEM.jpg"
jpeg(jpg, width = 1500, height = 1500)
vennDiagram(results[,c("ETOHRED12CvsCEMPK12C", "ETOHRED30CvsCEMPK30C", "ETOHRED39CvsCEMPK39C")], cex = 2.5, circle.col=c("turquoise", "salmon", "blue"))
dev.off()
jpg<-"vennDiagram_pangenome/ETHvsADY.jpg"
jpeg(jpg, width = 1500, height = 1500)
vennDiagram(results[,c("ETOHRED12CvsADY512C", "ETOHRED30CvsADY530C", "ETOHRED39CvsADY539C")], cex = 2.5, circle.col=c("turquoise", "salmon", "blue"))
dev.off()
jpg<-"vennDiagram_pangenome/ADYvsCEM.jpg"
jpeg(jpg, width = 1500, height = 1500)
vennDiagram(results[,c("ADY512CvsCEMPK12C", "ADY530CvsCEMPK30C", "ADY539CvsCEMPK39C")], cex = 3, circle.col=c("turquoise", "salmon", "blue"))
dev.off()

summary(results)

#PCA clustering samples
make.dir("plotMDS_pangenome")
jpg<-"plotMDS_pangenome/plotMDS.jpg"
jpeg(jpg, width = 1500, height = 1500)
par(mar=c(6,7,3,1)+.1)
plotMDS(d, cex=4, cex.lab=3.2, cex.axis=2)
dev.off()

#See differentially expressed genes in a table with relevant info for every
#contrast

contrastes<-c("ETOHRED12CvsCEMPK12C", "ETOHRED30CvsCEMPK30C", "ETOHRED39CvsCEMPK39C", "ETOHRED12CvsADY512C", "ETOHRED30CvsADY530C",
              "ETOHRED39CvsADY539C", "ADY512CvsCEMPK12C", "ADY530CvsCEMPK30C", "ADY539CvsCEMPK39C")

for (i in 1:length(contrastes)) {
  tabla<-topTable(fit_c, coef=i, number = length(fit_c$F.p.value), sort.by = "p", p = 0.05)
  assign(contrastes[i], tabla)
}


#See differentially expressed genes in every contrast
make.dir("plotMD_pangenome")
for (i in 1:length(contrastes)) {
  jpg<-paste0("plotMD_pangenome/plotMD_",contrastes[i],".jpg")
  jpeg(jpg, width = 1500, height = 1500)
  par(mar=c(7,8,5,1)+.1)
  plotMD(fit_c,coef=i,status=results[,contrastes[i]], values = c(-1, 1), hl.cex=1, cex.lab=4, cex.main=4.5, legend=FALSE, cex.axis=2)
  dev.off()
}

#Creamos una tabla con anotaciones. Despues generamos un plotMD interactivo en html.
make.dir("glMDPlot_pangenome")
for (i in 1:length(contrastes)) {
  geneid<-rownames(results[,contrastes[i]])
  genes<-select(org.Sc.sgd.db, keys=geneid, columns=c("GENENAME","DESCRIPTION"), 
                keytype="ORF")
  folderg<-paste0("glMDPlot_pangenome/","glMDPlot_", contrastes[i])
  glMDPlot(fit_c, coef=i, status=results[,contrastes[i]], main = contrastes[i], anno = genes, groups=group, launch=FALSE, folder= folderg)
}

#Es una prueba y no si servira de algo. Se generan heatmaps.
lcpm<-cpm(d)
make.dir("heatmaps_pangenome")
for (j in 1:length(contrastes)) {
  top<-get(contrastes[j])[1:100,]
  i <- which(genes$ORF %in% rownames(top))
  mycol <- colorpanel(1000,"blue","white","red")
  jpg<-paste0("heatmaps_pangenome/",contrastes[j],".jpg")
  print(jpg)
  jpeg(jpg, width = 1500, height = 1500)
  heatmap.2(lcpm[i,], scale="row",
            labRow=genes$ORF[i], labCol=group, 
            col=mycol, trace="none", density.info="none", 
            margin=c(10,8), lhei=c(2,10), dendrogram="column")
  dev.off()
}

#See distribution of differential expression, and differentially expressed 
#gene names for every contrast
make.dir("volcanoplot_pangenome")
for (j in 1:length(contrastes)) {
  jpg<-paste0("volcanoplot_pangenome/",contrastes[j],".jpg")
  print(jpg)
  jpeg(jpg, width = 1500, height = 1500)
  par(mar=c(5,6,4,1)+.1)
  volcanoplot(fit_c,coef=j,highlight=150,names=labels(fit_c$Amean), main= contrastes[j], cex.lab=3, cex.main=3.7, cex.axis=2)
  dev.off()
}

