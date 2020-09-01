#Limpieza de variables

rm (list = ls())

#Librerias

library(goseq)
library(limma)
library(edgeR)
library(statmod)

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
keep <- rowSums(d0$counts[,1:9])>3 & rowSums(d0$counts[,10:18])>3 & rowSums(d0$counts[,19:27])>3
d <- calcNormFactors(d0[keep,])

group <- factor(paste0(design$Strain, design$Temperature))
design2 <- model.matrix(~0+ group, data=design)

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

v <- voom(d, design2, plot = T, normalize = "quantile")

fit <- lmFit(v, design2)
fit_c <- contrasts.fit(fit,cont.matrix)
fit_c <- eBayes(fit_c)

contrastes<-c("ETOHRED12CvsCEMPK12C", "ETOHRED30CvsCEMPK30C", "ETOHRED39CvsCEMPK39C", "ETOHRED12CvsADY512C", "ETOHRED30CvsADY530C",
              "ETOHRED39CvsADY539C", "ADY512CvsCEMPK12C", "ADY530CvsCEMPK30C", "ADY539CvsCEMPK39C")

for (i in 1:length(contrastes)) {
  tabla<-topTable(fit_c, coef=i, number = length(fit_c$F.p.value), sort.by = "p")
  assign(contrastes[i], tabla)
}


make.dir("Genome_top")
make.dir("Genome_bottom")

for (i in 1:length(contrastes)) {
  
  data <- get(contrastes[i])
  data<-data[order(data$logFC, decreasing=TRUE),,drop=FALSE]
  
  gen_name_top<-rownames(data[data$logFC>1.5 & data$adj.P.Val<0.05,])
  
  gen_name_bottom<-rownames(data[data$logFC<(-1.5) & data$adj.P.Val<0.05,])
  
  file_top<-paste0("Genome_top/", contrastes[i],"_top")
  file_bottom<-paste0("Genome_bottom/",contrastes[i],"_bottom")
  
  write.table(gen_name_top, file=file_top, quote = FALSE, row.names= FALSE, col.names = FALSE)
  
  write.table(gen_name_bottom, file=file_bottom, quote = FALSE, row.names= FALSE, col.names = FALSE)

}

##########################################################################################################

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

#Script para el pangenoma

setwd("~/Documentos/TFM/Analisis_R_stranded")
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
keep <- rowSums(d0$counts[,1:9])>3 & rowSums(d0$counts[,10:18])>3 & rowSums(d0$counts[,19:27])>3
d <- calcNormFactors(d0[keep,])

group <- factor(paste0(design$Strain, design$Temperature))
design2 <- model.matrix(~0+ group, data=design)

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

v <- voom(d, design2, plot = T, normalize = "quantile")

fit <- lmFit(v, design2)
fit_c <- contrasts.fit(fit,cont.matrix)
fit_c <- eBayes(fit_c)

contrastes_p<-c("ETOHRED12CvsCEMPK12C_p", "ETOHRED30CvsCEMPK30C_p", "ETOHRED39CvsCEMPK39C_p", "ETOHRED12CvsADY512C_p", "ETOHRED30CvsADY530C_p",
                "ETOHRED39CvsADY539C_p", "ADY512CvsCEMPK12C_p", "ADY530CvsCEMPK30C_p", "ADY539CvsCEMPK39C_p")

for (i in 1:length(contrastes_p)) {
  tabla<-topTable(fit_c, coef=i, number = length(fit_c$F.p.value), sort.by = "p")
  assign(contrastes_p[i], tabla)
}


make.dir("Pangenome_top")
make.dir("Pangenome_bottom")

for (i in 1:length(contrastes_p)) {
  
  data <- get(contrastes_p[i])
  data<-data[order(data$logFC, decreasing=TRUE),,drop=FALSE]
  
  gen_name_top<-rownames(data[data$logFC>1.5 & data$adj.P.Val<0.05,])
  
  gen_name_bottom<-rownames(data[data$logFC<(-1.5) & data$adj.P.Val<0.05,])
  
  file_top<-paste0("Pangenome_top/", contrastes_p[i],"_top")
  file_bottom<-paste0("Pangenome_bottom/",contrastes_p[i],"_bottom")
  
  write.table(gen_name_top, file=file_top, quote = FALSE, row.names= FALSE, col.names = FALSE)
  
  write.table(gen_name_bottom, file=file_bottom, quote = FALSE, row.names= FALSE, col.names = FALSE)
  
}

