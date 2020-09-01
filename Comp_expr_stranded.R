#Limpieza de variables

rm (list = ls())

#Librerias

library(tidyverse)
library(edgeR)
library(ggplot2)
library(broom)

#Funciones

TPM<-function(genes, long){
  rpk<-c()
  tpm<- matrix(ncol = 1, nrow = length(genes[,1]))
  colnames(tpm)<-"TPM"
  rownames(tpm)<-rownames(genes)
  
  for (i in 1:length(genes[,1])) {
    id<-rownames(genes)[i]
    pos<-which(row.names(long)==id)
    rpk<-c(rpk, genes[i,1]/(long[pos]/1000))
  }
  suma_rpk<-sum(rpk)/1000000
  for (i in 1:length(rpk)){
    tpm[i,1]<-rpk[i]/suma_rpk
  }
  tpm<-as.data.frame(tpm)
  return(tpm)
}

make.dir <- function(fp) {
  if(!file.exists(fp)) {  # If the folder does not exist, create a new one
    make.dir(dirname(fp))
    dir.create(fp)
  } else {   # If it existed, delete and replace with a new one  
    unlink(fp, recursive = TRUE)
    dir.create(fp)
  }
} 

#Script

setwd("~/Documentos/TFM/Analisis_R_stranded")
filesG <- list.files(path = "./Genome/counts_genome")
filesP <- list.files(path = "./Pangenome/counts_pangenome")

filename<-c()
#Se obtienen los nombres para cada muestra. Estos se usaran como variables.
for (i in 1:length(filesG)) {
  a<-str_split(filesG[i],"\\.", simplify = TRUE)
  b<-str_replace_all(a[1],"\\-","")
  c<-str_replace_all(b,"\\_", "")
  filename<-c(filename,c)
}

#Se leen todos los archivos y se asignan a los nombres creados anteriormente.
for (i in 1:length(filesG)){
  tmp<-read.csv(paste("Genome/counts_genome/",filesG[i],sep=""),sep="\t",header=F)
  tmp2<-read.csv(paste("Pangenome/counts_pangenome/",filesP[i],sep=""),sep="\t",header=F)
  colnames(tmp)<-c("genes", "reads")
  colnames(tmp2)<-c("genes", "reads")
  malos <- grepl ("^_", tmp$genes)
  tmp[malos,]
  tmp <- tmp[!malos,]
  
  malos <- grepl ("^_", tmp2$genes)
  tmp2[malos,]
  tmp2 <- tmp2[!malos,]
  
  
  tmp<-tmp[(tmp$genes %in% tmp2$genes),]
  tmp2<-tmp2[(tmp2$genes %in% tmp$genes),]
  tmp<-tmp[order(tmp$genes),]
  tmp2<-tmp2[order(tmp2$genes),]
  
  tmp3<- cbind(tmp[,2],tmp2[,2])
  rownames(tmp3)<-tmp[,1]
  colnames(tmp3)<-c("genoma", "pangenoma")
  
  assign(filename[i], tmp3)
}

#Se lee el archivo con las distancias de cada gen en el genoma.
long<-read.csv("gen_length.csv",sep="\t",header=F)
colnames(long)<-c("genes", "longitud")
long<-long[(long$genes %in% tmp$genes),]
long<-long[order(long$genes), , drop=FALSE]
long2<-as.matrix(long[,2])
rownames(long2)<-long[,1]
colnames(long2)<-"genes"

#Se lee el archivo con las distancias de cada gen en el pangenoma.
long_pan<-read.csv("gen_length_pan.csv",sep="\t",header=F)
colnames(long_pan)<-c("genes", "longitud")
long_pan<-long_pan[(long_pan$genes %in% tmp$genes),]
long_pan<-long_pan[order(long_pan$genes), , drop=FALSE]
long2_pan<-as.matrix(long_pan[,2])
rownames(long2_pan)<-long_pan[,1]
colnames(long2_pan)<-"genes"

#En lista_regresion se almacenara cuantos valores quedan por encima, por debajo o en la recta de regresion.
#En estadisiticos se almacenan los p valores y las R**2.
lista_regresion<-c()
estadisticos<-matrix(ncol = 2, nrow = length(filesG))

make.dir("GenomevsPangenoma")
make.dir("LogGenomevsLogPangenoma")
make.dir("TPMGenomevsTPMPangenoma")
make.dir("LogTPMGenomevsLogTPMPangenoma")
make.dir("ggplot_regresion")

for (i in 1:length(filesG)){
  #Los archivos se filtran y se convierten en dataframes.
  file<-get(filename[i])
  keep <- rowSums(file[,1, drop=FALSE])>0 | rowSums(file[,2, drop=FALSE])>0
  file<-as.data.frame(file[keep,])
  file<-file[order(file$genoma), c(1,2)]
  d <- DGEList(file)
  #Crea una representacion normal.
  jpg<-paste0("GenomevsPangenoma/",filename[i],"_expresion_comp.jpg")
  print(jpg)
  jpeg(jpg, width = 1500, height = 1500)
  plot(x=d$counts[,1], y=d$counts[,2], xlab = "genoma", ylab = "pangenoma")
  text(x=d$counts[,1], y=d$counts[,2], labels=rownames(d$counts), cex= 1, pos=3)
  dev.off()
  #Crea la representacion log-log.
  jpg<-paste0("LogGenomevsLogPangenoma/",filename[i],"_expresion_comp_log.jpg")
  print(jpg)
  jpeg(jpg, width = 1500, height = 1500)
  plot(x=log(d$counts[,1]), y=log(d$counts[,2]), xlab = "genoma", ylab = "pangenoma")
  text(x=log(d$counts[,1]), y=log(d$counts[,2]), labels=rownames(d$counts), cex= 1, pos=3)
  dev.off()
  assign(filename[i], d)
  
  #Se usa la funcion TPM y se representan los resultados tal cual y en log-log.
  gen<-TPM(d$counts[,1, drop=FALSE], long2)
  pan<-TPM(d$counts[,2, drop=FALSE], long2_pan)
  
  jpg<-paste0("TPMGenomevsTPMPangenoma/",filename[i],"_expresion_comp_norm.jpg")
  print(jpg)
  jpeg(jpg, width = 1500, height = 1500)
  plot(x=gen[,1], y=pan[,1], xlab = "genoma", ylab = "pangenoma")
  text(x=gen[,1], y=pan[,1], labels=rownames(d$counts), cex= 1, pos=3)
  dev.off()
  
  jpg<-paste0("LogTPMGenomevsLogTPMPangenoma/",filename[i],"_expresion_comp_log_norm.jpg")
  print(jpg)
  jpeg(jpg, width = 1500, height = 1500)
  plot(x=log(gen[,1]), y=log(pan[,1]), xlab = "genoma", ylab = "pangenoma")
  text(x=log(gen[,1]), y=log(pan[,1]), labels=rownames(d$counts), cex= 1, pos=3)
  dev.off()
  assign(filename[i], d)
  #Se crea el modelo de regresion y se calcula su intervalo de confianza.
  
  log_pan<-log(pan)
  log_gen<-log(gen)
  log_genes<-cbind(log_gen, log_pan)
  colnames(log_genes)<-c("genoma", "pangenoma")
  keep<-log_genes[,1]!=-Inf & log_genes[,2]!=-Inf
  log_genes<-log_genes[keep,]
  
  lm.model <- lm(log_genes[,2] ~ log_genes[,1])
  conf_interval <- predict(lm.model, interval="confidence", level = 0.95)
  #Se observan que valores estan por debajo, por encima o en la recta.
  matriz<-matrix(nrow = length(log_genes[,2]), ncol=1)
  rownames(matriz)<-rownames(log_genes)
  
  for (j in 1:length(log_genes[,2])) {
    if (log_genes[j,2]>conf_interval[j,3]) {
      matriz[j,1]<-"up"
    }
    else if (log_genes[j,2]<conf_interval[j,2]) {
      matriz[j,1]<-"down"
      
    } 
    else{
      matriz[j,1]<-"in"
    }
    
  }
  
  regresion<-paste0(filename[i],"_regresion")
  assign(regresion, summary(matriz))
  
  lista_regresion<-c(lista_regresion,regresion)
  #Se calcula el intervalo de prediccion y se representa junto con el intervalo de confianza.
  pred_interval <- predict(lm.model, interval="prediction", level = 0.95)
  file2<-cbind(log_genes,pred_interval)
  colnames(file2)<-c("genoma", "pangenoma", "fit", "lwr", "upr")
  
  plot_regresion<-paste0("ggplot_regresion/",filename[i],"_regresion.png")
  p<-ggplot(file2, aes(x=genoma, y=pangenoma))+
    geom_point()+
    stat_smooth(method=lm)
  p + geom_line(aes(y = lwr), color = "red", linetype = "dashed")+
    geom_line(aes(y = upr), color = "red", linetype = "dashed")
  p + theme(text=element_text(size=16))
  
  ggsave(plot_regresion)
  
  #Se guardan los estadisticos en una tabla.
  colnames(estadisticos)<-c("p_value", "R_squared")
  rownames(estadisticos)<-filename
  estadisticos[i,1]<- glance(lm.model)$p.value
  estadisticos[i,2]<- glance(lm.model)$r.squared

}

#Se obtiene la media de los valores en relacion al intervalo de confianza.
down<-c()
in_list<-c()
up<-c()

for (i in lista_regresion) {
  dato<-get(i)
  down<-c(down, as.numeric(str_split(dato[1],":", simplify = TRUE)[2]))
  in_list<-c(in_list, as.numeric(str_split(dato[2],":", simplify = TRUE)[2]))
  up<-c(up, as.numeric(str_split(dato[3],":", simplify = TRUE)[2]))
}

mean(down)
mean(in_list)
mean(up)

#Se muestran los estadisticos.
estadisticos

write.table(estadisticos, file="estadisticos_comp_expresion.csv", sep="\t", col.names = NA)
