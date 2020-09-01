#Limpieza de variables

rm (list = ls())

#Librerias

library(limma)
library(edgeR)
library(tidyverse)
library(broom)

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

#Script_genoma

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

design<- as.matrix(read.csv2("layout.csv", sep ="\t"))
rownames(design) <- design[,1]
design <- design[,-1]
design <- as.data.frame(design)

oo<-order(colnames(datac))
datac<-datac[,oo]
oo<-order(rownames(design))
design <- design[oo,]
cbind(rownames(design),colnames(datac))

malos <- grepl ("^_", rownames (datac))
datac[malos,]
datac <- datac[!malos,]

###Differential expression with limma

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
  assign(contrastes[i], tabla[,c(1,5),drop=FALSE])
}

########################################################################################################################3

#Script_pangenoma

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

design<- as.matrix(read.csv2("layout.csv", sep ="\t"))
rownames(design) <- design[,1]
design <- design[,-1]
design <- as.data.frame(design)

oo<-order(colnames(datac))
datac<-datac[,oo]
oo<-order(rownames(design))
design <- design[oo,]
cbind(rownames(design),colnames(datac))

malos <- grepl ("^_", rownames (datac))
datac[malos,]
datac <- datac[!malos,]

###Differential expression with limma

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
  assign(contrastes_p[i], tabla[,c(1,5),drop=FALSE])
}

############################################################################################################################

#Regresion linear

make.dir("logFC_reg_plots")

lista_conteos<-c()
lista_regresion<-c()
lista_regresion_p_value<-c()
estadisticos<-matrix(ncol = 9, nrow = length(contrastes))
colnames(estadisticos)<-c("p_value", "R_squared", "no_significativo", "solo_pangenoma", "solo_genoma", "significativo", "down", "in", "up")
rownames(estadisticos)<-contrastes
estadisticos_pvalue<-matrix(ncol = 5, nrow = length(contrastes))
colnames(estadisticos_pvalue)<-c("p_value", "R_squared", "down", "in", "up")
rownames(estadisticos_pvalue)<-contrastes

for (i in 1:length(contrastes_p)) {
  
  tabla<-get(contrastes[i])
  tabla_p<-get(contrastes_p[i])
  
  tmp<- tabla[rownames(tabla) %in% rownames(tabla_p),]
  tmp2<- tabla_p[rownames(tabla_p) %in% rownames(tabla),]
  tmp<-tmp[order(rownames(tmp)),]
  tmp2<-tmp2[order(rownames(tmp2)),]
  g_tabla<-cbind(tmp,tmp2)
  colnames(g_tabla)<-c("logFC_genoma", "adj.P.Val_genoma", "logFC_pangenoma", "adj.P.Val_pangenoma")
  
  lm.model <- lm(g_tabla[,3] ~ g_tabla[,1])
  
  conf_interval <- predict(lm.model, interval="confidence", level = 0.95)
  matriz<-matrix(nrow = length(g_tabla[,3]), ncol=1)
  rownames(matriz)<-rownames(g_tabla)
  
  for (j in 1:length(g_tabla[,3])) {
    if (g_tabla[j,3]>conf_interval[j,3]) {
      matriz[j,1]<-"up"
    }
    else if (g_tabla[j,3]<conf_interval[j,2]) {
      matriz[j,1]<-"down"
      
    } 
    else{
      matriz[j,1]<-"in"
    }
    
  }
  
  regresion<-paste0(contrastes[i],"_regresion")
  assign(regresion, summary(matriz))
  
  lista_regresion<-c(lista_regresion,regresion)
  
  jpg<-paste0("logFC_reg_plots/",contrastes[i],"_logFC_reg.jpg")
  print(jpg)
  jpeg(jpg, width = 2500, height = 2500)
  par(mar=c(6,7,4,1)+.1)
  plot(x=g_tabla[,1], y=g_tabla[,3], cex = 4, pch = 16, cex.lab=6, xlab = "genoma_logFC", ylab = "pangenoma_logFC" ,col = ifelse (g_tabla[,2]>0.05 & g_tabla[,4] > 0.05, "black", ifelse(g_tabla[,2]<0.05 & g_tabla[,4] > 0.05, "red1", ifelse(g_tabla[,2]>0.05 & g_tabla[,4] < 0.05, "orange1", "green1"))))
  #text(x=g_tabla[,1], y=g_tabla[,3], labels=rownames(g_tabla), cex= 0.5, pos=3)
  abline(lm.model, col= rgb(0,0,0,alpha=0.3) )
  lines(x=g_tabla[,1], y=conf_interval[,3], col="grey1")
  lines(x=g_tabla[,1], y=conf_interval[,2], col="grey1")
  
  dev.off()
  
  estadisticos[i,1]<- glance(lm.model)$p.value
  estadisticos[i,2]<- glance(lm.model)$r.squared
   
  conteos<-matrix(ncol=1, nrow=length(g_tabla[,1]))
  rownames(conteos)<-rownames(g_tabla)
  nombre<-paste0(contrastes[i],"_conteos")
  
  for (j in 1:length(g_tabla[,1])) {
    if (g_tabla[j,2] > 0.05 & g_tabla[j,4] > 0.05) {
      conteos[j,1]<-"no significativo"
    }
    else if (g_tabla[j,2] > 0.05 & g_tabla[j,4] < 0.05) {
      conteos[j,1]<-"significativo en el pangenoma"
    }
    else if (g_tabla[j,2] < 0.05 & g_tabla[j,4] > 0.05) {
      conteos[j,1]<-"significativo en el genoma"
    }
    else if (g_tabla[j,2] < 0.05 & g_tabla[j,4] < 0.05) {
      conteos[j,1]<-"significativo en ambos"
    }
  }
  assign(nombre, conteos)
  
  lista_conteos<-c(lista_conteos,nombre)
  
  lm.model <- lm(log(g_tabla[,4]) ~ log(g_tabla[,2]))
  
  conf_interval <- predict(lm.model, interval="confidence", level = 0.95)
  matriz2<-matrix(nrow = length(g_tabla[,4]), ncol=1)
  rownames(matriz2)<-rownames(g_tabla)
  
  for (j in 1:length(g_tabla[,4])) {
    if (log(g_tabla[j,4])>conf_interval[j,3]) {
      matriz2[j,1]<-"up"
    }
    else if (log(g_tabla[j,4])<conf_interval[j,2]) {
      matriz2[j,1]<-"down"
      
    } 
    else{
      matriz2[j,1]<-"in"
    }
    
  }
  
  regresion<-paste0(contrastes[i],"_regresion_pvalue")
  assign(regresion, summary(matriz2))
  lista_regresion_p_value<-c(lista_regresion_p_value,regresion)
  
  jpg<-paste0("logFC_reg_plots/",contrastes[i],"_logFC_reg_pvalue.jpg")
  print(jpg)
  jpeg(jpg, width = 2500, height = 2500)
  plot(x=log(g_tabla[,2]), y=log(g_tabla[,4]), cex = 1.2, pch = 16, xlab = "genoma_pvalue", ylab = "pangenoma_value",col = ifelse (g_tabla[,2]>0.05 & g_tabla[,4] > 0.05, "black", ifelse(g_tabla[,2]<0.05 & g_tabla[,4] > 0.05, "red1", ifelse(g_tabla[,2]>0.05 & g_tabla[,4] < 0.05, "orange1", "green1"))) )
  #text(x=g_tabla[,2], y=g_tabla[,4], labels=rownames(g_tabla), cex= 0.5, pos=3)
  abline(lm.model, col= rgb(0,0,0,alpha=0.3) )
  lines(x=log(g_tabla[,2]), y=conf_interval[,3], col="grey1")
  lines(x=log(g_tabla[,2]), y=conf_interval[,2], col="grey1")
  
  dev.off()
  
  estadisticos_pvalue[i,1]<- glance(lm.model)$p.value
  estadisticos_pvalue[i,2]<- glance(lm.model)$r.squared
  
}


no_significativo<-c()
solo_pangenoma<-c()
solo_genoma<-c()
significativo<-c()

for (i in lista_conteos) {
  dato<-summary(get(i))
  no_significativo<-c(no_significativo, as.numeric(str_split(dato[1],":", simplify = TRUE)[2]))
  solo_pangenoma<-c(solo_pangenoma, as.numeric(str_split(dato[4],":", simplify = TRUE)[2]))
  solo_genoma<-c(solo_genoma, as.numeric(str_split(dato[3],":", simplify = TRUE)[2]))
  significativo<-c(significativo, as.numeric(str_split(dato[2],":", simplify = TRUE)[2]))
}

mean(no_significativo)
mean(solo_pangenoma)
mean(solo_genoma)
mean(significativo)

for (i in 1:length(no_significativo)) {
  estadisticos[i,3]<-no_significativo[i]
  estadisticos[i,4]<-solo_pangenoma[i]
  estadisticos[i,5]<-solo_genoma[i]
  estadisticos[i,6]<-significativo[i]
}

down<-c()
in_list<-c()
up<-c()

for (i in lista_regresion) {
  dato<-get(i)
  down<-c(down, as.numeric(str_split(dato[1],":", simplify = TRUE)[2]))
  in_list<-c(in_list, as.numeric(str_split(dato[2],":", simplify = TRUE)[2]))
  up<-c(up, as.numeric(str_split(dato[3],":", simplify = TRUE)[2]))
}

for (i in 1:length(down)) {
  estadisticos[i,7]<-down[i]
  estadisticos[i,8]<-in_list[i]
  estadisticos[i,9]<-up[i]
}

mean(down)
mean(in_list)
mean(up)

estadisticos

write.table(estadisticos, file="estadisticos_logFC.csv", sep="\t", col.names = NA)

down_p<-c()
in_list_p<-c()
up_p<-c()

for (i in lista_regresion_p_value) {
  dato<-get(i)
  down_p<-c(down_p, as.numeric(str_split(dato[1],":", simplify = TRUE)[2]))
  in_list_p<-c(in_list_p, as.numeric(str_split(dato[2],":", simplify = TRUE)[2]))
  up_p<-c(up_p, as.numeric(str_split(dato[3],":", simplify = TRUE)[2]))
}

for (i in 1:length(down_p)) {
  estadisticos_pvalue[i,3]<-down_p[i]
  estadisticos_pvalue[i,4]<-in_list_p[i]
  estadisticos_pvalue[i,5]<-up_p[i]
}

mean(down_p)
mean(in_list_p)
mean(up_p)

estadisticos_pvalue

write.table(estadisticos_pvalue, file="estadisticos_p_value.csv", sep="\t", col.names = NA)
