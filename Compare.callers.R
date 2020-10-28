#!/bin/R
# Author: Ismael
# Date: 10-2020

# Comparacion de variantes detectadas por:
# Strelka2 ILUMINA
# VarScan2
# Mutect2 GATK
# Instalacion de los paquetes necesarios
#packages <- c("xlsx", "gdata", "plyr", "png")
#if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
#  install.packages(setdiff(packages, rownames(installed.packages())))  
#}
#lapply(packages, require, character.only=T)

#install.packages("tidyverse")
library(tidyr)
library(plyr)
library(gdata)
library(png)

directory = c("~/.temporal/output/")
filesMutect = list.files(directory, pattern = "Mutect")
filesStrelka = list.files(directory, pattern = "Strelka")
filesVarScan = list.files(directory, pattern = "VarScan")

setwd(directory)
mutect2 = read.table(filesMutect[1], sep = "\t", header = T)
strelka2 = read.table(filesStrelka[1], sep = "\t", header = T)
varscan2 = read.table(filesVarScan[1], sep = "\t", header = T)

mutect2fil = mutect2[,c(1:22,35,34,42,41)]
strelka2fil = strelka2[,c(1:22,51,63,38,64)]
mutect2fil$sampleID = mutect2$sample
mutect2fil$AF = toString(mutect2fil$AF)
mutect2fil$gAF = toString(mutect2fil$gAF)
strelka2fil$sampleID =strelka2$SampleID
mutect2fil$Caller = c("Mutect2")
varscan2$Caller = c("Varscan2")
strelka2fil$Caller = c("Strelka2")

colnames(strelka2fil) = c(colnames(mutect2fil))

var = rbind(mutect2fil, strelka2fil)
var <- var[with(var, order(var$Start)), ]

varscan2$NAs = c(NA)
varscan2fil = varscan2[,c(1,2,2,3,4,30,27,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,26,11,25,7,28,30)]
colnames(varscan2fil) = colnames(var)
varscan2fil$Caller = c("Varscan2")

var = rbind(var, varscan2fil)
var <- var[with(var, order(var$Start, var$Alt) ), ]
setwd(directory)
write.table(var, "Variantes.Strelka.Varscan.Mutect.validation", sep = "\t", row.names = F)

k = c()

for (l in 2:nrow(var)) {
  x = l-1
  if (var$Start[l] == var$Start[x] && 
      var$Chr[l] == var$Chr[x] && 
      var$Alt[l] == var$Alt[x] &&
      var$Ref[l] == var$Ref[x] &&
      var$Caller[l] != var$Caller[x]) 
    {
    k = c(k, l, x)
  }
  }

var2 = var[unique(k),]
var2 <- var2[with(var2, order(var2$Start, var2$Alt) ),]

setwd(directory)
write.table(var2, "Variantes.Strelka.Varscan.Mutect.validation", sep = "\t", row.names = F)



# ------------------------------------------------



# Comparacion y actualizacion de base de datos

# Creando la funcion encargada de realizar la comparacion con la base de datos
confronti.multipli <- function(mut.file, path, out.file){
  path <- path
  filenames <- list.files(path, pattern = "")
  mut <<- read.table(mut.file, sep = "\t", header = T)
#  colnames(mut)[34] = "FREQ_AL"
  n = ncol(mut)
  for (j in filenames){
    print(filenames)
    confr <<- read.table(paste(path, j, sep = ""), sep = "\t", header = T)
    mut <- join(mut, confr[, c("Chr", "Start", "Ref", "Alt", "FREQ_AL")], 
                by = c("Chr", "Start", "Ref", "Alt"))
    colnames(mut)[ncol(mut)] <- strsplit(as.character(j), "\\.")[[1]][1]
  }
  
  mut[,c((n+1):ncol(mut))][is.na(mut[,c((n+1):ncol(mut))])] <- 0
  write.table(mut, out.file, row.names = F, dec = ".", sep = "\t")
  cat("Comparacion con base de datos correcta!")
}

# Lanzando la funcion en una sola muestra
confronti.multipli("~/.temporal/output/Variantes.Strelka.Varscan.Mutect.validation",  # Input file
                   "~/.temporal/DB/",  # Database compare
                   "~/.temporal/output/Variantes.Strelka.Varscan.Mutect.comparedDB.validation")


