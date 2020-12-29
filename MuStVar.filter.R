#!/bin/R
# Author: Ismael
# Date: 10-2019

# Mutect2

#install.packages("tidyverse")
library(tidyr)

bed = c("~/.temporal/Beato_OncoPanel_TE-91860463_hg38_mpileup.bed")
bedfile = read.table(bed, sep = "\t", header = F, quote= "")

data_base = c("~/.temporal/DB/")

final_output = c("~/.temporal/output/")
directorio = c("~/.temporal/Mutect2/")
output = c("~/.temporal/Mutect2/filter/")
cat("\nFiltering the Mutect2 annotation files...\n")
setwd(directorio)
files_1 = list.files(directorio, pattern = "_completa.txt")
files_2 = list.files(directorio, pattern = "_filtrada.txt")
t = NULL

for (f in 1:length(files_1)) {
  setwd(directorio)
  table1 = read.table(files_1[f], sep = "\t", header=T, quote="")
  setwd(directorio)
  table2 = read.table(files_2[f], sep = "\t")
  tableV12 = table2[(ncol(table2)-1)]
  table3 = separate(tableV12, 1, c("GT", "AD", "AF", "DP", "F1R2", "F2R1", "SB"), sep = ":")
  tableV32 = table3[2]
  table32 = separate(tableV32, 1, c("ADref", "ADvar1","ADvar2", "ADvar3","ADvar4", "ADvar5","ADvar6", "ADvar7","ADvar8", "ADvar9", "ADvar10"), sep = ",")
  tableV13 = table2[ncol(table2)]
  table4 = separate(tableV13, 1, c("gGT", "gAD", "gAF", "gDP", "gF1R2", "gF2R1", "gSB"), sep = ":")
  table5 = cbind(table2[1],table2[2],table2[3],table2[4],table2[5],table2[6],table2[7],table2[8],table2[9])
  names(table5)=c(".","pos","chr","posCHR",".","ref","alt",".","Q")
  table1= cbind(table1, table5, table3, table4, table32)
  table1$AF1 = c(0)
  table1$AF2 = c(0)
  table1$AF3 = c(0)
  table1$AF4 = c(0)
  table1$AF5 = c(0)
  table1$AF6 = c(0)
  table1$AF7 = c(0)
  table1$AF8 = c(0)
  table1$AF9 = c(0)
  table1$AF10 = c(0)
  
  x = c()
  y = c()
  z = c()
  for (i in 1:nrow(table1)) {
    if (table1$Func.refGene[i] == "exonic") {x = rbind(x, table1[i,])}
    else if (table1$Func.refGene[i] == "splicing") {x <- rbind(x, table1[i,])}
    else if (table1$Func.refGene[i] == "UTR5") {x <- rbind(x, table1[i,])}
    else if (table1$Func.refGene[i] == "UTR3") {x <- rbind(x, table1[i,])}}
  
  for (i in 1:nrow(table1)) {
    if (toString(x$Gene.refGene[i]) == "MUC1") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "MUC2") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "MUC3") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "MUC4") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "POTEJ") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "POTEM") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "POTEE") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "ZNF91") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "NOTCH2NLB") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "NOTCH2NLC") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "FRMPD3") {z = c(z, i)}}
  if (toString(is.null(z)) == "FALSE") {x = x[-c(z),]}
  
  genes = unique(as.vector(bedfile$V4))
  n = c()
  
  for (i in 1:nrow(x)) {
    gene = toString(x$Gene.refGene[i])
    p = which(genes %in% gene)
    if (length(p) == 1) {n = c(n, i)}
  }
  x = x[c(n),]
  
  z = NULL
  if (is.null(x) == FALSE) {for (i in 1:nrow(x)) {
    if (as.numeric(x$gDP[i]) > 29 && as.numeric(x$DP[i]) > 499) {y = rbind(y, x[i,])}}}
  
  BD = y
  
  x = c()
  
  for (l in 1:nrow(y)) {
    y$AF1 = as.numeric(y$ADvar1)/as.numeric(y$DP); if (y$AF1[l] > 0.005) {x = c(l, x)}
    y$AF2 = as.numeric(y$ADvar2)/as.numeric(y$DP); if (is.na(y$AF2[l]) == F && y$AF2[l] > 0.005) {x = c(l, x)}
    y$AF3 = as.numeric(y$ADvar3)/as.numeric(y$DP); if (is.na(y$AF3[l]) == F && y$AF3[l] > 0.005) {x = c(l, x)}
    y$AF4 = as.numeric(y$ADvar4)/as.numeric(y$DP); if (is.na(y$AF4[l]) == F && y$AF4[l] > 0.005) {x = c(l, x)}
    y$AF5 = as.numeric(y$ADvar5)/as.numeric(y$DP); if (is.na(y$AF5[l]) == F && y$AF5[l] > 0.005) {x = c(l, x)}
    y$AF6 = as.numeric(y$ADvar6)/as.numeric(y$DP); if (is.na(y$AF6[l]) == F && y$AF6[l] > 0.005) {x = c(l, x)}
    y$AF7 = as.numeric(y$ADvar7)/as.numeric(y$DP); if (is.na(y$AF7[l]) == F && y$AF7[l] > 0.005) {x = c(l, x)}
    y$AF8 = as.numeric(y$ADvar8)/as.numeric(y$DP); if (is.na(y$AF8[l]) == F && y$AF8[l] > 0.005) {x = c(l, x)}
    y$AF9 = as.numeric(y$ADvar9)/as.numeric(y$DP); if (is.na(y$AF9[l]) == F && y$AF9[l] > 0.005) {x = c(l, x)}
    y$AF10 = as.numeric(y$ADvar10)/as.numeric(y$DP); if (is.na(y$AF10[l]) == F && y$AF10[l] > 0.005) {x = c(l, x)}}
  
  x = unique(x)
  y = y[c(x),]
  
  sampleID = gsub("_somatic.hg38_multianno.txt_completa.txt","_Mutect2", files_1[f])
  if (is.null(y) == FALSE) {y$sample = c(sampleID)}
  setwd(output)
  if (is.null(y) == FALSE) {write.table(y, file = files_2[f], sep = "\t", row.names = F); t = rbind(t, y)}
  write.table(table1, file = files_1[f], sep = "\t", row.names = F)

  
  for (i in 1:nrow(BD)) {
    BD$AF[i] = as.numeric(BD$ADvar1[i]) /  as.numeric(BD$ADref[i])}
  
  # Actualizando la base de datos
  
  BD1 = BD[,c(1,2,4,5,34)]
  BDF = BD1
  sampleID = gsub("_somatic.hg38_multianno.txt_completa.txt","_Mutect2", files_1[f])
  setwd(data_base)
  colnames(BDF) = c("Chr", "Start", "Ref", "Alt", "FREQ_AL")
  BDF$FREQ_AL = as.numeric(BDF$FREQ_AL)
  write.table(BDF, toString(sampleID), sep = "\t", row.names = F)
  
  por = f*100/length(files_1)
  cat(round(por, 0) , "% completado\n")
}
setwd(output)
write.table(t, "Variantes.Mutect2", sep = "\t", row.names = F)
setwd(final_output)
write.table(t, "Variantes.Mutect2", sep = "\t", row.names = F)
cat("Done!\n")





# ----------------------------------------------------------------------------------------


# Strelka2


cat("\nFiltering the Strelka2 annotation files...\n")
directory = c("~/.temporal/Strelka2")
directory_output = c("~/.temporal/Strelka2/filter/")
setwd(directory)
files_1 = list.files(directory, pattern = "_snp.hg38_multianno.txt_completa.txt")
files_2 = list.files(directory, pattern = "_snp.hg38_multianno.txt_filtrada.txt")
t = NULL

for (f in 1:length(files_1)) {
  setwd(directory)
  table1 = read.table(files_1[f], sep = "\t", header=T, quote="")
  setwd(directory)
  table2 = read.table(files_2[f], sep = "\t")
  
  tableV11 = table2[11]
  table11 = separate(tableV11, 1, c("GTn", "DPn", "FGTn", "SGTn", "SUBGTn", "AUn", "CUn", "GUn", "TUn"), sep = ":")
  
  table111 = separate(table11, 6, c("AUn1", "AUn2"), sep = ",")
  table11 = separate(table111, 8, c("CUn1", "CUn2"), sep = ",")
  table111 = separate(table11, 10, c("GUn1", "GUn2"), sep = ",")
  table11 = separate(table111, 12, c("TUn1", "TUn2"), sep = ",")
  
  tableV12 = table2[12]
  table12 = separate(tableV12, 1, c("GTt", "DPt", "FGTt", "SGTt", "SUBGTt", "AUt", "CUt", "GUt", "TUt"), sep = ":")
  
  table112 = separate(table12, 6, c("AUt1", "AUt2"), sep = ",")
  table12 = separate(table112, 8, c("CUt1", "CUt2"), sep = ",")
  table112 = separate(table12, 10, c("GUt1", "GUt2"), sep = ",")
  table12 = separate(table112, 12, c("TUt1", "TUt2"), sep = ",")
  
  tableV9 = table2[9]
  table4 = separate(tableV9, 1, c("SOMATIC", "QSS", "TQSS", "NT", "QSS_NT", "TQSS_NT", "SGT", "DP", "MQ", "MQ0", "ReadPosRankSum", "SNVSB", "SomaticEVS"), sep = ";")
  
  table5 = table2[8]
  colnames(table5)=c("Q")
  table1= cbind(table1, table5, table4, table11, table12)
  table1$FA = c(0)
  table1$FAn = c(0)
  
  for (i in 1:nrow(table1)) {
    if (table1$Ref[i] == "A" && table1$Alt[i] == "C" && as.numeric(table1$AUt2[i]) != 0) {
      table1$FA[i] = as.numeric(table1$CUt2[i])/as.numeric(table1$AUt2[i])
      table1$FAn[i] = as.numeric(table1$CUn2[i])/as.numeric(table1$AUn2[i])}
    if (table1$Ref[i] == "A" && table1$Alt[i] == "T" && as.numeric(table1$AUt2[i]) != 0) {
      table1$FA[i] = as.numeric(table1$TUt2[i])/as.numeric(table1$AUt2[i])
      table1$FAn[i] = as.numeric(table1$TUn2[i])/as.numeric(table1$AUn2[i])}
    if (table1$Ref[i] == "A" && table1$Alt[i] == "G" && as.numeric(table1$AUt2[i]) != 0) {
      table1$FA[i] = as.numeric(table1$GUt2[i])/as.numeric(table1$AUt2[i])
      table1$FAn[i] = as.numeric(table1$GUn2[i])/as.numeric(table1$AUn2[i])}
    
    if (table1$Ref[i] == "C" && table1$Alt[i] == "A" && as.numeric(table1$CUt2[i]) != 0) {
      table1$FA[i] = as.numeric(table1$AUt2[i])/as.numeric(table1$CUt2[i])
      table1$FAn[i] = as.numeric(table1$AUn2[i])/as.numeric(table1$CUn2[i])}
    if (table1$Ref[i] == "C" && table1$Alt[i] == "T" && as.numeric(table1$CUt2[i]) != 0) {
      table1$FA[i] = as.numeric(table1$TUt2[i])/as.numeric(table1$CUt2[i])
      table1$FAn[i] = as.numeric(table1$TUn2[i])/as.numeric(table1$CUn2[i])}
    if (table1$Ref[i] == "C" && table1$Alt[i] == "G" && as.numeric(table1$CUt2[i]) != 0) {
      table1$FA[i] = as.numeric(table1$GUt2[i])/as.numeric(table1$CUt2[i])
      table1$FAn[i] = as.numeric(table1$GUn2[i])/as.numeric(table1$CUn2[i])}
    
    if (table1$Ref[i] == "T" && table1$Alt[i] == "A" && as.numeric(table1$TUt2[i]) != 0) {
      table1$FA[i] = as.numeric(table1$AUt2[i])/as.numeric(table1$TUt2[i])
      table1$FAn[i] = as.numeric(table1$AUn2[i])/as.numeric(table1$TUn2[i])}
    if (table1$Ref[i] == "T" && table1$Alt[i] == "C" && as.numeric(table1$TUt2[i]) != 0) {
      table1$FA[i] = as.numeric(table1$CUt2[i])/as.numeric(table1$TUt2[i])
      table1$FAn[i] = as.numeric(table1$CUn2[i])/as.numeric(table1$TUn2[i])}
    if (table1$Ref[i] == "T" && table1$Alt[i] == "G" && as.numeric(table1$TUt2[i]) != 0) {
      table1$FA[i] = as.numeric(table1$GUt2[i])/as.numeric(table1$TUt2[i])
      table1$FAn[i] = as.numeric(table1$GUn2[i])/as.numeric(table1$TUn2[i])}
    
    if (table1$Ref[i] == "G" && table1$Alt[i] == "A" && as.numeric(table1$GUt2[i]) != 0) {
      table1$FA[i] = as.numeric(table1$AUt2[i])/as.numeric(table1$GUt2[i])
      table1$FAn[i] = as.numeric(table1$AUn2[i])/as.numeric(table1$GUn2[i])}
    if (table1$Ref[i] == "G" && table1$Alt[i] == "T" && as.numeric(table1$GUt2[i]) != 0) {
      table1$FA[i] = as.numeric(table1$TUt2[i])/as.numeric(table1$GUt2[i])
      table1$FAn[i] = as.numeric(table1$TUn2[i])/as.numeric(table1$GUn2[i])}
    if (table1$Ref[i] == "G" && table1$Alt[i] == "C" && as.numeric(table1$GUt2[i]) != 0) {
      table1$FA[i] = as.numeric(table1$CUt2[i])/as.numeric(table1$GUt2[i])
      table1$FAn[i] = as.numeric(table1$CUn2[i])/as.numeric(table1$GUn2[i])}
  }
  
  x = c()
  y = c()
  z = c()
  for (i in 1:nrow(table1)) {
    if (table1$Func.refGene[i] == "exonic") {x = rbind(x, table1[i,])}
    else if (table1$Func.refGene[i] == "splicing") {x <- rbind(x, table1[i,])}
    else if (table1$Func.refGene[i] == "UTR5") {x <- rbind(x, table1[i,])}
    else if (table1$Func.refGene[i] == "UTR3") {x <- rbind(x, table1[i,])}
  }
  
  for (i in 1:nrow(table1)) {
    if (toString(x$Gene.refGene[i]) == "MUC1") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "MUC2") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "MUC3") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "MUC4") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "POTEJ") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "POTEM") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "POTEE") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "ZNF91") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "NOTCH2NLB") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "NOTCH2NLC") {z = c(z, i)}
    else if (toString(x$Gene.refGene[i]) == "FRMPD3") {z = c(z, i)}}
  if (toString(is.null(z)) == "FALSE") {x = x[-c(z),]}
  
  genes = unique(as.vector(bedfile$V4))
  n = c()
  
  for (i in 1:nrow(x)) {
    gene = toString(x$Gene.refGene[i])
    p = which(genes %in% gene)
    if (length(p) == 1) {n = c(n, i)}
  }
  x = x[c(n),]
  
  z = NULL
  if (is.null(x) == FALSE) {for (i in 1:nrow(x)) {
    if (as.numeric(x$DPn[i]) > 29 && as.numeric(x$DPt[i]) > 499) {y = rbind(y, x[i,])}}}
  
  BD = y
  
  x = NULL
  if (is.null(y) == FALSE) {for (i in 1:nrow(y)) {
    if (y$FA[i] > 0.005) {x = rbind(x, y[i,])}}}
  
  sampleID = gsub("_snp.hg38_multianno.txt_completa.txt","_Strelka", files_1[f])
  if (is.null(x) == FALSE) {x$SampleID = c(sampleID)}
  setwd(directory_output)
  if (is.null(x) == FALSE) {write.table(x, file = files_2[f], sep = "\t", row.names = F); t = rbind(t, x)}
  write.table(table1, file = files_1[f], sep = "\t", row.names = F)
  
  # Actualizando la base de datos
  
  sampleID = gsub("_snp.hg38_multianno.txt_completa.txt","_Strelka", files_1[f])
  BDF = BD[,c(1,2,4,5,(ncol(BD)-1))]
  colnames(BDF) = c("Chr", "Start", "Ref", "Alt", "FREQ_AL")
  BDF$FREQ_AL = as.numeric(BDF$FREQ_AL)
  setwd(data_base)
  write.table(BDF, toString(sampleID), sep = "\t", row.names = F)
  por = f*100/length(files_1)
  cat(round(por, 0) , "% completado\n")
}
setwd(directory_output)
write.table(t, "Variantes.SNVs.Strelka", sep = "\t", row.names = F)

setwd(final_output)
write.table(t, "Variantes.SNVs.Strelka", sep = "\t", row.names = F)
cat("Done!\n")





# ----------------------------------------------------------------------------------------



# VarScan2

# Filtrado de variantes de cfDNA de tumor y normal

cat("\nFiltering the VarScan2 annotation files...\n")
directorio = c("~/.temporal/VarScan2/")
directorio_output = c("~/.temporal/VarScan2/filter/")
setwd(directorio)
files_1 = list.files(directorio, pattern = ".snp")
files_2 = list.files(directorio, pattern = ".indel")

anotacion = c("anotacion")

PanelSize =  247433
Pvalue = 5/(PanelSize*4)

for (f in 1:length(files_1)) {
  
  # SNPs
  
  setwd(directorio)
  TD = read.table(files_1[f], sep = "\t", header = T, quote= "")
  TD$difVF = c(NA)
  TD$normal_var_freq = as.numeric(gsub("%", "", TD$normal_var_freq))
  TD$tumor_var_freq = as.numeric(gsub("%", "", TD$tumor_var_freq))
  x = c()
  y = c()
  t = c()
  r = c()
  b = c()
  for (l in 1:nrow(TD)) {
    TD$difVF[l] = TD$normal_var_freq[l] - TD$tumor_var_freq[l]
    TD$normalDP[l] = as.numeric(TD$normal_reads1[l] + TD$normal_reads2[l])
    TD$tumorDP[l] = as.numeric(TD$tumor_reads1[l] + TD$tumor_reads2[l])
    if (as.numeric(TD$tumor_var_freq[l]) > as.numeric(0.5)) {r = c(r, l)}
    if (TD$var[l] != "") {y = c(y, l)}
    if (TD$somatic_status[l] == "Somatic") {t = c(t, l)}
    if (TD$normalDP[l] > 29) {x = c(x, l)}
    if (TD$tumorDP[l] > 449) {b = c(b, l)}}
  
  g = intersect(y, x)
  z = intersect(g, t)
  g = intersect(z, b)
  validationSNPs = TD[g, ]
  k = intersect(g, r)
  TDfilter = TD[k, ]
  TDfilter$gene = c("")
  for (i in 1:nrow(TDfilter)) {
    for (b in 1:nrow(bedfile)) {
      cal1 = TDfilter$position[i] - bedfile$V2[b]
      cal2 = bedfile$V3[b] - TDfilter$position[i]
      if (TDfilter$chrom[i] == bedfile$V1[b] && cal1 > 0 && cal2 > 0) {
        TDfilter$gene[i] = toString(bedfile$V4[b])
      }
    }}
  
  setwd(directorio_output)
  write.table(TDfilter, toString(files_1[f]), row.names = F, sep = "\t")
  
  # Indels
  
  setwd(directorio)
  TD = read.table(files_2[f], sep = "\t", header = T, quote= "")
  TD$difVF = c(NA)
  TD$normal_var_freq = as.numeric(gsub("%", "", TD$normal_var_freq))
  TD$tumor_var_freq = as.numeric(gsub("%", "", TD$tumor_var_freq))
  x = c()
  y = c()
  t = c()
  r = c()
  b = c()
  for (l in 1:nrow(TD)) {
    TD$difVF[l] = TD$normal_var_freq[l] - TD$tumor_var_freq[l]
    TD$normalDP[l] = as.numeric(TD$normal_reads1[l] + TD$normal_reads2[l])
    TD$tumorDP[l] = as.numeric(TD$tumor_reads1[l] + TD$tumor_reads2[l])
    if (as.numeric(TD$tumor_var_freq[l]) > as.numeric(0.5)) {r = c(r, l)}
    if (TD$var[l] != "") {y = c(y, l)}
    if (TD$somatic_status[l] == "Somatic") {t = c(t, l)}
    if (TD$normalDP[l] > 29) {x = c(x, l)}
    if (TD$tumorDP[l] > 449) {b = c(b, l)}}
  
  g = intersect(y, x)
  z = intersect(g, t)
  g = intersect(z, b)
  validationIndel = TD[g, ]
  k = intersect(g, r)
  TDfilterIndel = TD[k, ]
  TDfilterIndel$gene = c("")
  for (i in 1:nrow(TDfilterIndel)) {
    for (b in 1:nrow(bedfile)) {
      cal1 = TDfilterIndel$position[i] - bedfile$V2[b]
      cal2 = bedfile$V3[b] - TDfilterIndel$position[i]
      if (TDfilterIndel$chrom[i] == bedfile$V1[b] && cal1 > 0 && cal2 > 0) {
        TDfilterIndel$gene[i] = toString(bedfile$V4[b])
      }
    }}
  
  setwd(directorio_output)
  write.table(TDfilterIndel, toString(files_2[f]), row.names = F, sep = "\t")
  
  # Actualizando base de datos
  
  validation = rbind(validationSNPs, validationIndel)
  validation = validation[,c(1,2,3,4,11)]
  colnames(validation) = c("Chr", "Start", "Ref", "Alt", "FREQ_AL")
  validation$FREQ_AL = as.numeric(validation$FREQ_AL)
  
  sampleID = gsub(".snp","_VarScan2", files_1[f])
  setwd(data_base)
  write.table(validation, toString(sampleID), sep = "\t", row.names = F)
  
  # Uniendo archivos en uno solo
  
  sampleID = gsub(".snp","_VarScan2", files_1[f])
  if (is.data.frame(library) == T) {
    TDfilter$sampleID = c(toString(sampleID))
    TDfilterIndel$sampleID = c(toString(sampleID))
    library = rbind(library, TDfilter, TDfilterIndel)}
  
  if (is.data.frame(library) == F) {
    TDfilter$sampleID = c(toString(sampleID))
    TDfilterIndel$sampleID = c(toString(sampleID))
    library = rbind(TDfilter, TDfilterIndel)}
  
  setwd(directorio_output)
  write.table(library, "Variantes.VarScan2", row.names = F, sep = "\t")
  setwd(final_output)
  write.table(library, "Variantes.VarScan2", row.names = F, sep = "\t")

  por = f*100/length(files_1)
  cat(round(por, 0) , "% completado")
}

cat("Done!\n")

