#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library('plotrix')
library("funr")

# library(phangorn)
# library(vcfR)
# library(phytools)
# library(dplyr)
# 
# library(foreach)
# library(doParallel)
current_folder = dirname(sys.script())

print(paste0("Number of arguments: ", length(args))  )
if (length(args) != 5) {
  print("gene.test.xxx.csv: pvalue of all gene")
  print("output.pdf: output name")
  stop(" Usage: hit.plot.R <gene.test.all.csv> <gene.test.sens_inh.csv> <gene.test.sens_mdr.csv> <gene.test.sens.csv> <output.pdf>", call.=FALSE)
}

gene_mapper = read.csv(paste0(current_folder, '/Candidate/candidate_genes.csv'))
gene_mapper = gene_mapper[which(gene_mapper$locus_tag != 'not_applicable'),]
candiate_gene = as.character(gene_mapper$locus_tag)

GetDotSize <- function(stepNum){
  if(stepNum <= 10){
    return(2)
  }else{
    if(stepNum <= 20){
      return(4)
    }else{
      if(stepNum <= 30){
        return(6)
      }else{
        return(8)
      }
    }
    
  }
}

map <<- read.csv(paste0(current_folder, '/H37RV/H37rv.intervals.genbank.csv'),header = F, sep = "\t") 
GetGeneInfo <- function(locustag){
  idx = which(map$V2==locustag)
  return(c(as.character(map$V2[idx]), as.character(map$V3[idx]), as.numeric(map$V4[idx]), as.numeric(map$V5[idx])))
} 

setEPS()
postscript(args[5], width = 40, height = 20)
# pdf(file = args[8], width = 40, height = 20)

par(mfrow=c(2,3))
for (kk in 1:4) {
  plot(-1,-1,xlim = c(0,4.5), ylim = c(-0.5, 13),xlab ='Genomic position (Mbp)',ylab='-log(P-value)')
  
  data = read.csv(args[kk])
  Threshold = 0.9/200000
  data$Pvalue[which(data$Pvalue == 0)] = Threshold
  data$Pvalue = -log(data$Pvalue)
  
  data2 = data[which(data$Gene %in% candiate_gene), ]
  data1 = data[which(!data$Gene %in% candiate_gene), ]
  
  
  
  for(i in 1:dim(data1)[1]){
    locustag = as.character(data1$Gene[i])
    sitenum = data1$SiteNum[i]
    stepnum = data1$StepNum[i]
    info = GetGeneInfo(locustag)
    r = GetDotSize(stepnum)
    mycolor = c('gray', 'white')
    
    x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
    y = data1$Pvalue[i]
    
    
    border = as.character(data1$border[i])
    points(x,y,pch=21,bg=mycolor[1], col = mycolor[2], cex=r)
  }
  
  
  
  for(i in 1:dim(data2)[1]){
    locustag = as.character(data2$Gene[i])
    sitenum = data2$SiteNum[i]
    stepnum = data2$StepNum[i]
    info = GetGeneInfo(locustag)
    r = GetDotSize(stepnum)
    mycolor = c('red', 'black')
    
    x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
    y = data2$Pvalue[i]
    
    
    border = as.character(data2$border[i])
    points(x,y,pch=21,bg=mycolor[1], col = mycolor[2], cex=r)
  }
  
  
  
  for(i in 1:dim(data2)[1]){
    locustag = as.character(data2$Gene[i])
    info = GetGeneInfo(locustag)
    t = info[2]
    x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
    y = data2$Pvalue[i]
    text(x,y,t)
  }
  
  points(4,10,pch=21, cex=8)
  points(4,10,pch=21, cex=6)
  points(4,10,pch=21, cex=4)
  points(4,10,pch=21, cex=2)
}


# 
# data = read.csv(args[1])
# 
# Threshold = 0.9/200000
# data$Pvalue[which(data$Pvalue == 0)] = Threshold
# data$Pvalue = -log(data$Pvalue)
# 
# plot(-1,-1,xlim = c(0,4.5), ylim = c(-0.5, 13),xlab ='Genomic position (Mbp)',ylab='-log(P-value)')
# 
# 
# 
# for(i in 1:dim(data)[1]){
#   locustag = as.character(data$Gene[i])
#   sitenum = data$SiteNum[i]
#   stepnum = data$StepNum[i]
#   info = GetGeneInfo(locustag)
#   r = GetDotSize(stepnum)
#   if(i <= 10){
#     mycolor = c('red', 'black')
#   }else{
#     mycolor = c('gray', 'white')
#   }
#   
#   x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
#   y = data$Pvalue[i]
#   
#   
#   border = as.character(data$border[i])
#   points(x,y,pch=21,bg=mycolor[1], col = mycolor[2], cex=r)
#   
# }
# 
# for(i in 1:10){
#   
#   locustag = as.character(data$Gene[i])
#   info = GetGeneInfo(locustag)
#   t = info[2]
#   x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
#   y = data$Pvalue[i]
#   
#   text(x,y,t)
#   
# }
# 
# points(4,10,pch=21, cex=4)
# points(4,10,pch=21, cex=3)
# points(4,10,pch=21, cex=2)
# points(4,10,pch=21, cex=1)
# 
# 
# 
# 
# 
# data = read.csv(args[2])
# 
# Threshold = 0.9/200000
# data$Pvalue[which(data$Pvalue == 0)] = Threshold
# data$Pvalue = -log(data$Pvalue)
# 
# plot(-1,-1,xlim = c(0,4.5), ylim = c(-0.5, 13),xlab ='Genomic position (Mbp)',ylab='-log(P-value)')
# 
# 
# 
# for(i in 1:dim(data)[1]){
#   locustag = as.character(data$Gene[i])
#   sitenum = data$SiteNum[i]
#   stepnum = data$StepNum[i]
#   info = GetGeneInfo(locustag)
#   r = GetDotSize(stepnum)
#   if(i <= 10){
#     mycolor = c('red', 'black')
#   }else{
#     mycolor = c('gray', 'white')
#   }
#   
#   x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
#   y = data$Pvalue[i]
#   
#   
#   border = as.character(data$border[i])
#   points(x,y,pch=21,bg=mycolor[1], col = mycolor[2], cex=r)
#   
# }
# 
# for(i in 1:10){
#   
#   locustag = as.character(data$Gene[i])
#   info = GetGeneInfo(locustag)
#   t = info[2]
#   x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
#   y = data$Pvalue[i]
#   
#   text(x,y,t)
#   
# }
# 
# 
# 
# 
# 
# 
# data = read.csv(args[3])
# 
# Threshold = 0.9/200000
# data$Pvalue[which(data$Pvalue == 0)] = Threshold
# data$Pvalue = -log(data$Pvalue)
# 
# plot(-1,-1,xlim = c(0,4.5), ylim = c(-0.5, 13),xlab ='Genomic position (Mbp)',ylab='-log(P-value)')
# 
# 
# 
# for(i in 1:dim(data)[1]){
#   locustag = as.character(data$Gene[i])
#   sitenum = data$SiteNum[i]
#   stepnum = data$StepNum[i]
#   info = GetGeneInfo(locustag)
#   r = GetDotSize(stepnum)
#   if(i <= 10){
#     mycolor = c('red', 'black')
#   }else{
#     mycolor = c('gray', 'white')
#   }
#   
#   x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
#   y = data$Pvalue[i]
#   
#   
#   border = as.character(data$border[i])
#   points(x,y,pch=21,bg=mycolor[1], col = mycolor[2], cex=r)
#   
# }
# 
# for(i in 1:10){
#   
#   locustag = as.character(data$Gene[i])
#   info = GetGeneInfo(locustag)
#   t = info[2]
#   x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
#   y = data$Pvalue[i]
#   
#   text(x,y,t)
#   
# }
# 
# 
# 
# data = read.csv('gene.test.1139.csv')
# 
# Threshold = 0.9/200000
# data$Pvalue[which(data$Pvalue == 0)] = Threshold
# data$Pvalue = -log(data$Pvalue)
# 
# plot(-1,-1,xlim = c(0,4.5), ylim = c(-0.5, 13),xlab ='Genomic position (Mbp)',ylab='-log(P-value)')
# 
# 
# 
# for(i in 1:dim(data)[1]){
#   locustag = as.character(data$Gene[i])
#   sitenum = data$SiteNum[i]
#   stepnum = data$StepNum[i]
#   info = GetGeneInfo(locustag)
#   r = GetDotSize(stepnum)
#   if(i <= 10){
#     mycolor = c('red', 'black')
#   }else{
#     mycolor = c('gray', 'white')
#   }
#   
#   x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
#   y = data$Pvalue[i]
#   
#   
#   border = as.character(data$border[i])
#   points(x,y,pch=21,bg=mycolor[1], col = mycolor[2], cex=r)
#   
# }
# 
# for(i in 1:10){
#   
#   locustag = as.character(data$Gene[i])
#   info = GetGeneInfo(locustag)
#   t = info[2]
#   x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
#   y = data$Pvalue[i]
#   
#   text(x,y,t)
#   
# }
# 
# 
# 
# data = read.csv('gene.test.351.csv')
# 
# Threshold = 0.9/200000
# data$Pvalue[which(data$Pvalue == 0)] = Threshold
# data$Pvalue = -log(data$Pvalue)
# 
# plot(-1,-1,xlim = c(0,4.5), ylim = c(-0.5, 13),xlab ='Genomic position (Mbp)',ylab='-log(P-value)')
# 
# 
# 
# for(i in 1:dim(data)[1]){
#   locustag = as.character(data$Gene[i])
#   sitenum = data$SiteNum[i]
#   stepnum = data$StepNum[i]
#   info = GetGeneInfo(locustag)
#   r = GetDotSize(stepnum)
#   if(i <= 10){
#     mycolor = c('red', 'black')
#   }else{
#     mycolor = c('gray', 'white')
#   }
#   
#   x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
#   y = data$Pvalue[i]
#   
#   
#   border = as.character(data$border[i])
#   points(x,y,pch=21,bg=mycolor[1], col = mycolor[2], cex=r)
#   
# }
# 
# for(i in 1:10){
#   
#   locustag = as.character(data$Gene[i])
#   info = GetGeneInfo(locustag)
#   t = info[2]
#   x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
#   y = data$Pvalue[i]
#   
#   text(x,y,t)
#   
# }
# 
# 
# data = read.csv('gene.test.292.csv')
# 
# Threshold = 0.9/200000
# data$Pvalue[which(data$Pvalue == 0)] = Threshold
# data$Pvalue = -log(data$Pvalue)
# 
# plot(-1,-1,xlim = c(0,4.5), ylim = c(-0.5, 13),xlab ='Genomic position (Mbp)',ylab='-log(P-value)')
# 
# 
# 
# for(i in 1:dim(data)[1]){
#   locustag = as.character(data$Gene[i])
#   sitenum = data$SiteNum[i]
#   stepnum = data$StepNum[i]
#   info = GetGeneInfo(locustag)
#   r = GetDotSize(stepnum)
#   if(i <= 10){
#     mycolor = c('red', 'black')
#   }else{
#     mycolor = c('gray', 'white')
#   }
#   
#   x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
#   y = data$Pvalue[i]
#   
#   
#   border = as.character(data$border[i])
#   points(x,y,pch=21,bg=mycolor[1], col = mycolor[2], cex=r)
#   
# }
# 
# for(i in 1:10){
#   
#   locustag = as.character(data$Gene[i])
#   info = GetGeneInfo(locustag)
#   t = info[2]
#   x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
#   y = data$Pvalue[i]
#   
#   text(x,y,t)
#   
# }
# 
# 
# 
# 
# data = read.csv('gene.test.59.csv')
# 
# Threshold = 0.9/200000
# data$Pvalue[which(data$Pvalue == 0)] = Threshold
# data$Pvalue = -log(data$Pvalue)
# 
# plot(-1,-1,xlim = c(0,4.5), ylim = c(-0.5, 13),xlab ='Genomic position (Mbp)',ylab='-log(P-value)')
# 
# 
# 
# for(i in 1:dim(data)[1]){
#   locustag = as.character(data$Gene[i])
#   sitenum = data$SiteNum[i]
#   stepnum = data$StepNum[i]
#   info = GetGeneInfo(locustag)
#   r = GetDotSize(stepnum)
#   if(i <= 10){
#     mycolor = c('red', 'black')
#   }else{
#     mycolor = c('gray', 'white')
#   }
#   
#   x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
#   y = data$Pvalue[i]
#   
#   
#   border = as.character(data$border[i])
#   points(x,y,pch=21,bg=mycolor[1], col = mycolor[2], cex=r)
#   
# }
# 
# for(i in 1:10){
#   
#   locustag = as.character(data$Gene[i])
#   info = GetGeneInfo(locustag)
#   t = info[2]
#   x = (as.numeric(info[3]) + as.numeric(info[4]))/2000000
#   y = data$Pvalue[i]
#   
#   text(x,y,t)
#   
# }
# 
# 
