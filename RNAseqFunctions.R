require(DESeq2)
require(plyr)
require(ggplot2)

class.perc <- function(df, colname){
  
  df.count <- count(df, colname)
  df.count <- data.frame(df.count,
                         'percentage' = df.count$freq / sum(df.count$freq))
  return(df.count)
}

class.enrichment <- function(exp.rel.count, bg.rel.count){
  comp.df <- NULL
  for(class in exp.rel.count$Class){
    dat.e <- subset(exp.rel.count, Class == class)
    dat.bg <- subset(bg.rel.count, Class == class)
    
    enr <- dat.e$percentage / dat.bg$percentage
    enr.vec <- data.frame('Class' = class, 'Enrichment' = enr)
    comp.df <- rbind(comp.df, enr.vec)
  }
  return(comp.df)
}

run.deseq2 <- function(count.matrix, colDescriptionVector, contrast){
  colIDs <- data.frame('sample' = colnames(count.matrix),
                       'condition' = colDescriptionVector)
  dds <- DESeqDataSetFromMatrix(count.matrix, colData = colIDs, design=~condition)
  dds <- DESeq(dds)
  rld <- rlog(dds)
  plotPCA(rld)
  plotMA(dds)
  
  contrast <- c('condition', contrast)
  res <- results(dds, contrast = contrast)
  
  res.df <- as.data.frame(res)
  res.df <- data.frame('GeneID' = row.names(res.df), res.df, row.names = NULL)
  return(res.df)
}

slideFunct <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- median(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}