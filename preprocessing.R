preprocessing <- function(dir, dataset, words_percent, docs_percent){

  stopwords1 <- read.table(paste(dir,'/Recover/code/stopwords.txt',sep=''))
  stopwords1 <- as.matrix(stopwords1)
  
  stopwords2 <- read.table(paste(dir,'/TSVD/code/stopwords.txt',sep=''))
  stopwords2 <- as.matrix(stopwords2)
  
  stopwords <- unique(c(stopwords1, stopwords2))
  
  setwd(paste(dir,'/real_data',sep=''))
  
  vocab <- read.table(paste(dataset, '.vocab.txt', sep=''))
  vocab <- as.matrix(vocab)
  data <- read.table(paste(dir,'/real_data/',dataset,'.txt',sep=''),header=FALSE,sep=" ")
  data <- as.matrix(data)
  
  # eliminate the stop words from the dataset.txt and dataset.vacab.txt
  del_stopwords <- rep(0,length(vocab))
  t <- 1
  for (i in 1:length(vocab)){
    if (!vocab[i]%in%stopwords){
      del_stopwords[i] <- t
      t <- t+1
    }
  }
  
  vocab <- vocab[del_stopwords!=0]
  data <- data[del_stopwords[data[,2]]!=0,]
  data[,2] <- del_stopwords[data[,2]]
  
  n <- max(data[,1])
  p <- max(data[,2])
  D <- matrix(0,p,n)
  for (t in 1:dim(data)[1]){
    D[data[t,2], data[t,1]] <-data[t,3]
  }
  
  # keep top words_percent most frequent words
  D_rowsum <- rowSums(D)
  del_low_words <- rep(0,dim(D)[1])
  threshold <- quantile(D_rowsum, 1-words_percent)
  t <- 1
  for (i in 1:dim(D)[1]){
    if (D_rowsum[i] >= threshold){
      del_low_words[i] <- t
      t <- t+1
    }
  }
  
  vocab <- vocab[del_low_words!=0]
  D <- D[del_low_words!=0,] 
  
  # keep top docs_percent largest documents
  D_colsum <- colSums(D)
  threshold <- quantile(D_colsum, 1-docs_percent)
  D <- D[,D_colsum >= threshold]
  
  # write the data after preprocessing into the files
  fileConn<-file(paste("prepro_", words_percent,"_",docs_percent,'_', dataset, ".txt", sep=""),'w')
  close(fileConn)
  fileConn<-file(paste("prepro_", words_percent,"_",docs_percent,'_', dataset, ".txt", sep=""),'a')
  for (i in 1:dim(D)[1]){
    for (j in 1:dim(D)[2]){
      if (D[i,j]>0){
        writeLines(paste(j,i,D[i,j]), fileConn)
      }
    }
  }
  close(fileConn)
  
  fileConn<-file(paste("prepro_", words_percent,"_",docs_percent,'_', dataset, '.vocab', ".txt", sep=""),'w')
  close(fileConn)
  fileConn<-file(paste("prepro_", words_percent,"_",docs_percent,'_', dataset, '.vocab', ".txt", sep=""),'a')
  for (i in 1:length(vocab)){
    writeLines(vocab[i], fileConn)
  }
  close(fileConn)
}


