semi_synthetic_generation <- function(dir, dataset, K, N=NULL, seed = NULL){
  # p is the number of words in the dictionary.
  # n is the number of documents.
  # N is a vector of length n, with ith entry as the total number of words in the ith documen.
  
  library(tm)
  library(topicmodels)
  library(R.matlab)
  setwd(dir)
  
  # read from sparse matrix txt file 
  data = read.table(paste(dir,'/real_data','/',dataset,'.txt', sep=''),header=FALSE,sep="");
  data <- as.matrix(data)
  
  # specify the output file names
  if (is.null(N)){
    out_txt <- paste("semi_synthetic_", K,"_", dataset, ".txt", sep="")
    out_A <- paste("semi_synthetic_", K,"_", dataset, ".A", sep="")
  } else{
    out_txt <- paste("semi_synthetic_", K,"_", N,"_", dataset, ".txt", sep="")
    out_A <- paste("semi_synthetic_", K,"_", N,"_", dataset, ".A", sep="")
  }

  n <- max(data[,1])
  p <- max(data[,2])
  D <- matrix(0, n, p)
  for (i in 1:dim(data)[1]){
    D[data[i,1], data[i,2]] <-data[i,3]
  }
  
  if (is.null(N)){
    N <- rowSums(D)
  } else if(length(N)==1){
    N <- rep(N,n)
  }
  
  # use LDA to fit the model
  if (is.null(seed)){
    LDA_obj <- LDA(D, K, method = 'VEM')
  } else{
    LDA_obj <- LDA(D, K, method = 'VEM', control=list(seed=seed))
  }
  set.seed(seed)
  D_sim <- LDA_obj@gamma%*%exp(LDA_obj@beta)
  D_synth <- matrix(0, n, p)
  for (i in 1:length(N)){
    D_synth[i,] <- rmultinom(1, N[i], D_sim[i,])
  }
  
  #write the sparse semi-synthetic data matrix into txt file
  setwd(paste(dir,'/semi_synthetic_data', sep=""));
  fileConn<-file(out_txt,'w')
  close(fileConn)
  fileConn<-file(out_txt,'a')
  for (i in 1:n){
    for (j in 1:p){
      if (D_synth[i,j]>0){
        writeLines(paste(i,j,D_synth[i,j]), fileConn)
      }
    }
  }
  close(fileConn)
  
  writeMat(out_A, A=exp(t(LDA_obj@beta)))
}
