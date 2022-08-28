# synthetic data generation

# m_n is the number of pure documents for each topic
# m_p is the number of anchor words for each topic
# delta_p is the magnitude of the probability mass of each anchor words in each topic

# synthetic_generation_1: Anchor words + homogeneity
synthetic_generation_1 <- function(dir, dataset, p, n, N, K, m_p, delta_p, m_n, W_opt = "puredoc", seed = NULL){
  set.seed(seed)
  #Generate A0
  A0 <- matrix(0,p,K)
  for (k in 1:K){
    A0[((k-1)*m_p+1):(k*m_p),k] <- rep(delta_p,m_p)
  }
  
  A0[(K*m_p+1):p,] <- (1-m_p*delta_p)*
    apply(matrix(runif((p-K*m_p)*K), p-K*m_p, K), 2, function(x) x/sum(x))
  
  #Generate W0
  if (W_opt == "puredoc"){
    W0 <- matrix(0,K,n)
    if (m_n>0){
      for (k in 1:K){
        W0[k,((k-1)*m_n+1):(k*m_n)] <- rep(1,m_n)
      }
    }
    W0[,(K*m_n+1):n] <- apply(matrix(runif((n-K*m_n)*K), K, n-K*m_n), 2, function(x) x/sum(x))
  }else if (W_opt == "nopuredoc"){
    W0 <- apply(matrix(runif(n*K), K, n), 2, function(x) x/sum(x))
  }else{
    stop('Incorrect imput for W_opt')
  }
      
  #Generate D0
  D0 <- A0%*%W0

  if (length(N)!=n){
    N <- rep(N, n)
  }
  #Generate X
  X <- matrix(0,p,n)
  for (j in 1:n){
    X[,j] <- as.numeric(rmultinom(1, N[j], D0[,j]))
  }
  
  #Compute D
  D <- t(t(X)/N)
  
  #write the sparse semi-synthetic data matrix into txt file
  setwd(paste(dir,'/synthetic_data', sep=""));
  fileConn<-file(paste(dataset, ".txt", sep=""),'w')
  close(fileConn)
  fileConn<-file(paste(dataset, ".txt", sep=""),'a')
  for (i in 1:p){
    for (j in 1:n){
      if (X[i,j]>0){
        writeLines(paste(j,i,X[i,j]), fileConn)
      }
    }
  }
  close(fileConn)
  
  writeMat(paste(dataset, ".A", sep=""), A=A0)
  
  #Generate synthetic vocabulary file, which is required for the Recover method.
  fileConn<-file(paste(dataset, '.vocab.txt', sep=""),'w')
  close(fileConn)
  fileConn<-file(paste(dataset, '.vocab.txt', sep=""),'a')
  for (j in 1:p){
    writeLines(toString(j), fileConn)
  }
  close(fileConn)
}

# synthetic_generation_2: Anchor words + heterogeneity(power law)
synthetic_generation_2 <- function(dir, dataset, p, n, N, K, m_p, delta_p, m_n, stopnum, W_opt = "puredoc"){
  set.seed(seed)
  #Generate A0
  A0 <- matrix(0,p,K)
  for (k in 1:K){
    A0[((k-1)*m_p+1):(k*m_p),k] <- rep(delta_p,m_p)
  }

  for (i in (K*m_p+1):p){
    A0[i,] <- rexp(K,(i+stopnum)^1.07)
  }

  A0[(K*m_p+1):p,] <- (1-m_p*delta_p)*apply(A0[(K*m_p+1):p,], 2, function(x) x/sum(x))

  #Generate W0
  if (W_opt == "puredoc"){
    W0 <- matrix(0,K,n)
    if (m_n>0){
      for (k in 1:K){
        W0[k,((k-1)*m_n+1):(k*m_n)] <- rep(1,m_n)
      }
    }
    W0[,(K*m_n+1):n] <- apply(matrix(runif((n-K*m_n)*K), K, n-K*m_n), 2, function(x) x/sum(x))
  }else if (W_opt == "nopuredoc"){
    W0 <- apply(matrix(runif(n*K), K, n), 2, function(x) x/sum(x))
  }else{
    stop('Incorrect imput for W_opt')
  }
  
  #Generate D0
  D0 <- A0%*%W0
  
  if (length(N)!=n){
    N <- rep(N, n)
  }
  #Generate X
  X <- matrix(0,p,n)
  for (j in 1:n){
    X[,j] <- as.numeric(rmultinom(1, N[j], D0[,j]))
  }
  
  #Compute D
  D <- t(t(X)/N)
  
  #write the sparse semi-synthetic data matrix into txt file
  setwd(paste(dir,'/synthetic_data', sep=""));
  fileConn<-file(paste(dataset, ".txt", sep=""),'w')
  close(fileConn)
  fileConn<-file(paste(dataset, ".txt", sep=""),'a')
  for (i in 1:p){
    for (j in 1:n){
      if (X[i,j]>0){
        writeLines(paste(j,i,X[i,j]), fileConn)
      }
    }
  }
  close(fileConn)
  
  writeMat(paste(dataset, ".A", sep=""), A=A0)
  
  #Generate synthetic vocabulary file, which is required for the Recover method.
  fileConn<-file(paste(dataset, '.vocab.txt', sep=""),'w')
  close(fileConn)
  fileConn<-file(paste(dataset, '.vocab.txt', sep=""),'a')
  for (j in 1:p){
    writeLines(toString(j), fileConn)
  }
  close(fileConn)
}

# synthetic_generation_3: Anchor words + heterogeneity(Unif)
synthetic_generation_3 <- function(dir, dataset, p, n, N, K, m_p, delta_p, m_n, hmax, W_opt = "puredoc"){
  set.seed(seed)
  #Generate A0
  A0 <- matrix(0,p,K)
  for (k in 1:K){
    A0[((k-1)*m_p+1):(k*m_p),k] <- rep(delta_p,m_p)
  }
  
  nmax <- floor((1-m_p*delta_p)/(2*hmax))
  nmin <- p-K*m_p-nmax
  hmin <- (1-m_p*delta_p-hmax*nmax)/(nmin)
  
  #rest <- c(rep(hmax, nmax), rep(hmin, nmin))  
  #for (k in 1:K){
  #  A0[(K*m_p+1):p,k] <- rest[sample(p-K*m_p,p-K*m_p)]
  #}
  A0[(K*m_p+1):(K*m_p+nmax),] <- hmax*matrix(runif(K*nmax), nmax, K)
  A0[(K*m_p+nmax+1):p,] <- hmin*matrix(runif(K*nmin), nmin, K)
  A0[(K*m_p+1):p,] <- (1-m_p*delta_p)*apply(A0[(K*m_p+1):p,], 2, function(x) x/sum(x))
  
  #Generate W0
  if (W_opt == "puredoc"){
    W0 <- matrix(0,K,n)
    if (m_n>0){
      for (k in 1:K){
        W0[k,((k-1)*m_n+1):(k*m_n)] <- rep(1,m_n)
      }
    }
    W0[,(K*m_n+1):n] <- apply(matrix(runif((n-K*m_n)*K), K, n-K*m_n), 2, function(x) x/sum(x))
  }else if (W_opt == "nopuredoc"){
    W0 <- apply(matrix(runif(n*K), K, n), 2, function(x) x/sum(x))
  }else{
    stop('Incorrect imput for W_opt')
  }
  
  #Generate D0
  D0 <- A0%*%W0
  
  if (length(N)!=n){
    N <- rep(N, n)
  }
  #Generate X
  X <- matrix(0,p,n)
  for (j in 1:n){
    X[,j] <- as.numeric(rmultinom(1, N[j], D0[,j]))
  }
  
  #Compute D
  D <- t(t(X)/N)
  
  #write the sparse semi-synthetic data matrix into txt file
  setwd(paste(dir,'/synthetic_data', sep=""));
  fileConn<-file(paste(dataset, ".txt", sep=""),'w')
  close(fileConn)
  fileConn<-file(paste(dataset, ".txt", sep=""),'a')
  for (i in 1:p){
    for (j in 1:n){
      if (X[i,j]>0){
        writeLines(paste(j,i,X[i,j]), fileConn)
      }
    }
  }
  close(fileConn)
  
  writeMat(paste(dataset, ".A", sep=""), A=A0)
  
  #Generate synthetic vocabulary file, which is required for the Recover method.
  fileConn<-file(paste(dataset, '.vocab.txt', sep=""),'w')
  close(fileConn)
  fileConn<-file(paste(dataset, '.vocab.txt', sep=""),'a')
  for (j in 1:p){
    writeLines(toString(j), fileConn)
  }
  close(fileConn)
}

# synthetic_generation_4: No anchor words + homogeneity
synthetic_generation_4 <- function(dir, dataset, p, n, N, K, m_p, delta_p, m_n, decay, W_opt = "puredoc"){
  set.seed(seed)
  #Generate A0
  A0 <- matrix(0,p,K)
  for (k in 1:K){
    A0[((k-1)*m_p+1):(k*m_p),k] <- delta_p
    A0[((k-1)*m_p+1):(k*m_p),-k] <- delta_p*decay
  }
  
  A0[(K*m_p+1):p,] <- (1-m_p*delta_p*(1+(K-1)*decay))*
    apply(matrix(runif((p-K*m_p)*K), p-K*m_p, K), 2, function(x) x/sum(x))
  
  #Generate W0
  if (W_opt == "puredoc"){
    W0 <- matrix(0,K,n)
    if (m_n>0){
      for (k in 1:K){
        W0[k,((k-1)*m_n+1):(k*m_n)] <- rep(1,m_n)
      }
    }
    W0[,(K*m_n+1):n] <- apply(matrix(runif((n-K*m_n)*K), K, n-K*m_n), 2, function(x) x/sum(x))
  }else if (W_opt == "nopuredoc"){
    W0 <- apply(matrix(runif(n*K), K, n), 2, function(x) x/sum(x))
  }else{
    stop('Incorrect imput for W_opt')
  }
  
  #Generate D0
  D0 <- A0%*%W0
  
  if (length(N)!=n){
    N <- rep(N, n)
  }
  #Generate X
  X <- matrix(0,p,n)
  for (j in 1:n){
    X[,j] <- as.numeric(rmultinom(1, N[j], D0[,j]))
  }
  
  #Compute D
  D <- t(t(X)/N)
  
  #write the sparse semi-synthetic data matrix into txt file
  setwd(paste(dir,'/synthetic_data', sep=""));
  fileConn<-file(paste(dataset, ".txt", sep=""),'w')
  close(fileConn)
  fileConn<-file(paste(dataset, ".txt", sep=""),'a')
  for (i in 1:p){
    for (j in 1:n){
      if (X[i,j]>0){
        writeLines(paste(j,i,X[i,j]), fileConn)
      }
    }
  }
  close(fileConn)
  
  writeMat(paste(dataset, ".A", sep=""), A=A0)
  
  #Generate synthetic vocabulary file, which is required for the Recover method.
  fileConn<-file(paste(dataset, '.vocab.txt', sep=""),'w')
  close(fileConn)
  fileConn<-file(paste(dataset, '.vocab.txt', sep=""),'a')
  for (j in 1:p){
    writeLines(toString(j), fileConn)
  }
  close(fileConn)
}

# synthetic_generation_5: No anchor words + heterogeneity(power law)
synthetic_generation_5 <- function(dir, dataset, p, n, N, K, m_p, m_n, decay, stopnum, W_opt = "puredoc"){
  set.seed(seed)
  #Generate A0
  A0 <- matrix(0,p,K)
  for (i in 1:p){
    A0[i,] <- rexp(K,(i+stopnum)^1.07)
  }
  
  colmax <- apply(A0, 1, max)
  for (k in 1:K){
    anchor_k_index <- sample(which(A0[,k]==colmax), m_p)
    A0[anchor_k_index, -k] <- A0[anchor_k_index, -k]*decay
  }
  
  A0 <- apply(A0, 2, function(x) x/sum(x))

  #Generate W0
  if (W_opt == "puredoc"){
    W0 <- matrix(0,K,n)
    if (m_n>0){
      for (k in 1:K){
        W0[k,((k-1)*m_n+1):(k*m_n)] <- rep(1,m_n)
      }
    }
    W0[,(K*m_n+1):n] <- apply(matrix(runif((n-K*m_n)*K), K, n-K*m_n), 2, function(x) x/sum(x))
  }else if (W_opt == "nopuredoc"){
    W0 <- apply(matrix(runif(n*K), K, n), 2, function(x) x/sum(x))
  }else{
    stop('Incorrect imput for W_opt')
  }
  
  #Generate D0
  D0 <- A0%*%W0
  
  if (length(N)!=n){
    N <- rep(N, n)
  }
  #Generate X
  X <- matrix(0,p,n)
  for (j in 1:n){
    X[,j] <- as.numeric(rmultinom(1, N[j], D0[,j]))
  }
  
  #Compute D
  D <- t(t(X)/N)
  
  #write the sparse semi-synthetic data matrix into txt file
  setwd(paste(dir,'/synthetic_data', sep=""));
  fileConn<-file(paste(dataset, ".txt", sep=""),'w')
  close(fileConn)
  fileConn<-file(paste(dataset, ".txt", sep=""),'a')
  for (i in 1:p){
    for (j in 1:n){
      if (X[i,j]>0){
        writeLines(paste(j,i,X[i,j]), fileConn)
      }
    }
  }
  close(fileConn)
  
  writeMat(paste(dataset, ".A", sep=""),A=A0)
  
  #Generate synthetic vocabulary file, which is required for the Recover method.
  fileConn<-file(paste(dataset, '.vocab.txt', sep=""),'w')
  close(fileConn)
  fileConn<-file(paste(dataset, '.vocab.txt', sep=""),'a')
  for (j in 1:p){
    writeLines(toString(j), fileConn)
  }
  close(fileConn)
}

