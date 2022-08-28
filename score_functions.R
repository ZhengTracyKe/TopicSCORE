score <- function(K, K0, m, D, scatterplot=FALSE){
  library('rARPACK')
  p <- dim(D)[1]
  obj <- svds(D,K)
  Xi <- obj$u
  
  #Step 1
  Xi[,1] <- abs(Xi[,1])
  R <- apply(as.matrix(Xi[,2:K]),2,function(x) x/Xi[,1])
  
  #Step 2
  vertices_est_obj <- vertices_est(R,K0,m)
  V <- vertices_est_obj$V
  theta <- vertices_est_obj$theta
  
  if (scatterplot){
    par(mar=c(1,1,1,1))
    plot(R[,1],R[,2])
    points(V[,1],V[,2],col=2,lwd=5)
  }
  
  #Step 3
  Pi <- cbind(R, rep(1,p))%*%solve(cbind(V,rep(1,K)))
  Pi <- pmax(Pi,matrix(0,dim(Pi)[1],dim(Pi)[2]))
  temp <- rowSums(Pi)
  Pi <- apply(Pi,2,function(x) x/temp)
  
  #Step 4
  A_hat <- Xi[,1]*Pi
  
  #Step 5
  temp <- colSums(A_hat)
  A_hat <- t(apply(A_hat,1,function(x) x/temp))
  
  return(list(A_hat=A_hat, R=R,V=V, Pi=Pi, theta=theta))
}

norm_score_N <- function(K, K0, m, N, D, scatterplot=FALSE){
  library('rARPACK')
  p <- dim(D)[1]
  n <- dim(D)[2]
  M <- rowSums(D/t(matrix(rep(N,p),n,p)))
  obj <- svds(sqrt(M^(-1))*D, K)
  Xi <- obj$u
  
  #Step 1
  Xi[,1] <- abs(Xi[,1])
  R <- apply(Xi[,2:K],2,function(x) x/Xi[,1])
  
  #Step 2
  vertices_est_obj <- vertices_est(R,K0,m)
  V <- vertices_est_obj$V
  theta <- vertices_est_obj$theta
  
  if (scatterplot){
    par(mar=c(1,1,1,1))
    plot(R[,1],R[,2])
    points(V[,1],V[,2],col=2,lwd=5)
  }
  
  #Step 3
  Pi <- cbind(R, rep(1,p))%*%solve(cbind(V,rep(1,K)))
  Pi <- pmax(Pi,matrix(0,dim(Pi)[1],dim(Pi)[2]))
  temp <- rowSums(Pi)
  Pi <- apply(Pi,2,function(x) x/temp)
  
  #Step 4
  A_hat <- sqrt(M)*Xi[,1]*Pi
  #A_hat <- diag(sqrt(M))%*%diag(Xi[,1])%*%Pi
  
  #Step 5
  temp <- colSums(A_hat)
  A_hat <- t(apply(A_hat,1,function(x) x/temp))
  
  return(list(A_hat=A_hat, R=R,V=V, theta=theta))
}

norm_score <- function(K, K0, m, D, Mquantile=0, scatterplot=FALSE, VHMethod = 'SVS'){
  library('rARPACK')
  library('nnls')
  p <- dim(D)[1]
  n <- dim(D)[2]
  M <- rowMeans(D)
  M_trunk <- pmin(M,quantile(M,Mquantile))
  
  obj <- svds(sqrt(M_trunk^(-1))*D, K)
  Xi <- obj$u
  
  #Step 1
  Xi[,1] <- abs(Xi[,1])
  R <- apply(Xi[,2:K],2,function(x) x/Xi[,1])
  
  #Step 2
  if (VHMethod == 'SVS'){
    vertices_est_obj <- vertices_est(R,K0,m)
    V <- vertices_est_obj$V
    theta <- vertices_est_obj$theta
  } else if (VHMethod == 'SP'){
    vertices_est_obj <- successiveProj(R, K)
    V <- vertices_est_obj$V
    theta <- NULL
  } else if (VHMethod == 'SVS-SP'){
    vertices_est_obj <- vertices_est_SP(R, m)
    V <- vertices_est_obj$V
    theta <- NULL
  }
  
  
  if (scatterplot){
    print('wrong!')
    par(mar=c(1,1,1,1))
    plot(R[,1],R[,2])
    points(V[,1],V[,2],col=2,lwd=5)
  }
  
  #Step 3
  Pi <- cbind(R, rep(1,p))%*%solve(cbind(V,rep(1,K)))
  Pi <- pmax(Pi,0)
  temp <- rowSums(Pi)
  Pi <- apply(Pi,2,function(x) x/temp)
  
  #Step 4
  A_hat <- sqrt(M_trunk)*Xi[,1]*Pi
  
  #Step 5
  temp <- colSums(A_hat)
  A_hat <- t(apply(A_hat,1,function(x) x/temp))
  
  return(list(A_hat=A_hat, R=R,V=V, Pi=Pi, theta=theta))
}

score_W <- function(K, K0, m, D){
  library('rARPACK')
  p <- dim(D)[1]
  obj <- svds(D,K)
  Xi <- obj$v
  
  #Step 1
  Xi[,1] <- abs(Xi[,1])
  R <- apply(as.matrix(Xi[,2:K]),2,function(x) x/Xi[,1])
  
  #Step 2
  vertices_est_obj <- vertices_est(R,K0,m)
  V <- vertices_est_obj$V
  
  #Step 3
  Pi <- cbind(R, rep(1,n))%*%solve(cbind(V,rep(1,K)))
  Pi <- pmax(Pi,0)
  temp <- rowSums(Pi)
  Pi <- apply(Pi,2,function(x) x/temp)
  
  #Step 4
  W_tilde <- Pi*Xi[,1]
  W_tilde_norm <- apply(W_tilde, 2, function(x){x-mean(x)})
  Q_direction <- eigen(t(W_tilde_norm)%*%W_tilde_norm)$vectors[,K]
  
  #Step 5
  W_hat <- apply(W_tilde*abs(Q_direction), 1, function(x){x/sum(x)})
  
  #Step 6: nnls to recover A
  library(nnls)
  A_hat <- matrix(0, p, K)
  for (i in 1:p){
    A_hat[i,] <- nnls(t(W_hat), D[i,])$x
  }
  for (k in 1:K){
    A_hat[,k] <- A_hat[,k]/sum(A_hat[,k])
  }
  
  return(list(A_hat=A_hat, W_hat=W_hat, R=R,V=V, Pi=Pi))
}

debias_score_N <- function(K, K0, m, N, D, scatterplot=FALSE){
  library('rARPACK')
  p <- dim(D)[1]
  n <- dim(D)[2]
  M <- rowSums(D/t(matrix(rep(N,p),n,p)))
  obj <- eigs(D%*%t(D) - diag(M), K)
  Xi <- obj$vectors
  
  #Step 1
  Xi[,1] <- abs(Xi[,1])
  R <- apply(as.matrix(Xi[,2:K]),2,function(x) x/Xi[,1])
  
  #Step 2
  vertices_est_obj <- vertices_est(R,K0,m)
  V <- vertices_est_obj$V
  theta <- vertices_est_obj$theta
  
  if (scatterplot){
    par(mar=c(1,1,1,1))
    plot(R[,1],R[,2])
    points(V[,1],V[,2],col=2,lwd=5)
  }
  
  #Step 3
  Pi <- cbind(R, rep(1,p))%*%solve(cbind(V,rep(1,K)))
  Pi <- pmax(Pi,matrix(0,dim(Pi)[1],dim(Pi)[2]))
  temp <- rowSums(Pi)
  Pi <- apply(Pi,2,function(x) x/temp)
  
  #Step 4
  A_hat <- Xi[,1]*Pi
  
  #Step 5
  temp <- colSums(A_hat)
  A_hat <- t(apply(A_hat,1,function(x) x/temp))
  
  return(list(A_hat=A_hat, R=R,V=V,Pi=Pi, theta=theta))
}

debias_score <- function(K, K0, m, N, D, scatterplot=FALSE){
  library('rARPACK')
  p <- dim(D)[1]
  n <- dim(D)[2]
  M <- rowMeans(D)
  obj <- eigs(D%*%t(D) - n/mean(N)*diag(M), K)
  Xi <- obj$vectors
  
  #Step 1
  Xi[,1] <- abs(Xi[,1])
  R <- apply(as.matrix(Xi[,2:K]),2,function(x) x/Xi[,1])
  
  #Step 2
  vertices_est_obj <- vertices_est(R,K0,m)
  V <- vertices_est_obj$V
  theta <- vertices_est_obj$theta
  
  if (scatterplot){
    par(mar=c(1,1,1,1))
    plot(R[,1],R[,2])
    points(V[,1],V[,2],col=2,lwd=5)
  }
  
  #Step 3
  Pi <- cbind(R, rep(1,p))%*%solve(cbind(V,rep(1,K)))
  Pi <- pmax(Pi,matrix(0,dim(Pi)[1],dim(Pi)[2]))
  temp <- rowSums(Pi)
  Pi <- apply(Pi,2,function(x) x/temp)
  
  #Step 4
  A_hat <- Xi[,1]*Pi
  
  #Step 5
  temp <- colSums(A_hat)
  A_hat <- t(apply(A_hat,1,function(x) x/temp))
  
  return(list(A_hat=A_hat, R=R,V=V,Pi=Pi,theta=theta))
}

vertices_est <- function(R,K0,m){
  library(quadprog)
  K <- dim(R)[2] + 1

  #Step 2a
  obj <- kmeans(R,m,iter.max=K*100,nstart = K*10)
  
  theta <- as.matrix(obj$centers)
  theta_original <- theta
  plot(R[,1],R[,2])
  points(theta[,1], theta[,2], col=2,lwd=4)
  
  #Step 2b'
  inner <- theta%*%t(theta)
  distance <- diag(inner)%*%t(rep(1,length(diag(inner)))) + rep(1,length(diag(inner)))%*%t(diag(inner)) - 2*inner
  top2 <- which(distance==max(distance),arr.ind=TRUE)[1,]
  theta0 <- as.matrix(theta[top2,])
  theta <- as.matrix(theta[-top2,])
  
  if (K0 > 2){
    for (k0 in 3:K0){
      inner <- theta%*%t(theta)
      distance <- rep(1,k0-1)%*%t(diag(inner))-2*theta0%*%t(theta)
      ave_dist <- colMeans(distance)
      index <- which(ave_dist==max(ave_dist))[1]
      theta0 <- rbind(theta0, theta[index,])
      theta <- as.matrix(theta[-index,])
    }
    theta <- theta0
  }
  
  #Step 2b
  comb <- combn(1:K0, K)
  max_values <- rep(0, dim(comb)[2])
  for (i in 1:dim(comb)[2]){
    for (j in 1:K0){
      max_values[i] <- max(simplex_dist(as.matrix(theta[j,]), as.matrix(theta[comb[,i],])), max_values[i])
    }
  }
  
  min_index <- which(max_values == min(max_values))
  
  #plot(theta[,1],theta[,2])
  #points(theta[comb[,min_index],1],theta[comb[,min_index],2],col=2,pch=2)
  
  return(list(V=theta[comb[,min_index[1]],], theta=theta_original))
}

simplex_dist <- function(theta, V){
  library(Matrix)
  VV <- cbind(diag(rep(1,dim(V)[1]-1)), -rep(1,dim(V)[1]-1))%*%V
  D <- VV%*%t(VV)
  d <- VV%*%(theta-V[dim(V)[1],])
  
  A <- cbind(diag(rep(1,dim(V)[1]-1)), -rep(1,dim(V)[1]-1))
  b0 <- c(rep(0,dim(V)[1]-1),-1)
  
  # D <- matrix(nearPD(D)$mat, nrow(D), ncol(D))
  # D <- nearPD(D)
  obj <- solve.QP(D, d, A, b0)
  return(sum((theta-V[dim(V)[1],]) ^2)+ 2*obj$value)
}

error1_A <- function(A, A_hat){
  library(combinat)
  K <- dim(A)[2]
  
  all_perm <- permn(1:K)
  error <- Inf
  
  for (i in 1:length(all_perm)){
    error <- min(error, sum(colSums(abs(A[,all_perm[[i]]]-A_hat))))
  }
  
  return(error)
}

error2_A <- function(A, A_hat){
  K <- dim(A)[2]
  used <- rep(1,K)
  A_perm <- matrix(0,dim(A)[1],dim(A)[2])
  
  for (k in 1:K){
    dis <- colSums(abs(A-A_hat[,k]))*used
    index <- which(dis == min(dis))
    index <- index[1]
    A_perm[,k] <- A[,index]
    used[index] <- Inf
  }
  
  return(sum(colSums(abs(A_perm-A_hat))))
}

compute_W_from_AD <- function(A_hat, D){
  library(Matrix)
  # can also be used as computing A from W and D
  K <-dim(A_hat)[2]
  n <- dim(D)[2]
  
  W_hat <- matrix(0, K, n)
  M <- rbind(diag(K-1), rep(-1,K-1))
  bM <- diag(K)[,K]
  Dmat <- 2*t(A_hat%*%M)%*%(A_hat%*%M)
  Amat <- t(M)
  bvec <- -bM
  
  AM <- A_hat%*%M
  AbM <- A_hat%*%bM
  for (i in 1:n){
    dvec <- 2*t(D[,i]-AbM)%*%AM
    # Dmat <- matrix(nearPD(Dmat)$mat, nrow(Dmat), ncol(Dmat))
    # Dmat <- nearPD(Dmat)
    qp_sol <- solve.QP(Dmat, dvec, Amat, bvec)$solution
    W_hat[,i] <- c(qp_sol, 1-sum(qp_sol))
  }
  W_hat <- pmax(W_hat,0)
  
  return(W_hat)
}

replaceWithLeastPositive <- function(vec){
  vec[vec<=0] = min(vec[vec>0])
  return(vec)
}

nearPD <- function(mat){
  mat <- (mat+t(mat))/2
  eigenObj <- eigen(mat)
  values <- eigenObj$values
  vectors <- eigenObj$vectors
  values <- replaceWithLeastPositive(values)
  return(vectors%*%diag(values)%*%t(vectors))
}

successiveProj <- function(R, K){
  # succesive projection on rows of R
  n <- dim(R)[1]

  Y <- cbind(rep(1,n),R)
  indexSet = c()
  while (length(indexSet) < K){
    l2Norms <- apply(Y,1,function(x) sqrt(sum(x^2)))
    index <- which(l2Norms == max(l2Norms))
    indexSet <- c(indexSet, index)
    u <- Y[index,] / sqrt(sum(Y[index,]^2))
    Y <- t(apply(Y,1,function(x) x-sum(x*u)*u))
  }
  return(list(V=R[indexSet,], indexSet=indexSet))
}

vertices_est_SP <- function(R,m){
  library(quadprog)
  K <- dim(R)[2] + 1

  obj <- kmeans(R,m,iter.max=K*100,nstart = K*10)
  theta <- as.matrix(obj$centers)
  return(successiveProj(theta, K))
}
  
