library(R.matlab)
# library("dplyr")
# library("ggplot2")
dir <- ''
source('./score_functions.R')
source('/real_data/preprocessing.R')
source('/synthetic_data/synthetic_generation_functions.R')
source('/SCORE/code/synthetic_data_implementation_SCORE.R')
source('/LDA/code/synthetic_data_implementation_LDA.R')
options(scipen=999)

## For synthetic_generation_1: Anchor words + homogeneity
# base case p=1000,n=1000,N=2000,K=5,m_p=p/100,delta_p=1/p,m_n=n/100
# change p to be c(100, 500, 1000, 2000, 4000) while fix the remaining
# change n to be c(500, 1000, 2000, 5000,10000) while fix the remaining
# change N to be c(500, 1000, 2000, 5000, 10000) while fix the remaining
# change K to be c(3, 5, 10, 15, 20) while fix the remaining
# change m_p to be c(p/500, p/200, p/100, p/50, p/10) while fix the remaining
# change delta_p to be c(1/(10*p), 1/(5*p), 1/p, 5/p, 10/p) while fix the remaining
# change m_n to be c(0,n/1000, n/500, n/200, n/100, n/50, n/10) while fix the remaining

set_base <- function(){
  p <<- 1000
  n <<- 1000
  N <<- 2000
  K <<- 5
  m_p <<- p/100
  delta_p <<- 1/p
  m_n <<- n/100
  W_opt <<- 'puredoc'
  print(paste("set base: ", p,n,N,K,m_p,delta_p,m_n, sep=", "))
}

set_base()
p_seq <- c(100, 500, 1000, 2000, 4000)
n_seq <- c(500, 1000, 2000, 5000,10000)
N_seq <- c(500, 1000, 2000, 5000, 10000)
K_seq <- c(3, 5, 10, 15, 20)
m_p_seq <- c(p/500, p/200, p/100, p/50, p/10)
delta_p_seq <- c(1/(10*p), 1/(5*p), 1/p, 5/p, 10/p)
m_n_seq <- c(0,n/1000, n/500, n/200, n/100, n/50, n/10)

get_dataset <- function(){
  # change this function for different synthetic experiment settings
  return(paste("synthetic_1_", p,',', n,',', mean(N),'-',length(unique(N)),',',
               K, ',', m_p, ',', delta_p, ',', m_n,',', W_opt, sep=""))
}

fitOnce <- function(){
  dataset <- get_dataset()
  A_synth <- readMat(paste(dir,'/synthetic_data/',dataset,'.A',sep=''))$A
  data <- read.table(paste(dir,'/synthetic_data/',dataset,'.txt',sep=''),header=FALSE,sep=" ")
  for (seed in seq(0,19)){
    print(c(seed,K))
    # norm_SCORE with SVS
    # obj <- synthetic_data_implementation_norm_SCORE(data, K, ceiling(1.5*K), 10*K, seed = seed, VHMethod = 'SVS')
    # A_hat_norm_SCORE <- obj$A_hat
    # time_norm_SCORE <- obj$run_time
    # write.table(A_hat_norm_SCORE,
    #             paste(dir, '/SCORE/results/',seed,'/SCORE.',dataset,'.A', sep=''))
    # print(paste("norm_SCORE runtime: ",as.numeric(time_norm_SCORE[1]), ", error: ", error2_A(A_synth, A_hat_norm_SCORE),sep=''))
    
    # norm_SCORE with SP
    obj <- synthetic_data_implementation_norm_SCORE(data, K, ceiling(1.5*K), 10*K, seed = seed, VHMethod = 'SP')
    A_hat_norm_SCORE <- obj$A_hat
    time_norm_SCORE <- obj$run_time
    write.table(A_hat_norm_SCORE,
                paste(dir, '/SCORE/results/',seed,'/SCORE_SP.',dataset,'.A', sep=''))
    print(paste("norm_SCORE SP runtime: ",as.numeric(time_norm_SCORE[1]), ", error: ", error2_A(A_synth, A_hat_norm_SCORE),sep=''))
    
    # norm_SCORE with SVS-SP
    obj <- synthetic_data_implementation_norm_SCORE(data, K, ceiling(1.5*K), 10*K, seed = seed, VHMethod = 'SVS-SP')
    A_hat_norm_SCORE <- obj$A_hat
    time_norm_SCORE <- obj$run_time
    write.table(A_hat_norm_SCORE,
                paste(dir, '/SCORE/results/',seed,'/SCORE_SVS_SP.',dataset,'.A', sep=''))
    print(paste("norm_SCORE SVS SP runtime: ",as.numeric(time_norm_SCORE[1]), ", error: ", error2_A(A_synth, A_hat_norm_SCORE),sep=''))
    
    # # LDA
    # obj <- synthetic_data_implementation_LDA(data, K, seed = seed)
    # A_hat_LDA <- obj$A_hat
    # time_LDA <- obj$run_time
    # write.table(A_hat_LDA,
    #             paste(dir, '/LDA/results/',seed,'/LDA.',dataset,'.A', sep=''))
    # print(paste("LDA runtime: ",as.numeric(time_LDA[1]), ", error: ", error2_A(A_synth, A_hat_LDA),sep=''))
  }
}

genResultDf <- function(datasets, change, values){
  # load and evaluate the errors
  df <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(df) <- c('method', 'change', 'value', 'seed', 'error')
  for (i in seq(length(datasets))){
    dataset <-  datasets[i]
    value <- values[i]
    data <- readMat(paste(dir,'/synthetic_data/',dataset,'.A',sep=''))
    A_synth <- data$A
    for (seed in seq(0,19)){
      # # norm_SCORE with SP
      # A_hat <- read.table(paste(dir, '/SCORE/results/',seed,'/SCORE_SP.',dataset,'.A', sep=''))
      # df[nrow(df) + 1,] <- c('T-SCORE-SP', change, value, seed, error2_A(A_synth, A_hat))
      
      # norm_SCORE with SVS-SP
      A_hat <- read.table(paste(dir, '/SCORE/results/',seed,'/SCORE_SVS_SP.',dataset,'.A', sep=''))
      df[nrow(df) + 1,] <- c('T-SCORE', change, value, seed, error2_A(A_synth, A_hat))
      
      # # LDA
      # A_hat <- read.table(paste(dir, '/LDA/results/',seed,'/LDA.',dataset,'.A', sep=''))
      # df[nrow(df) + 1,] <- c('LDA', change, value, seed, error2_A(A_synth, A_hat))
      
      # pLSI
      A_hat <- read.table(paste(dir, '/plsa/results/',seed,'/plsa.',dataset,'.A', sep=''))
      df[nrow(df) + 1,] <- c('EM', change, value, seed, error2_A(A_synth, A_hat))
      
      # pLSI_SCORE_SVS_SP
      A_hat <- read.table(paste(dir, '/plsa/results/',seed,'/plsa_SCORE_SVS_SP.',dataset,'.A', sep=''))
      df[nrow(df) + 1,] <- c('T-SCORE+EM', change, value, seed, error2_A(A_synth, A_hat))
      
      # # Recover
      # A_hat <- read.table(paste(dir, '/Recover/results/',seed,'/Recover.',dataset,'.A', sep=''))
      # df[nrow(df) + 1,] <- c('AWR', change, value, seed, error2_A(A_synth, A_hat))
    }
    
    # # TSVD
    # A_hat <- readMat(paste(dir, '/TSVD/results/TSVD.',dataset,'.A', sep=''))
    # A_hat <-A_hat[[1]]
    # df[nrow(df) + 1,] <- c('TSVD', change, value, NA, error2_A(A_synth, A_hat))
  }
  df['value'] <- as.numeric(df[['value']])
  df['seed'] <- as.numeric(df[['seed']])
  df['error'] <- as.numeric(df[['error']])
  # df['error'] <- log(as.numeric(df[['error']]))
  return(df)
}

plotResultDf <- function(df, change, methods = c("T-SCORE-SVS-SP","LDA","pLSI","pLSI_T-SCORE-SVS-SP","AWR","TSVD")){
  library("dplyr")
  library("ggplot2")
  df <- df[(df['change'] == change),]
  df <- df[df[['method']] %in% methods,]
  data_group <- df %>%
    group_by(method, value) %>%
    dplyr::summarize(mean = mean(error), sd = sd(error)) %>% 
    as.data.frame()
  data_group[is.na(data_group)] <- 0
  data_group$method <- factor(data_group$method, levels=methods)
  ggplotObj <- ggplot(data_group, aes(x=value, y=mean, group=method, color=method)) + 
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(0.05))+
    xlab(change)
  return(ggplotObj)
}


# generate synthetic data
set_base()
for (p in p_seq){
  print(paste("p = ", p))
  dataset <- get_dataset()
  synthetic_generation_1(dir, dataset, p, n, N, K, m_p, delta_p, m_n,W_opt,seed = 0)
}
set_base()
for (n in n_seq){
  print(paste("n = ", n))
  dataset <- get_dataset()
  synthetic_generation_1(dir, dataset, p, n, N, K, m_p, delta_p, m_n,W_opt,seed = 0)
}
set_base()
for (N in N_seq){
  print(paste("N = ", N))
  dataset <- get_dataset()
  synthetic_generation_1(dir, dataset, p, n, N, K, m_p, delta_p, m_n,W_opt,seed = 0)
}
set_base()
for (K in K_seq){
  print(paste("K = ", K))
  dataset <- get_dataset()
  synthetic_generation_1(dir, dataset, p, n, N, K, m_p, delta_p, m_n,W_opt,seed = 0)
}
set_base()
for (m_p in m_p_seq){
  print(paste("m_p = ", m_p))
  dataset <- get_dataset()
  synthetic_generation_1(dir, dataset, p, n, N, K, m_p, delta_p, m_n,W_opt,seed = 0)
}
set_base()
for (delta_p in delta_p_seq){
  print(paste("delta_p = ", delta_p))
  dataset <- get_dataset()
  synthetic_generation_1(dir, dataset, p, n, N, K, m_p, delta_p, m_n,W_opt,seed = 0)
}
set_base()
for (m_n in m_n_seq){
  print(paste("m_n = ", m_n))
  dataset <- get_dataset()
  synthetic_generation_1(dir, dataset, p, n, N, K, m_p, delta_p, m_n,W_opt,seed = 0)
}


# learn all norm_SCORE and LDA A_hat
set_base()
for (p in p_seq){
  print(paste("p = ", p))
  fitOnce()
}
set_base()
for (n in n_seq){
  print(paste("n = ", n))
  fitOnce()
}
set_base()
for (N in N_seq){
  print(paste("N = ", N))
  fitOnce()
}
set_base()
for (K in K_seq){
  print(paste("K = ", K))
  fitOnce()
}
set_base()
for (m_p in m_p_seq){
  print(paste("m_p = ", m_p))
  fitOnce()
}
set_base()
for (delta_p in delta_p_seq){
  print(paste("delta_p = ", delta_p))
  fitOnce()
}
set_base()
for (m_n in m_n_seq){
  print(paste("m_n = ", m_n))
  fitOnce()
}



# change p
set_base()
K <- 5
datasets <- c()
for (p in p_seq){
  datasets <- c(datasets, get_dataset())
}
df <- genResultDf(datasets, 'p', p_seq)
pdf(file = "/Users/greenlink/Desktop/Some_related_existing_methods/synthetic_data/p.3methods.20rep.legend.pdf",
    width = 3, height = 3)
plotResultDf(df, 'p', methods = c("EM","T-SCORE","T-SCORE+EM"))+
  theme(legend.position = c(0.7, 0.6))
  # theme(legend.position = "none")
dev.off()

# change n
set_base()
K <- 5
datasets <- c()
for (n in n_seq){
  datasets <- c(datasets, get_dataset())
}
df <- genResultDf(datasets, 'n', n_seq)
pdf(file = "/Users/greenlink/Desktop/Some_related_existing_methods/synthetic_data/n.3methods.20rep.legend.pdf",
    width = 3, height = 3)
plotResultDf(df, 'n', methods = c("EM","T-SCORE","T-SCORE+EM"))+
  theme(legend.position = c(0.7, 0.5))
  # theme(legend.position = "none")
dev.off()


