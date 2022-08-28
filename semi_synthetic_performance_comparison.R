library(R.matlab)
library("dplyr")
library("ggplot2")
dir <- ''
source('./score_functions.R')
source('./preprocessing.R')
source('./semi_synthetic_generation_functions.R')
source('./SCORE/code/semi_synthetic_data_implementation_SCORE.R')
source('./LDA/code/semi_synthetic_data_implementation_LDA.R')

## For semi-synthetic_generation
# base case K=5, N=NULL or 2000
# change K to be {3, 5, 8, 12} while fix N=ULL
# change N to be {100, 200, 500,1000, 2000} while fix K=5
K_seq <- c(3, 4, 5, 6, 7, 8, 9, 12)
N_seq <- c(100, 200, 500, 1000, 2000)

# For nips.proc, generate semi-synthetic data using the following settings
# K = 6, seed = 1
# K = 7, seed = 1
# K = 8, seed = 9
# K = 9, seed = 9

# dataset_real <- 'nips.proc'
dataset_real <- 'ap'
words_percent <- 0.5
docs_percent <- 0.95
dataset <- paste("prepro_", words_percent,"_",docs_percent,'_', dataset_real, sep="")

# generate all semi_synthetic data first
preprocessing(dir, dataset_real, words_percent, docs_percent)

for (K in K_seq){
  print(paste("K = ", K, ", N = NULL"))
  semi_synthetic_generation(dir, dataset, K, seed = 0)
}
for (N in N_seq){
  print(paste("K = 5, N = ", N))
  semi_synthetic_generation(dir, dataset, 5, N, seed = 0)
}

# learn all norm_SCORE and LDA A_hat
# change K
for (K in c(12)){
  A_synth <- readMat(paste(dir,'/semi_synthetic_data/semi_synthetic_',K,'_',dataset,'.A',sep=''))$A
  data <- read.table(paste(dir,'/semi_synthetic_data/semi_synthetic_',K,'_',dataset,'.txt',sep=''),header=FALSE,sep=" ")
  for (seed in seq(17,19)){
    print(c(seed,K))
    # # norm_SCORE with SVS
    # obj <- semi_synthetic_data_implementation_norm_SCORE(data, K, ceiling(1.5*K), 10*K, seed = seed, VHMethod = 'SVS')
    # A_hat_norm_SCORE <- obj$A_hat
    # time_norm_SCORE <- obj$run_time
    # write.table(A_hat_norm_SCORE,
    #             paste(dir, '/SCORE/results/',seed,'/SCORE.semi_synthetic_',K,'_',dataset,'.A', sep=''))
    # print(paste("norm_SCORE runtime: ",as.numeric(time_norm_SCORE[1]), ", error: ", error2_A(A_synth, A_hat_norm_SCORE),sep=''))

    # norm_SCORE with SP
    obj <- semi_synthetic_data_implementation_norm_SCORE(data, K, ceiling(1.5*K), 10*K, seed = seed, VHMethod = 'SP')
    A_hat_norm_SCORE <- obj$A_hat
    time_norm_SCORE <- obj$run_time
    write.table(A_hat_norm_SCORE,
                paste(dir, '/SCORE/results/',seed,'/SCORE_SP.semi_synthetic_',K,'_',dataset,'.A', sep=''))
    print(paste("norm_SCORE SP runtime: ",as.numeric(time_norm_SCORE[1]), ", error: ", error2_A(A_synth, A_hat_norm_SCORE),sep=''))
  
    # norm_SCORE with SVS-SP
    obj <- semi_synthetic_data_implementation_norm_SCORE(data, K, ceiling(1.5*K), 10*K, seed = seed, VHMethod = 'SVS-SP')
    A_hat_norm_SCORE <- obj$A_hat
    time_norm_SCORE <- obj$run_time
    write.table(A_hat_norm_SCORE,
                paste(dir, '/SCORE/results/',seed,'/SCORE_SVS_SP.semi_synthetic_',K,'_',dataset,'.A', sep=''))
    print(paste("norm_SCORE SP runtime: ",as.numeric(time_norm_SCORE[1]), ", error: ", error2_A(A_synth, A_hat_norm_SCORE),sep=''))
  
    # LDA
    obj <- semi_synthetic_data_implementation_LDA(data, K, seed = seed)
    A_hat_LDA <- obj$A_hat
    time_LDA <- obj$run_time
    write.table(A_hat_LDA,
                paste(dir, '/LDA/results/',seed,'/LDA.semi_synthetic_',K,'_',dataset,'.A', sep=''))
    print(paste("LDA runtime: ",as.numeric(time_LDA[1]), ", error: ", error2_A(A_synth, A_hat_LDA),sep=''))
  }
}

# change N
K <- 5
for (N in N_seq){
  A_synth <- readMat(paste(dir,'/semi_synthetic_data/semi_synthetic_5_',N,'_',dataset,'.A',sep=''))$A
  data <- read.table(paste(dir,'/semi_synthetic_data/semi_synthetic_5_',N,'_',dataset,'.txt',sep=''),header=FALSE,sep=" ")
  for (seed in seq(5,19)){
    print(c(seed,N))
    # norm_SCORE
    obj <- semi_synthetic_data_implementation_norm_SCORE(data, K, ceiling(1.5*K), 10*K, seed = seed, VHMethod = 'SVS')
    A_hat_norm_SCORE <- obj$A_hat
    time_norm_SCORE <- obj$run_time
    write.table(A_hat_norm_SCORE,
                paste(dir, '/SCORE/results/',seed,'/SCORE.semi_synthetic_5_',N,'_',dataset,'.A', sep=''))
    print(paste("norm_SCORE runtime: ",as.numeric(time_norm_SCORE[1]), ", error: ", error2_A(A_synth, A_hat_norm_SCORE),sep=''))
    
    # norm_SCORE with SP
    obj <- semi_synthetic_data_implementation_norm_SCORE(data, K, ceiling(1.5*K), 10*K, seed = seed, VHMethod = 'SP')
    A_hat_norm_SCORE <- obj$A_hat
    time_norm_SCORE <- obj$run_time
    write.table(A_hat_norm_SCORE,
                paste(dir, '/SCORE/results/',seed,'/SCORE_SP.semi_synthetic_5_',N,'_',dataset,'.A', sep=''))
    print(paste("norm_SCORE SP runtime: ",as.numeric(time_norm_SCORE[1]), ", error: ", error2_A(A_synth, A_hat_norm_SCORE),sep=''))
    
    # norm_SCORE with SVS-SP
    obj <- semi_synthetic_data_implementation_norm_SCORE(data, K, ceiling(1.5*K), 10*K, seed = seed, VHMethod = 'SVS-SP')
    A_hat_norm_SCORE <- obj$A_hat
    time_norm_SCORE <- obj$run_time
    write.table(A_hat_norm_SCORE,
                paste(dir, '/SCORE/results/',seed,'/SCORE_SVS_SP.semi_synthetic_5_',N,'_',dataset,'.A', sep=''))
    print(paste("norm_SCORE SVS-SP runtime: ",as.numeric(time_norm_SCORE[1]), ", error: ", error2_A(A_synth, A_hat_norm_SCORE),sep=''))
  
    # LDA
    obj <- semi_synthetic_data_implementation_LDA(data, K, seed = seed)
    A_hat_LDA <- obj$A_hat
    time_LDA <- obj$run_time
    write.table(A_hat_LDA,
                paste(dir, '/LDA/results/',seed,'/LDA.semi_synthetic_5_',N,'_',dataset,'.A', sep=''))
    print(paste("LDA runtime: ",as.numeric(time_LDA[1]), ", error: ", error2_A(A_synth, A_hat_LDA),sep=''))
  }
}



# load and evaluate the errors
df <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df) <- c('method', 'change', 'value', 'seed', 'error')

# change K
for (K in K_seq){
  data <- readMat(paste(dir,'/semi_synthetic_data/semi_synthetic_',K,'_',dataset,'.A',sep=''))
  A_synth <- data$A
  for (seed in seq(0,19)){
    # # norm_SCORE
    # A_hat_norm_SCORE <- read.table(paste(dir, '/SCORE/results/',seed,'/SCORE.semi_synthetic_',K,'_',dataset,'.A', sep=''))
    # df[nrow(df) + 1,] <- c('T-SCORE-SVS', 'K', K, seed, error2_A(A_synth, A_hat_norm_SCORE))

    # norm_SCORE with SP
    A_hat_norm_SCORE <- read.table(paste(dir, '/SCORE/results/',seed,'/SCORE_SP.semi_synthetic_',K,'_',dataset,'.A', sep=''))
    df[nrow(df) + 1,] <- c('T-SCORE-SP', 'K', K, seed, error2_A(A_synth, A_hat_norm_SCORE))
    
    # norm_SCORE with SVS-SP
    A_hat_norm_SCORE <- read.table(paste(dir, '/SCORE/results/',seed,'/SCORE_SVS_SP.semi_synthetic_',K,'_',dataset,'.A', sep=''))
    df[nrow(df) + 1,] <- c('T-SCORE', 'K', K, seed, error2_A(A_synth, A_hat_norm_SCORE))

    # LDA
    A_hat_LDA <- read.table(paste(dir, '/LDA/results/',seed,'/LDA.semi_synthetic_',K,'_',dataset,'.A', sep=''))
    df[nrow(df) + 1,] <- c('LDA', 'K', K, seed, error2_A(A_synth, A_hat_LDA))
    
    # pLSI
    A_hat_plsa <- read.table(paste(dir, '/plsa/results/',seed,'/plsa.semi_synthetic_',K,'_',dataset,'.A', sep=''))
    df[nrow(df) + 1,] <- c('EM', 'K', K, seed, error2_A(A_synth, A_hat_plsa))
    
    # pLSI_SCORE_SVS_SP
    A_hat_plsa <- read.table(paste(dir, '/plsa/results/',seed,'/plsa_SCORE_SVS_SP.semi_synthetic_',K,'_',dataset,'.A', sep=''))
    df[nrow(df) + 1,] <- c('T-SCORE+EM', 'K', K, seed, error2_A(A_synth, A_hat_plsa))
    
    # Recover
    A_hat_Recover <- read.table(paste(dir, '/Recover/results/',seed,'/Recover.semi_synthetic_',K,'_',dataset,'.A', sep=''))
    df[nrow(df) + 1,] <- c('AWR', 'K', K, seed, error2_A(A_synth, A_hat_Recover))
  }

  # TSVD
  A_hat_TSVD <- readMat(paste(dir, '/TSVD/results/TSVD.semi_synthetic_',K,'_',dataset,'.A', sep=''))
  A_hat_TSVD <-A_hat_TSVD[[1]]
  df[nrow(df) + 1,] <- c('TSVD', 'K', K, NA, error2_A(A_synth, A_hat_TSVD))
}

# change N
for (N in N_seq){
  data <- readMat(paste(dir,'/semi_synthetic_data/semi_synthetic_5_',N,'_',dataset,'.A',sep=''))
  A_synth <- data$A
  for (seed in seq(0,19)){
    # # norm_SCORE
    # A_hat_norm_SCORE <- read.table(paste(dir, '/SCORE/results/',seed,'/SCORE.semi_synthetic_5_',N,'_',dataset,'.A', sep=''))
    # df[nrow(df) + 1,] <- c('T-SCORE-SVS', 'N', N, seed, error2_A(A_synth, A_hat_norm_SCORE))
    
    # norm_SCORE with SP
    A_hat_norm_SCORE <- read.table(paste(dir, '/SCORE/results/',seed,'/SCORE_SP.semi_synthetic_5_',N,'_',dataset,'.A', sep=''))
    df[nrow(df) + 1,] <- c('T-SCORE-SP', 'N', N, seed, error2_A(A_synth, A_hat_norm_SCORE))
    
    # norm_SCORE with SVS-SP
    A_hat_norm_SCORE <- read.table(paste(dir, '/SCORE/results/',seed,'/SCORE_SVS_SP.semi_synthetic_5_',N,'_',dataset,'.A', sep=''))
    df[nrow(df) + 1,] <- c('T-SCORE', 'N', N, seed, error2_A(A_synth, A_hat_norm_SCORE))

    # LDA
    A_hat_LDA <- read.table(paste(dir, '/LDA/results/',seed,'/LDA.semi_synthetic_5_',N,'_',dataset,'.A', sep=''))
    df[nrow(df) + 1,] <- c('LDA', 'N', N, seed, error2_A(A_synth, A_hat_LDA))
    
    # pLSI
    A_hat_plsa <- read.table(paste(dir, '/plsa/results/',seed,'/plsa.semi_synthetic_5_',N,'_',dataset,'.A', sep=''))
    df[nrow(df) + 1,] <- c('EM', 'N', N, seed, error2_A(A_synth, A_hat_plsa))
    
    # pLSI_SCORE_SVS_SP
    A_hat_plsa <- read.table(paste(dir, '/plsa/results/',seed,'/plsa_SCORE_SVS_SP.semi_synthetic_5_',N,'_',dataset,'.A', sep=''))
    df[nrow(df) + 1,] <- c('T-SCORE+EM', 'N', N, seed, error2_A(A_synth, A_hat_plsa))
    
    # Recover
    A_hat_Recover <- read.table(paste(dir, '/Recover/results/',seed,'/Recover.semi_synthetic_5_',N,'_',dataset,'.A', sep=''))
    df[nrow(df) + 1,] <- c('AWR', 'N', N, seed, error2_A(A_synth, A_hat_Recover))
  }
  
  # TSVD
  A_hat_TSVD <- readMat(paste(dir, '/TSVD/results/TSVD.semi_synthetic_5_',N,'_',dataset,'.A', sep=''))
  A_hat_TSVD <-A_hat_TSVD[[1]]
  df[nrow(df) + 1,] <- c('TSVD', 'N', N, NA, error2_A(A_synth, A_hat_TSVD))
}

df['value'] <- as.numeric(df[['value']])
df['seed'] <- as.numeric(df[['seed']])
df['error'] <- as.numeric(df[['error']])
# df['error'] <- log(as.numeric(df[['error']]))

median(subset(df, (method == 'T-SCORE' & K == 5))[,'error'])


# plot change K
pdf(file = paste("/Users/greenlink/Desktop/Some_related_existing_methods/semi_synthetic_data/", dataset_real, ".K.20rep.pdf", sep = ''),
    width = 3, height = 3)
methods <- c("T-SCORE","LDA","EM","AWR","TSVD")
dfK <- df[(df['change'] == 'K'),]
dfK <- dfK[dfK[['method']] %in% methods,]
data_group <- dfK %>%
  group_by(method, value) %>%
  dplyr::summarize(mean = mean(error), sd = sd(error)) %>% 
  as.data.frame()
data_group[is.na(data_group)] <- 0
data_group$method <- factor(data_group$method, levels=methods)
p<- ggplot(data_group, aes(x=value, y=mean, group=method, color=method)) + 
  geom_line() +
  geom_point()+
  # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
  #               position=position_dodge(0.05))+
  geom_errorbar(aes(ymin=mean, ymax=mean), width=.2,
                position=position_dodge(0.05))+
  xlab('K')+
  # theme(legend.position = c(0.22, 0.67))
  theme(legend.position = "none")
print(p)
dev.off()

# plot change N
pdf(file = paste("/Users/greenlink/Desktop/Some_related_existing_methods/semi_synthetic_data/", dataset_real, ".N.20rep.pdf", sep = ''),
    width = 3, height = 3)
methods <- c("T-SCORE","LDA","EM","AWR","TSVD")
dfN <- df[(df['change'] == 'N'),]
dfN <- dfN[dfN[['method']] %in% methods,]
data_group <- dfN %>%
  group_by(method, value) %>%
  dplyr::summarize(mean = mean(error), 
                   sd = sd(error),
                   median = median(error)) %>% 
  as.data.frame()
data_group[is.na(data_group)] <- 0
data_group$method <- factor(data_group$method, levels=methods)
p<- ggplot(data_group, aes(x=value, y=mean, group=method, color=method)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=median, ymax=median), width=.2,
                position=position_dodge(0.05))+
  xlab('N')+
  # theme(legend.position = c(0.75, 0.67))
  theme(legend.position = "none")
print(p)
dev.off()

# dfN$value <- as.factor(dfN$value)
# ggplot(dfN[(dfN['change'] == 'N'),], aes(x=value, y=error, fill=method)) + 
#  geom_boxplot()

