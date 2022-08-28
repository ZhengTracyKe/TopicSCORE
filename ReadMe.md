# TopicSCORE
Code and data for the paper "Using SVD for Topic Modeling" (https://arxiv.org/abs/1704.07016).

The main file is TSCOREfunctions.R. It contains all the R functions needed for running the Topic-SCORE algorithm. Given the word-document-frequency matrix $D$ and the number of topics $K$, use the command $\texttt{score(K, K0, m, D)}$ to obtain the estimated topic matrix $\hat{A}$ and other output. Here, $\texttt{K0}$ and $\texttt{m}$ are two tuning parameters, whose default choices are given in the paper.

We also include the Associated Press (AP) data set and the Statistical Literature Abstracts (SLA) data set analyzed in the paper. 


