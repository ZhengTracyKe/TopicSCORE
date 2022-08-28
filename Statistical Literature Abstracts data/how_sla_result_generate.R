stopwords <- unique(c('a','about','above','after','again','against','all','am',
                               'an','and','any','are','arent','as','at','be','because','been','before','being',
                               'below','between','both','but','by','cant','cannot','could','couldnt','did','didnt','do',
                               'does','doesnt','doing','dont','down','during','each','few','for','from','further','had',
                               'hadnt','has','hasnt','have','havent','having','he','hed','hes','her','here','heres',
                               'hers','herself','him','himself','his','how','hows','i','id','im','ive','if','in',
                               'into','is','isnt','it','its','itself','lets','me','more','most','mustnt','my',
                               'myself','no','nor','not','of','off','on','once','only','or','other','ought',
                               'our','ours','ourselves','out','over','own','same','shant','she','shes','should',
                               'shouldnt','so','some','such','than','that','thats','the','their','theirs','them',
                               'themselves','then','there','theres','these','they','theyd','theyll','theyre','theyve',
                               'this','those','through','to','too','under','until','up','very','was','wasnt','we',
                               'weve','were','werent','what','whats','when','whens','where','wheres','which','while',
                               'who','whos','whom','why','whys','with','wont','would','wouldnt','you','youd','youll',
                               'youre','youve','your','yours','yourself','yourselves','percent','proposed', 'using', 
                               'also', 'results', 'approach', 'study', 'problem', 'propose', 'consider', 'often', 
                               'however', 'several', 'derived', 'result', 'resulting', 'obtain', 'literature','statements',
                               letters))

# Read in data
data <- read.table('abstract.txt')
data <- as.matrix(data)
vocab <- read.table('abstract.vocab.txt',colClasses = "character")[[1]]
D_count <- matrix(0,max(data[,2]),max(data[,1]))
for (t in 1:dim(data)[1]){
  D_count[data[t,2], data[t,1]] <-data[t,3]
}

p <- dim(D_count)[1]
n <- dim(D_count)[2]

# Set some thresholds
w_num <- 3000        #number of words to keep
d_percent <- 0.6     #percentage of docs to keep
Mquantile <- 1     #Truncate quantile of M

#Only keep d_percent% longest documents
doc_count <- colSums(D_count)
doc_keep <- which(rank(-doc_count, ties.method = 'first')<=round(d_percent*n))
D_count <- D_count[,doc_keep]

#Only keep top w_num most frequent words
word_count <- rowSums(D_count)
word_keep <- which(rank(-word_count, ties.method = 'first')<=w_num)
D_count <- D_count[word_keep,]
vocab <- vocab[word_keep]

#Generate D from D_count
D <- apply(D_count, 2, function(x) x/sum(x))

## Normalized SCORE
p <- dim(D)[1]
n <- dim(D)[2]
M <- rowMeans(D)
M_trunk <- pmin(M,quantile(M,Mquantile))

obj <- svd(sqrt(M_trunk^(-1))*D)
Xi <- obj$u