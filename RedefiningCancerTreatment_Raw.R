setwd("C:\\Users\\SROY\\Documents\\CodeBase\\Datasets\\Redefining Cancer Treatment")
rm(list=ls())
seedVal = 17869
#load("RDC3.RData")
#save.image("RDC3.RData")

# Import Data
test_variants <- read.csv2('test_variants', sep=',')
train_variants <- read.csv2('training_variants', sep=',')

library('tibble')
library('readr')
library('dplyr')
library('tidyr')
train_txt_dump <- tibble(text = read_lines('training_text', skip = 1))
train_txt <- train_txt_dump %>%
  separate(text, into = c("ID", "txt"), sep = "\\|\\|")
train_txt$ID <- as.integer(train_txt$ID)

test_txt_dump <- tibble(text = read_lines('test_text', skip = 1))
test_txt <- test_txt_dump %>%
  separate(text, into = c("ID", "txt"), sep = "\\|\\|")
test_txt$ID <- as.integer(test_txt$ID)

# Merge data
test <- merge(x = test_txt, y = test_variants, by = 'ID')
train <- merge(x = train_txt, y = train_variants, by='ID')

# Append test with train
test$Class <- -1
train <- rbind(train, test)

# remove variables to save space.
rm(list=c('test','train_txt_dump','train_txt', 'test_txt_dump',
          'test_txt', 'train_variants', 'test_variants'))

# Convert into datatype
str(train)
train$Variation <- as.character(train$Variation)
train$Class <- as.factor(train$Class)
train$Gene <- as.character(train$Gene)

library(ggplot2)
ggplot(train,aes(Class)) + geom_bar()
# Majority is class 7 fowwed by 4
train[is.na(train$Class),]
apply(train, 2, function(x) length(which(is.na(x) | x == "NA")))

################### Natural Language Processing ###############
#
# minimize distance vector, cosine similarity or latent semantic analysis
# tokenize, stemmer, lemmatize, remove stopwords, Clean, detokenize
# CountVectorizer, tfidf, word2vec, doc2vec

# # too many alpha numeric issues. Let's remove all
# library(stringr)
# train$txt <- sapply(train$txt, function(x) str_replace_all(x, "[^[:alnum:]]", " "))

# Use tm packge for cleanup
library(NLP)
library(tm)

train_txtCorp <- Corpus(VectorSource(train$txt))

# Clean data
toSpace <- content_transformer(function(x, pattern) gsub(pattern, " ", x, fixed = FALSE))
# remove url
train_txtCorp <- tm_map(train_txtCorp, toSpace, "(ht|f)tp(s?)\\:\\/\\/[0-9a-zA-Z]([-.\\w]*[0-9a-zA-Z])*(:(0-9)*)*(\\/?)([a-zA-Z0-9\\-\\.\\?\\,\\'\\/\\\\+&amp;%\\$#_]*)?")
# Remove partial URLs:
train_txtCorp <- tm_map(train_txtCorp, toSpace, "[0-9a-zA-Z]*\\.*\\.[0-9a-zA-Z]*")
# Remove href tags
train_txtCorp <- tm_map(train_txtCorp, toSpace, "<a.+/a>")
# remove alpha numeric values
train_txtCorp <- tm_map(train_txtCorp, toSpace, "[^[:alnum:]]")
# Remove stopwords and stem
train_txtCorp <- tm_map(train_txtCorp, content_transformer(tolower))
train_txtCorp <- tm_map(train_txtCorp, removeWords, stopwords("en"))
train_txtCorp <- tm_map(train_txtCorp, stemDocument, language="english")

#** saved 1

# Remove more

# Let's verify corpus what else we should remove
inspect(train_txtCorp[3])
# remove short length words of length 2 or less
removeShortWords <- content_transformer(function(x) gsub('\\b\\w{1,2}\\b' , " ", x))
train_txtCorp <- tm_map(train_txtCorp, removeShortWords)

# analyse numbers
# Let's remove only numbers and not s80n etc. Analyse
removeDigitsOnly <- content_transformer(function(x) gsub('\\b\\d+\\b' , "", x))
train_txtCorp <- tm_map(train_txtCorp, removeDigitsOnly)

# Remove white spaces to organize the text a bit
train_txtCorp <- tm_map(train_txtCorp, stripWhitespace)

inspect(train_txtCorp[10])

#train_txtCorp <- tm_map(train_txtCorp, removeNumbers)
# remove words since lemmatize is not working

#-saved 2

#need this?
# train_txtCorp <- tm_map(train_txtCorp, PlainTextDocument)
# summary(train_txtCorp)

# Let's convert into dtm without idf for analyzing
train_dtm <- DocumentTermMatrix(train_txtCorp)
train_dtm
inspect(train_dtm)

# Frequency.  occur at least 100000 times
findFreqTerms(train_dtm, lowfreq = 100000)
#findFreqTerms(removeSparseTerms(train_dtm, 0.2), lowfreq = 100000)
#doesn't help much

# show and shown occur so many times. Lemmatizing is required.
# So many useless terms in the top results:
# also, fig, use, includ ... needs filtering

# We can plot the frequency occurance as the whole freq matrix doesn't fit
#highFreqTerms <- data.frame(findMostFreqTerms(train_dtm))

# Let's find high assicoation terms
findAssocs(train_dtm, terms = "cancer", corlimit = 0.4)
# breast has the only association. 2-gram required??

# # ngram
# library(RWeka)
# BigramTokenizer <- function(x) NGramTokenizer(x, Weka_control(min = 2, max = 2))
# TrigramTokenizer <- function(x) NGramTokenizer(x, Weka_control(min = 3, max = 3))
# FourgramTokenizer <- function(x) NGramTokenizer(x, Weka_control(min = 4, max = 3))
# 
# train_dtm_Bigram <- DocumentTermMatrix(train_txtCorp, 
#                                        control = list(tokenize = BigramTokenizer))
# findFreqTerms(train_dtm_Bigram)
# inspect(train_dtm_Bigram)
# train_dtm_Bigram$dimnames$Terms
# 
# train_dtm_Trigram <- DocumentTermMatrix(train_txtCorp, 
#                                         control = list(tokenize = TrigramTokenizer))
# findFreqTerms(train_dtm_Trigram)
# 
# train_dtm_Fourgram <- DocumentTermMatrix(train_txtCorp, 
#                                           control = list(tokenize = FourgramTokenizer))
# findFreqTerms(train_dtm_Fourgram)
# # ngrams can be rolled up using slam package.
# # ngrams not working because the data has been cleaned. Need to do it beforehand
# # giving this a pass
  
# Recreate dtm sides clip
# remove frequency
ndocs <- length(train_txtCorp)
# ignore overly sparse terms (appearing in less than 2% of the documents)
minDocFreq <- ndocs * 0.05
# ignore overly common terms (appearing in more than 95% of the documents)
maxDocFreq <- ndocs * 0.99
# train_dtm_Clip <- DocumentTermMatrix(train_txtCorp, 
#                          control = list(bounds = list(global = c(minDocFreq, maxDocFreq))))
# findFreqTerms(train_dtm_Clip, lowfreq = 100000)
# inspect(train_dtm_Clip) # Maximum term length is reduced drastically

# TFIDF custom function
TFIDF <- function(vector) {
  # tf 
  news_corpus  <- Corpus( VectorSource(vector) )
  control_list <- list(removePunctuation = TRUE, stopwords = TRUE, tolower = TRUE)
  tf <- TermDocumentMatrix(news_corpus, control = control_list) %>% as.matrix()
  
  # idf
  idf <- log( ncol(tf) / ( 1 + rowSums(tf != 0) ) ) %>% diag()
  return( crossprod(tf, idf) )
}

idfWeight <- function(x) weightTfIdf(x, normalize = TRUE)
# remove min and max freq words and apply tfidf
# train_dtm_tfidf <- DocumentTermMatrix(train_txtCorp, 
#                                 control = list(bounds = list(global = c(minDocFreq, maxDocFreq)), 
#                                                weighting = idfWeight))

train_dtm_tfidf <- DocumentTermMatrix(train_txtCorp, 
                                      control = list(weighting = idfWeight))

# -- save 3
# This makes a matrix that is p% empty space, maximum.i.e removes > p% empty
library(NLP)
library(tm)
train_dtm_tfidf_ns <- removeSparseTerms(train_dtm_tfidf, 0.80)

# cosine similarity
library(proxy)
library(dplyr)
# cosine_dist_mat <- as.matrix(dist(t(train_dtm_tfidf_ns_set), method = "cosine"))
# diag(cosine_dist_mat) <- NA
# cosine_dist <- apply(cosine_dist_mat, 2, mean, na.rm=TRUE)

Cosine <- function(x, y) {
  similarity <- sum(x * y) / ( sqrt( sum(y ^ 2) ) * sqrt( sum(x ^ 2) ) )
  # given the cosine value, use acos to convert back to degrees
  # acos returns the radian, multiply it by 180 and divide by pi to obtain degrees
  return( acos(similarity) * 180 / pi )
}

pr_DB$set_entry(FUN = Cosine, names = c("Cosines"))
cosDM <- as.matrix(dist(as.matrix(train_dtm_tfidf_ns), method = "Cosine"))
pr_DB$delete_entry("Cosines")

Hiercluster <- hclust(cosDM, method = "ward.D")
plot(Hiercluster)
rect.hclust(Hiercluster, 17)


# Apply LSI
# https://www.kaggle.com/danyalkhaliq/lsa-correlation-rf/code
library(SnowballC)
library(lsa)
myLSAspace <- lsa(train_dtm_tfidf_ns, dims=dimcalc_raw())
myLSAspaceTM <- as.matrix(train_dtm_tfidf_ns)
lsa::cosine(myLSAspaceTM[1,], myLSAspaceTM[10,])

# Done data cleaning
# Bind it to original dataframe
train <- cbind(train, as.matrix(train_dtm_tfidf_ns))
# Backup
#bkp_train <- train
#train <- bkp_train

# Need to use make.names as text is picking reserved keywords
names(train) <- make.names(names(train))

# Now since the dtm is consistent across test train, we split now
train$Class <- as.numeric(as.character(train$Class))
test <- train[train$Class == -1,]
train <- train[train$Class != -1,]

# Start modelling
yVar <- 'Class'
excludeList <- c('ID', 'Gene', 'Variation', 'txt')
includeList <- names(train)[!names(train) %in% c(excludeList,yVar)]

# Create stratified sampling due to class imbalance
library(caret)
set.seed(seedVal)
trainPct <- .8
testPct <- 1 - trainPct
inTrain <- createDataPartition(y = train[,c(yVar)], p = trainPct, list = FALSE)
traindata <- train[inTrain, ]
testdata <- train[-inTrain, ]
stopifnot(nrow(traindata) + nrow(testdata) == nrow(train))

library("xgboost")
# Class needs to be from 0 to 8
traindata[,yVar] <- as.integer(traindata[,yVar]) - 1
testdata[,yVar] <- as.integer(testdata[,yVar]) - 1

# Prepare matrix
mtrain <- model.matrix(~.+0,data = traindata[,includeList]) 
mtest <- model.matrix(~.+0,data = testdata[,includeList])
mtestNew <- model.matrix(~.+0,data = test[,includeList])

dtrain <- xgb.DMatrix(data = mtrain,label = traindata[,yVar])
dtest <- xgb.DMatrix(data = mtest,label=testdata[,yVar])
dtestNew <- xgb.DMatrix(data = mtestNew)


######################## Grid Tune #########################

library(mlr)
tm <- proc.time()

# Since we have best rounds using default params, we can use grid tune to improve performance
# convert characters to factors
traindata$Class <- as.factor(traindata$Class)
testdata$Class <- as.factor(testdata$Class)

# create tasks
traintask <- makeClassifTask(data = traindata[,c(includeList,yVar)],target = "Class")
testtask <- makeClassifTask(data = testdata[,c(includeList,yVar)],target = "Class")

# One hot encoding
traintask <- createDummyFeatures (obj = traintask)
testtask <- createDummyFeatures (obj = testtask)

#create learner
lrn <- makeLearner("classif.xgboost",predict.type = "prob")
lrn$par.vals <- list( objective="multi:softprob", eval_metric="mlogloss", nrounds=100L)

params <- makeParamSet( makeNumericParam("eta", lower = 0, upper = 0.3, trafo = function(x) x+0.1),
                        makeDiscreteParam("booster", values = c("gbtree","gblinear")),
                        makeIntegerParam("max_depth",lower = 3L, upper = 20L),
                        makeNumericParam("min_child_weight",lower = 1L, upper = 20L),
                        makeNumericParam("subsample", lower = 0.4, upper = 1, trafo = function(x) x+0.5),
                        makeNumericParam("colsample_bytree", lower = 0.4, upper = 1, trafo = function(x) x+0.5))


#set resampling strategy
rdesc <- makeResampleDesc("CV",stratify = T,iters=5L)

#search strategy
ctrl <- makeTuneControlRandom()

#set parallel backend # Seriously I have to work on others as well.
# library(parallel)
# library(parallelMap)
# parallelStartSocket(cpus = detectCores())

#parameter tuning
tuned <- tuneParams(learner = lrn, task = traintask,
                    resampling = rdesc, measures = acc,
                    par.set = params, control = ctrl, show.info = T)
print(tuned)
elapsedtm <- proc.time() - tm
print(elapsedtm[3])

######################## Grid Tune End #########################

#Tune result: 
#Op. pars: eta=0.177; booster=gbtree; max_depth=5; min_child_weight=17.5; 
#subsample=0.933; colsample_bytree=0.933
#acc.test.mean=0.631 at 70 best round 48

#Tune result:
#Op. pars: eta=0.177; booster=gbtree; max_depth=5; min_child_weight=17.5; 
#subsample=0.933; colsample_bytree=0.933
#acc.test.mean=0.635 at 95 best round 42

#Tune result: (Best performance till now)
#Op. pars: eta=0.177; booster=gbtree; max_depth=5; min_child_weight=17.5; 
#subsample=0.933; colsample_bytree=0.933
# at 80 best round 73

######################## Re run using tuned results ########

params <- list(booster = "gbtree", objective = "multi:softprob", 
               eta=0.177, gamma=0, max_depth=5, min_child_weight=17.5,
               seed = seedVal, nthread = 4, num_class = 9,
               subsample=0.933, colsample_bytree=0.933)

# Calculate the best round
xgbcv <- xgb.cv( params = params, data = dtrain, nrounds = 200, nfold = 5,
                  showsd = T, stratified = T, print_every_n = 10,
                  early_stopping_rounds = 20, maximize = F)

# Train model with the best round
xgbmodel <- xgb.train (params = params, data = dtrain, nrounds = 73, 
                       watchlist = list(val=dtest,train=dtrain), 
                       print_every_n = 10, early_stopping_rounds = 10, 
                       maximize = F , eval_metric = "mlogloss")

summary(xgbmodel)

# Compute feature importance matrix
impFeatures <- xgb.importance(includeList, model = xgbmodel)
# Nice graph
xgb.plot.importance(impFeatures[1:20,])

# Predict test instances
pred_prob <- predict(xgbmodel, dtest)
pred <- matrix(pred_prob, ncol=9, byrow=TRUE)
pred_labels <- max.col(pred) - 1

# Measure performance
confusion <- confusionMatrix(pred_labels, testdata[,yVar])
confusion

# Predict new instances
pred_prob_new <- predict(xgbmodel, dtestNew)
df_pred_labels_new <- data.frame(t(matrix(pred_prob_new, nrow=9, ncol=length(pred_prob_new)/9)))
names(df_pred_labels_new) <- c("class1","class2","class3","class4","class5","class6","class7","class8","class9")


################################ Create o/p File #################################

# Append values to test
submission <- test[,'ID']
submission <- data.frame(cbind(submission, df_pred_labels_new))
names(submission)[1] <- 'ID'

# Export to file
write.table(submission, file = "Output.csv", quote = FALSE, row.names=FALSE, sep=",")
