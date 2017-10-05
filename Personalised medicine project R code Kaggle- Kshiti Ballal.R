
# edWisor- Data Science Project 2
# Kaggle - Personalized Medicine: Redefining Cancer Treatment

# Objective: 
# To develop a Machine Learning algorithm that automatically classifies genetic variations.

# Once sequenced, a cancer tumor can have thousands of genetic mutations.
# But the challenge is distinguishing the mutations that contribute to 
# tumor growth (drivers) from the neutral mutations (passengers). 
 

# Submitted by Kshiti Ballal


rm(list=ls(all=TRUE))  # Clear all the objects in the console

setwd("C:/MyFiles/Kshiti/Edwisor/PROJECT 2")  # Set current working directory

# Install required pacakges

#-----------------------------------------------------------------------------------------------------------------

# Exploratory data Analysis


library(data.table)
library(readr)

# Load CSV files

#train_text = read_csv("C:/MyFiles/Kshiti/Edwisor/PROJECT 2/training_text")
#test_text = read_csv("C:/MyFiles/Kshiti/Edwisor/PROJECT 2/test_text") 

# The above two data frames need to be formated and cleaned

train_text = do.call(rbind,strsplit(readLines('C:/MyFiles/Kshiti/Edwisor/PROJECT 2/training_text'),'||',fixed=T))
test_text = do.call(rbind,strsplit(readLines('C:/MyFiles/Kshiti/Edwisor/PROJECT 2/test_text'),'||',fixed=T))

# Convert data into dataframe from matrix
train_text = as.data.table(train_text)
test_text = as.data.table(test_text)

#Remove the extra first row
train_text = train_text[-1,] 
test_text = test_text[-1,] 

#Rename the row headers
colnames(test_text) = c("ID", "Text")
colnames(train_text) = c("ID", "Text") 

train = fread("C:/MyFiles/Kshiti/Edwisor/PROJECT 2/training_variants", sep=",", stringsAsFactors = T)
test = fread("C:/MyFiles/Kshiti/Edwisor/PROJECT 2/test_variants", sep=",", stringsAsFactors = T)


# Structures of our data

library(dplyr)

glimpse(train)
#Observations: 3,321, Variables: 4, 
#ID and class are integers Gene and variatons are factors

prop.table(table(train$Class))
# Most observations belong to class 7

summary(train,maxsum = 8)
#Most frequent Gene- BRCA1: 264 times 
#Most frequent Variation- Truncating Mutations: 93 times           
   

# Type conversions before merge
train_text$ID = as.numeric(train_text$ID)
test_text$ID = as.numeric(test_text$ID)

# Joining test and train variants and text by ID
train = merge(train,train_text,by="ID")
test = merge(test,test_text,by="ID")

test$Class = -1  #Assign -1 to all test classes since class field is absent in test data
AllData= rbind(train,test) # Combine train and test data

rm(test_text,train_text,test) # Remove the unwanted data frames
gc()

#----------------------------------------------------------------------------------------------------------------

# # Missing value analyis
# 
# sum(is.na(AllData))
# #No missing values in the data
# 
# #Grubbs test for outliers
# 
# library(outliers)
# grubbs.test(AllData$Class, type = 10) 
# # No outliers present

#----------------------------------------------------------------------------------------------------

# Visualisation

library(ggplot2)

# Frequencies of each class

Class_freq = data.frame(table(train$Class))
colnames(Class_freq) = c('Class','Class Frequencies')
ggplot(data=Class_freq,aes(x=Class_freq[,1],y=Class_freq[,2]))+
  geom_bar(stat="identity",fill="antiquewhite",size=1, colour = "steelblue") + 
  ggtitle('Class frequencies')+ theme_bw()+
  xlab('Classes') + ylab('Frequencies')

# Most observations belong to class 7, 2nd most- class 4,3rd most- class 1


# Gene grouped by Class

gene_per_class = aggregate(train$Gene, by=list(train$Class),FUN= n_distinct)                                                        
colnames(gene_per_class) = c('Class','Genes')
ggplot(Class_freq, aes_string(x = Class_freq[,1], y = gene_per_class$Genes))+
  geom_bar(stat="identity",fill="springgreen1",size=1, colour = "indianred2") + 
  ggtitle('Genes per Class')+ theme_bw()+
  xlab('Classes') + ylab('Genes')

# Most genes belong to class 1


# Variations grouped by Class

var_per_class = aggregate(train$Variation, by=list(train$Class),FUN= n_distinct)                                                        
colnames(var_per_class) = c('Class','Variations')
ggplot(Class_freq, aes_string(x = Class_freq[,1], y = var_per_class$Variations))+
  geom_bar(stat="identity",fill="coral",size=1, colour = "yellowgreen") + 
  ggtitle(' Variations per Class')+ theme_bw()+
  xlab('Classes') + ylab('Variations')
 

# Most variatioons belong to class 7


# Merge the frequencies grouped by class
Freq_table = merge(var_per_class,gene_per_class,by="Class")
Freq_table = merge(Class_freq,Freq_table,by="Class")
rm(var_per_class,gene_per_class,Class_freq)


#Grouping AllData by gene

gene_count = AllData %>%
  group_by(Gene) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  filter(counts>50)

#Plotting the top 20 gene mutations
ggplot(gene_count, aes_string(x =  gene_count$Gene, y = gene_count$counts))+
  geom_bar(stat="identity",fill="steelblue",size=1, colour = "pink") + 
  ggtitle('Top 20 Gene mutations')+ theme_bw()+
  xlab('Genes') + ylab('Count')

#Most frequent Gene- BRCA1: 264 times 

#Grouping AllData by variation

var_count = AllData %>%
  group_by(Variation) %>%
  summarise(counts = n())%>%
  arrange(desc(counts))%>%
  filter(counts>2)

#Plotting the top 14 variations
ggplot(var_count, aes_string(x = var_count$Variation, y = var_count$counts))+
  geom_bar(stat="identity",fill="goldenrod1",size=1, colour = "maroon") + 
  ggtitle(' Top 14 Variations')+ theme_bw()+
  xlab('Variations') + ylab('Count')

#Most frequent Variation- Truncating Mutations: 93 times

#--------------------------------------------------------------------------------------------------------------

# Text Mining

#library(stringr)
#library(plyr)
#detach(package='ggplot2',unload=FALSE)
library(tm)

AllDataCorpus = Corpus(VectorSource(AllData$Text))

#case folding
AllDataCorpus = tm_map(AllDataCorpus, content_transformer(tolower))

#remove stop word
AllDataCorpus = tm_map(AllDataCorpus, removeWords, c('may','previous','also','us','use','used','using','well','yes','say','can','take','one','two','three','four',stopwords('english')))

#remove punctuation marks
AllDataCorpus = tm_map(AllDataCorpus, removePunctuation)

#remove numbers
AllDataCorpus = tm_map(AllDataCorpus, removeNumbers)

#remove unnecesary spaces
AllDataCorpus = tm_map(AllDataCorpus, stripWhitespace)

#stemming
AllDataCorpus = tm_map(AllDataCorpus, stemDocument, language="english")

tdm = TermDocumentMatrix(AllDataCorpus)

# Word frequency

library(slam)

words_freq = rollup(tdm, 2, na.rm=TRUE, FUN = sum)
words_freq = as.matrix(words_freq)
words_freq = data.frame(words_freq)

words_freq$words = row.names(words_freq)
row.names(words_freq) = NULL
names(words_freq) = c("Frequency","Words")
words_freq = words_freq[,c(2,1)]
words_freq= subset(words_freq,Frequency>35000)

head(arrange(words_freq,desc(Frequency)), n = 15)

# Word cloud for top 120 words
library(wordcloud) 
wordcloud(AllDataCorpus, max.words = 120, scale=c(5, .3), colors=brewer.pal(6, "Dark2"))

#-----------------------------------------------------------------------------------------------------------

# Building Frequency-Inverse Document Frequency from Corpus

#calculate the terms frequency

findFreqTerms(tdm, 100)
dtm = DocumentTermMatrix(AllDataCorpus, control = list(weighting = weightTfIdf))
dtm = removeSparseTerms(dtm, 0.95)
 
AllDataDTM = cbind(AllData, as.matrix(dtm))
 
# Sparse matrix

library(Matrix)

varnames = setdiff(colnames(AllDataDTM), c("ID", "Class", "Text"))
train_sparse = Matrix(as.matrix(sapply(AllDataDTM[Class > -1, varnames, with=FALSE],as.numeric)), sparse=TRUE)
test_sparse = Matrix(as.matrix(sapply(AllDataDTM[Class == -1, varnames, with=FALSE],as.numeric)), sparse=TRUE)
y_train = AllDataDTM[Class > -1,Class]-1
test_ids = AllDataDTM[Class == -1,ID]
dtrain = xgb.DMatrix(data=train_sparse, label=y_train)
dtest = xgb.DMatrix(data=test_sparse)


# Prediction

library(xgboost)
#library(caret)
#library(syuzhet) 

# Params for xgboost

param = list(booster = "gbtree",
              objective = "multi:softprob",
              eval_metric = "mlogloss",
              num_class = 9,
              eta = .2,
              gamma = 1,
              max_depth = 5,
              min_child_weight = 1,
              subsample = .7,
              colsample_bytree = .7)

#rounds <- which.min(xgb_cv$dt[, test.mlogloss.mean])

# Train the model
xgb_model = xgb.train(data = dtrain,
                       params = param,
                       watchlist = list(train = dtrain),
                       nrounds = 250,
                       verbose = 1,
                       print.every.n = 5)


# Feature importance matrix
names = dimnames(train_sparse)[[2]]
importance_matrix = xgb.importance(names,model=xgb_model)
head(importance_matrix)
xgb.plot.importance(importance_matrix[1:100,],10, xlab = "Relative importance")


# Final prediction
Drivers = as.data.table(t(matrix(predict(xgb_model, dtest), nrow=9, ncol=nrow(dtest))))
colnames(Drivers) = c("class1","class2","class3","class4","class5","class6","class7","class8","class9")

# Write the output in a csv file 
write.table(data.table(ID=test_ids, Drivers), "C:/MyFiles/Kshiti/Edwisor/PROJECT 2/submission.csv", sep=",", dec=".", quote=FALSE, row.names=FALSE)

submission <- read_csv("C:/MyFiles/Kshiti/Edwisor/PROJECT 2/submission.csv")
View(submission)