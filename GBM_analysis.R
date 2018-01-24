##Reanalysis of gbm phospho

datadir <- "Z:/Documents/Dropbox/PhD/pKSEA_Chapter3/Data/GBM"

##import phosphodata and perform prep
library(readr)
library(readxl)
data <- read_excel(file.path(datadir,"PhosphoStats4L4S.xlsx"))

colnames(data)
data <- data[,c(4:11,13:20,25:27)]
names(data)[names(data)== "Modified Peptide Sequence"] <- "seq"

dataonly <- data[,1:16]
##imputing missing values (R3.4.1)
library(norm)
library(imputeLCMD)
# source("https://bioconductor.org/biocLite.R")
# biocLite("impute")
# dataonly <- as.matrix(mean_ints)
dataonly <- as.matrix(dataonly)
sum(is.na(dataonly))
rephist <- hist(log2(dataonly), breaks = 20)
rephist$breaks

min(dataonly)
table(dataonly==0)
rephist$counts
table(log2(dataonly) <13)

#replace all low scores with NA
dataonly[dataonly==min(dataonly[,1])] <- NA

mis_ind <- which(is.na(dataonly), arr.ind=T)

table(mis_ind[,2])

boxplot(log2(dataonly))

mistable <- table(mis_ind[,1])
mistable <- as.data.frame(mistable)
table(mistable$Freq)
rows_mis_nonrand <- as.numeric(levels(mistable$Var1)[mistable$Freq >=2])
rows_mis_rand <- as.numeric(levels(mistable$Var1)[mistable$Freq <2])

m.s <- model.Selector(dataonly)
m.s[[1]][rows_mis_nonrand] <- 0
table(m.s[[1]])

interestingrows <- dataonly[rows_mis_nonrand,]

dataonly_comp <- impute.MAR.MNAR(dataonly,m.s, method.MAR ="KNN", method.MNAR= "MinDet")

# #Using min.det
# dataonly_comp <-impute.MinDet(dataonly)
#
# #using min.prob
# set.seed(1234)
# dataonly_comp <- impute.MinProb(dataonly, tune.sigma = 0.02)
# dataonly_comp[dataonly_comp <0] <- NA
# dataonly_comp <-impute.MinDet(dataonly_comp)
#
# #Using KNN
# dataonly_comp <-impute.knn(dataonly)[[1]]
# sum(dataonly_comp < 0)
#
# interestingrows_imp <- dataonly_comp[rows_mis_nonrand,]

rephist <- hist(log10(dataonly_comp), breaks = 20, main="new")
summary(rephist)
rephist$breaks
rephist$counts

colnames(data)
datafinal <- cbind(dataonly_comp, data[,c(17:19)])

boxplot(log10(datafinal[,1:16]))

##Average technical replicates
dataonly <- datafinal[,1:16]
mean_ints <- data.frame(matrix(nrow = nrow(dataonly), ncol = ncol(dataonly)/2))

for(i in 1:4){
  mean_ints[,i] <- (dataonly[,i*2]+dataonly[,i*2-1])/2
  message(i*2, i*2-1)
}
for(i in 9:12){
  idx1 <-
  mean_ints[,i-4] <- dataonly[,i]+dataonly[i+4]
  message(i," ", i+4)
}

colnames(mean_ints)
colnames(mean_ints) <- c("STS1", "STS2", "STS3", "STS4", "LTS1", "LTS2", "LTS3", "LTS4")

means_final <- cbind(mean_ints,datafinal[17:19])

write.table(means_final, file = file.path(datadir, "GBM_mean_ints.txt"), 
            quote = T, row.names = F, sep= "\t")

##########################
# Summary Stat Calcs#####
#########################

data <- read.table(file.path(datadir, "GBM_mean_ints.txt"), header = T, sep= "\t",
                   stringsAsFactors = F)

colnames(data)
sts <- c(1:4)
lts <- c(5:8)

getcalcs <- function (data, colEx, colCon){
  
  calcs <- as.data.frame(matrix(nrow=nrow(data), ncol= 5,
                                dimnames = list(c(1:nrow(data)),c("fc", "t", "pval",
                                                                    "seq", "zscore"))))
  
  for (i in 1:nrow(data)){
    calcs$fc[i] <- rowMeans(data[i,colEx])/rowMeans(data[i, colCon])
    calcs$t[i] <- t.test(t(data[i,colEx]), t(data[i,colCon]))[["statistic"]]
    calcs$pval[i] <- t.test(t(data[i,colEx]), t(data[i,colCon]))[["p.value"]]
    
    calcs$seq[i] <- data$seq[i]
  }
  # calcs$pval <- p.adjust(calcs$pval, method = "fdr")
  
  meanlog2fc <- mean(log2(calcs$fc))
  calcs$zscore <- (log2(calcs$fc)-meanlog2fc)/sd(log2(calcs$fc))
  return(calcs)
}

calc_GBMsurv <- getcalcs(data, sts, lts)

write.csv(calc_GBMsurv, file.path(datadir, "GBM_sumstats.csv"), row.names = F)
###################
# Prep for pKSEA########
####################
data <- read.table(file.path(datadir, "GBM_mean_ints.txt"), header = T, sep= "\t")
sumstats <- read.csv(file.path(datadir, "GBM_sumstats.csv"), header = T)

#extract descriptions
library(stringr)
getGN <- function(fulldesc, initial.pattern, final){
  pattern <- paste(initial.pattern, ".*?", final, sep="")
  gene <- str_extract(fulldesc, pattern)
  gene <- sub(initial.pattern, "", gene)
  gene <- sub(final, "", gene)
  return(gene)
}
# str_extract(data$Protein.Description[1], "GN=.*? ")
# getGN(data$Protein.Description[1], "GN=", " ")

data$GN <- getGN(data$`Protein.Description`, "GN=", " ")

# Annotate with official names?..
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# biocLite("Biostrings")
# library("biomartr")
# library("biomaRt")
# peptides <- read.delim("C:/Users/pll21/Dropbox/PhD/pKSEA_Chapter3/Data/GBM/peptides.tdv")
# 
# listMarts()
# ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
# 
# View(listDatasets(ensembl))

##Annotate with PX info
peptides <- read.delim("C:/Users/pll21/Dropbox/PhD/pKSEA_Chapter3/Data/GBM/peptides.tdv",
                       stringsAsFactors = F)

table(data$seq %in% peptides$peptide)
data_sub <- data[data$seq %in% peptides$peptide,]
data_sub_merge <- merge(x = data, y= peptides, by.x = "seq", by.y = "peptide")
colnames(data_sub_merge)
seq_GN <- data_sub_merge[,c(1,12:14)]

sumstat_merge <- merge(x= sumstats, y= seq_GN, by.x = "seq", by.y= "seq")

sumstat_merge$Phosphosites <- NA
# combine residue lists
df <- sumstat_merge

combineRes <- function(df){
  df$Phosphosites <- NA
  df$Known.Sites[df$Known.Sites== ""] <- NA
  df$Unknown.Sites[df$Unknown.Sites== ""] <- NA
  
  for (i in 1:nrow(df)){
    if (!is.na(df$Known.Sites[i]) & !is.na(df$Unknown.Sites[i])) {
      df$Phosphosites[i] <- paste(df$Known.Sites[i], df$Unknown.Sites[i], sep=",")
    } else if (is.na(df$Unknown.Sites[i])){
      df$Phosphosites[i] <- df$Known.Sites[i]
    } else if (is.na(df$Known.Sites[i])){
      df$Phosphosites[i] <- df$Unknown.Sites[i]
    }
  }
  return(df)
}
sumstat_merge_res <- combineRes(sumstat_merge)

colnames(sumstat_merge_res)[1] <- "Peptide"

################################
# Run pKSEA #################
################################
library(pKSEA)

# sumstat_merge_res[duplicated(sumstat_merge_res$Peptide) |
#                     duplicated(sumstat_merge_res$Peptide, fromLast=T), ]
# sumstat_merge_res <- sumstat_merge_res[!duplicated(sumstat_merge_res$Peptide),]
# 
# write.csv(sumstat_merge_res, file.path(datadir, "sumstat_GBMsurv_full.csv"),
#           row.names = F)

sumstat_merge_res <- read.csv(file.path(datadir, "sumstat_GBMsurv_full.csv"), header = T)
table(duplicated(sumstat_merge_res))


db_matched <- get_matched_data(sumstat_merge_res, NetworKINPred_db)
GBM_comp <- compare(matched_data = db_matched, 
                    predictionDB = NetworKINPred_db, kseadb = KSEAdb)
results_write(GBM_comp, outputpath= datadir ,outputname = "demo.csv", singlefolder= "GBMout")
