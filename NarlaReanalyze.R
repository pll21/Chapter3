library(pKSEA)
library(plyr)

datadir <- "C:/Users/pll21/Dropbox/PhD/SystemPhos2/Data/H358Comb/witht"

# datalist <- list.files(path = datadir, pattern = ".txt")
# 
# i <- 1
# for(i in 1:length(datalist)){
#   df <- read.delim(file.path(datadir, datalist[i]), stringsAsFactors = F)
#   colnames(df)[colnames(df)=="sites"] <- "Phosphosites"
#   df$Peptide <- rownames(df)
#   
#   write.csv(df, file = file.path(datadir, gsub(".txt", ".csv",datalist[i])), row.names = F)
# }

###
#Run pKSEA
###

batchrun(summaryfiledir = datadir, commonfilestring = ".csv", predictionDB = NetworKINPred_db, results_folder = "170119test",
         kseadb =KSEAdb)

##
#Get Substrate info and permutation test correlations
##

datalist <- list.files(path = datadir, pattern = ".csv")
datalist
df <- read.csv(file.path(datadir, datalist[2]))
mdf <- get_matched_data(df, predictionDB = NetworKINPred_db)
sublist <- getsubs(mdf)
sublist_summary<- count(df = sublist, vars = "kinase_id")
write.csv(sublist_summary, file= "Results/170119test/Narlasublist_summary.csv")

mdf_kseaonly <- KSEAfilter(mdf, KSEAdb)
mdf_ko_calc <- calc_contribution(mdf_kseaonly)

sublist <- getsubs(mdf_ko_calc)


mdf_noksea <- KSEAfilter(mdf, KSEAdb, reverse = T)
mdf_nk_calc <- calc_contribution(mdf_noksea)


