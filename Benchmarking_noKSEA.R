
library(readr)

# conditions_df <- read_excel("Data/benchmark/Supplementary_Table2.xlsx",
#                             sheet = "Conditions")
# write.csv(conditions_df, "Data/benchmark/Conditions.csv", row.names = F)
# sites_df <- read_excel("Data/benchmark/Supplementary_Table2.xlsx",
#                             sheet = "Phosphosites")
# write.csv(sites_df, "Data/benchmark/Phosphosites.csv", row.names = F)

conditions_df <- read.csv("Data/benchmark/Conditions.csv")
sites_df <- read.csv("Data/benchmark/Phosphosites.csv")
fc_df <- read.csv("Data/benchmark/Foldchanges.csv")

## code results such that each experiment output can be expressed as condition_kin : pKSEA
library(pKSEA)
sites_df$Phosphosites <- paste( sites_df$residues,sites_df$positions, sep = "")
names(sites_df)[names(sites_df)== "gene_name"] <- "GN"

#sites_df

##initialize results file
compiled_results <- data.frame(Kinase = levels(factor(NetworKINPred_db$kinase_id)))

for (i in 1:ncol(fc_df)){
  sumstat <- cbind.data.frame(Phosphosites = sites_df$Phosphosites, 
                              GN= sites_df$GN, 
                              t = fc_df[,i],
                              Peptide = row.names(sites_df)
                              )
  names(sumstat)[names(sumstat)== names(fc_df)[i]] <- "t"
  sumstat <- sumstat[complete.cases(sumstat),]
  
  sumstat$fc <- NA
  sumstat$pval <- NA
  
  mdf <- get_matched_data(sumstat, NetworKINPred_db)
  
  # filter KSEA
  mdf <- KSEAfilter(mdf, KSEAdb, reverse = F)
  
  mdf_calc <- calc_contribution(mdf)
  mdf_scores <- getscores(mdf_calc)
  mdf_perm <- permtest(mdf_calc, perms= 1000, seed=123)
  
  mdf_results <- perc.permutation(mdf_scores, mdf_perm)
  mdf_results$Kinase <- row.names(mdf_results)
  
  tmpdf <- data.frame(mdf_results[,2:3])
  names(tmpdf)[names(tmpdf)== "permutationScore"] <- names(fc_df)[i]
  
  compiled_results <- merge(compiled_results, tmpdf, by= "Kinase", all = T)
  message(paste(i,".."))
}

write.csv(compiled_results, file = "Results/180125_Benchmark_pred_only.csv", row.names = F)

## combine all results
## define positive set- condition_kin, and negative set
## assign positive/negative assignments, run analysis (ROC, precision-recall)
