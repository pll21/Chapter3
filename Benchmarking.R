library(readxl)
library(readr)

conditions_df <- read_excel("Data/benchmark/Supplementary_Table2.xlsx",
                            sheet = "Conditions")
sites_df <- read_excel("Data/benchmark/Supplementary_Table2.xlsx",
                            sheet = "Phosphosites")
fc_df <- read_csv("Data/benchmark/Foldchanges.csv")
names(fc_df) <- paste("e", names(fc_df), sep = "")

head(fc_df)
## code results such that each experiment output can be expressed as condition_kin : pKSEA
library(pKSEA)
sites_df$Phosphosites <- paste( sites_df$residues,sites_df$positions, sep = "")
names(sites_df)[names(sites_df)== "gene_name"] <- "GN"

sites_df

##initialize results file
compiled_results <- data.frame(Kinase = levels(factor(NetworKINPred_db$kinase_id)))

for (i in 1:ncol(fc_df)){
  sumstat <- cbind.data.frame(Phosphosites = sites_df$Phosphosites, 
                              GN= sites_df$GN, 
                              t = fc_df[,i],
                              Peptide = row.names(sites_df)
                              )
  names(sumstat)[names(sumstat)== names(fc_df[,i])] <- "t"
  sumstat <- sumstat[complete.cases(sumstat),]
  
  sumstat$fc <- NA
  sumstat$pval <- NA
  
  mdf <- get_matched_data(sumstat, NetworKINPred_db)
  mdf_calc <- calc_contribution(mdf)
  mdf_scores <- getscores(mdf_calc)
  mdf_perm <- permtest(mdf_calc, perms= 1000, seed=123)
  
  mdf_results <- perc.permutation(mdf_scores, mdf_perm)
  mdf_results$Kinase <- row.names(mdf_results)
  
  tmpdf <- data.frame(mdf_results[,2:3])
  names(tmpdf)[names(tmpdf)== "permutationScore"] <- names(fc_df[,i])
  
  compiled_results <- merge(compiled_results, tmpdf, by= "Kinase", all = T)
}

## combine all results
## define positive set- condition_kin, and negative set
## assign positive/negative assignments, run analysis (ROC, precision-recall)