
library(readr)
library(plyr)
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
can_v_pred <- data.frame(canonical = numeric(),
                         predictions = numeric())

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
  
  all_subs <- length(levels(factor(mdf$siteid)))
  # filter KSEA
  mdf <- KSEAfilter(mdf, KSEAdb, reverse = T)
  predsubs <- length(levels(factor(mdf$siteid)))
  
  can_subs <- all_subs-predsubs
  tmpdf <- data.frame(can_subs,predsubs)
  
  row.names(tmpdf) <- names(fc_df)[i]
  
  
  can_v_pred <- rbind(can_v_pred, tmpdf)
  message(i)
}

can_v_pred_nodup <- can_v_pred[!duplicated(can_v_pred),]

can_v_pred_nodup$perc <- can_v_pred_nodup$can_subs/(can_v_pred_nodup$can_subs + can_v_pred_nodup$predsubs)

mean(can_v_pred_nodup$perc)
sd(can_v_pred_nodup$perc)

write.csv(can_v_pred_nodup, file = "Results/allconds_substratecounts.csv", row.names = T)

## combine all results
## define positive set- condition_kin, and negative set
## assign positive/negative assignments, run analysis (ROC, precision-recall)
