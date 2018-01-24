library(readxl)

conditions_df <- read_excel("Data/benchmark/Supplementary_Table2.xlsx",
                            sheet = "Conditions")
sites_df <- read_excel("Data/benchmark/Supplementary_Table2.xlsx",
                            sheet = "Phosphosites")
fc_df <- read_excel("Data/benchmark/Supplementary_Table2.xlsx",
                       sheet = "Foldchanges")

## code results such that each experiment output can be expressed as condition_kin : pKSEA
## combine all results
## define positive set- condition_kin, and negative set
## assign positive/negative assignments, run analysis (ROC, precision-recall)