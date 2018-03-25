library(dplyr)
library(plyr)

load("benchmarkingpermscores.RData")
load("kseaonly_permscores.RData")
kinconpairs_df <- read.csv("Data/benchmark/kinconpairs.csv", stringsAsFactors = F)

# pKSEAresults <- read.csv("Results/180125_Benchmark_all.csv", stringsAsFactors = F)
# pKSEAresults <- read.csv("Results/180125_Benchmark_pred_only.csv", stringsAsFactors = F)
# pKSEAresults <- read.csv("Results/180125_Benchmark_KSEA_only.csv", stringsAsFactors = F)

kincononly_df <- kinconpairs_df[,c(1,6,7)]

#missing in pKSEA- "ALK", "FLT4", "LRRK2", "MTOR", "PIK3CA"; recode MTOR for p70S6K
removelist <- c("ALK", "FLT4", "LRRK2", "PIK3CA")
kincononly_df <- kincononly_df[!(kincononly_df$Kinase %in% removelist),]

# Generate negative set
# factor conditions and count number of kinase-condition pairs for each
factor(kinconpairs_df$Condition)

counts <- kincononly_df %>%
  group_by(Condition) %>%
  summarise(no_rows= length(Condition))

#get unique conditions
kc_uniq <- kincononly_df[,2:3]
kc_uniq <- unique(kc_uniq)

#check condition counts
i <- 48

getnegset <- function(poscounts, kincononly_df, pos_pairs, kc_uniq){
  negset <- data.frame(Condition = character(),
                       Kinase = character(),
                       Regulation = character(),
                       stringsAsFactors = F)
  
  for(i in 1:nrow(poscounts)){
    npos <- poscounts$no_rows[i]
    pos_pairs <- kincononly_df[kincononly_df$Condition == poscounts$Condition[i],]
    kc_other <- kc_uniq[!((kc_uniq$Kinase %in% pos_pairs$Kinase) & (kc_uniq$Regulation %in% pos_pairs$Regulation)),]
    
    randomneg <- sample(row.names(kc_other), npos, replace = F)
    neg <- kincononly_df[row.names(kincononly_df) %in% randomneg,]
    neg$Condition <- poscounts$Condition[i]
    negset <- rbind(negset, neg)
  }
  return(negset)
}

#recode kinases
levels(factor(kincononly_df$Kinase))


#need recode- akt1 to PKBalpha, aurka/b to AuroraA/B, BRAF to RAF1, CAMk2a to CaMKIIalpha, LCK to Lck, PDGFRA/B 
  #to PDGFRalpha/beta, PRKACA to PKAalpha, PRKCA to PKCalpha, RPS6KB1 to p70S6Kb, use p70S6K as surrogate for MTOR,
  #ABL to Abl


dict <- list(AKT1 = "PKBalpha", AURKA = "AuroraA", AURKB = "AuroraB", BRAF = "RAF1", CAMK2A = "CaMKIIalpha", LCK= "Lck",
             PDGFRA= 'PDGFRalpha', PDGFRB = 'PDGFRbeta', PRKACA = "PKCalpha", PRKCA = "PKCalpha", RPS6KB1 = "p70S6Kb",
             MTOR = "p70S6K", ABL= "Abl")

remap<- function(kinconset, dict, removelist){
  kinconset <- kinconset[!(kinconset$Kinase %in% removelist),]
  for(i in 1:length(dict)){
    kinconset <- replace(kinconset, kinconset == names(dict[i]), as.character(dict[i]))
  }
  return(kinconset)
}

#evaluate against data
i <- 1
scoring_df <- data.frame(TP = 0,
                         TN = 0,
                         FP = 0,
                         FN = 0)
#find relevant column
threshold <- 5
upthresh <- 100-5

score_interations <- function(remap_pos, remap_neg, pKSEAresults, threshold){
  scoring_df <- data.frame(TP = 0,
                           TN = 0,
                           FP = 0,
                           FN = 0)
  posreadout_df <- data.frame(Condition= character(),
                           Kinase = character(),
                           Regulation = character(),
                           pKSEA = numeric(),
                           Hit = character())
  negreadout_df <- data.frame(Condition= character(),
                              Kinase = character(),
                              Regulation = character(),
                              pKSEA = numeric(),
                              Hit = character())
  upthresh <- 100-threshold
  for(i in 1:nrow(remap_pos)){
    bool <- T
    usecol <- paste("e",remap_pos$Condition[i], sep= "")
    if(!is.na(pKSEAresults[,usecol][pKSEAresults$Kinase == remap_pos$Kinase[i]])){
      if(remap_pos$Regulation[i] == "up"){
        if(pKSEAresults[,usecol][pKSEAresults$Kinase == remap_pos$Kinase[i]] >= upthresh){
          scoring_df$TP = scoring_df$TP + 1
        } else {
          scoring_df$FN = scoring_df$FN + 1
          bool <- F
        }
      } else if(remap_pos$Regulation[i] == "down"){
        if(pKSEAresults[,usecol][pKSEAresults$Kinase == remap_pos$Kinase[i]] < threshold){
          scoring_df$TP = scoring_df$TP + 1
        } else {
          scoring_df$FN = scoring_df$FN + 1
          bool <- F
        }
      }
    } else bool <- NA
    # message(paste("Condition: ", remap_pos$Condition[i],
    #               "\nKin: ", remap_pos$Kinase[i], ", Reg:", remap_pos$Regulation[i],
    #               "\npKSEA: ", pKSEAresults[,usecol][pKSEAresults$Kinase == remap_pos$Kinase[i]], sep = ""))
  posreadout_line <- data.frame(Condition = remap_pos$Condition[i],
                             Kinase = remap_pos$Kinase[i],
                             Regulation = remap_pos$Regulation[i],
                             pKSEA = pKSEAresults[,usecol][pKSEAresults$Kinase == remap_pos$Kinase[i]],
                             Hit = bool)
  posreadout_df <- rbind(posreadout_df, posreadout_line)
    
  }
  
  for(i in 1:nrow(remap_neg)){
    bool <- T
    usecol <- paste("e",remap_neg$Condition[i], sep= "")
    if(!is.na(pKSEAresults[,usecol][pKSEAresults$Kinase == remap_neg$Kinase[i]])){
      if(remap_neg$Regulation[i] == "up"){
        if(pKSEAresults[,usecol][pKSEAresults$Kinase == remap_neg$Kinase[i]] >= upthresh){
          scoring_df$FP = scoring_df$FP + 1
        } else {
          scoring_df$TN = scoring_df$TN + 1
          bool <- F
        }
      } else if(remap_neg$Regulation[i] == "down"){
        if(pKSEAresults[,usecol][pKSEAresults$Kinase == remap_neg$Kinase[i]] < threshold){
          scoring_df$FP = scoring_df$FP + 1
        } else {
          scoring_df$TN = scoring_df$TN + 1
          bool <- F
        }
      }
    } else bool <- NA
    negreadout_line <- data.frame(Condition = remap_neg$Condition[i],
                                  Kinase = remap_neg$Kinase[i],
                                  Regulation = remap_neg$Regulation[i],
                                  pKSEA = pKSEAresults[,usecol][pKSEAresults$Kinase == remap_neg$Kinase[i]],
                                  Hit = bool)
    negreadout_df <- rbind(negreadout_df, negreadout_line) 
  }

  output <- list(scoring_df = scoring_df, posreadout= posreadout_df, negreadout= negreadout_df)
  # return(output)
  return(scoring_df)
}


remap_pos <- remap(kincononly_df, dict= dict, removelist = removelist)


permscores <- list()

for(i in 1:60){
  negset <- getnegset(poscounts = counts, kincononly_df = kincononly_df, pos_pairs = pos_pairs, kc_uniq = kc_uniq)
  remap_neg <- remap(negset, dict = dict, removelist = removelist)
  
  score_df <- data.frame(TP = numeric(),
                         TN = numeric(),
                         FP = numeric(),
                         FN = numeric())
  for(j in -1:101){
    score <- score_interations(remap_pos = remap_pos, remap_neg = remap_neg, pKSEAresults = pKSEAresults, threshold = j)
    score_df <- rbind(score_df, score)
  }
  permscores[[i]] <- score_df
  message(paste(i, "..\n", sep = ""))
}

#take average of lists?
calcstats <- function(df){
  df$TPR <- df$TP/(df$TP+df$FN)
  df$FPR <- 1-df$TN/(df$TN+df$FP)
  df$pres <- df$TP/(df$TP+df$FP)
  return(df)
}

permscores_calc <- lapply(permscores, calcstats)

resultmeans <- aaply(laply(permscores_calc, as.matrix), c(2,3), mean)
resultsd <- aaply(laply(permscores_calc, as.matrix), c(2,3), sd)
resultmedian <- aaply(laply(permscores_calc, as.matrix), c(2,3), median)

resultmeans <- as.data.frame(resultmeans)
resultmedian<- as.data.frame(resultmedian)

plot(resultmedian$TPR~resultmedian$FPR,type="S", ylim = c(0,1),
     main = "pKSEA ROC Curve", xlab= "False Positive Rate", ylab= "True Positive Rate")
abline(0,1)
auc <- sum(resultmedian$TPR*diff(c(0,resultmedian$FPR)))

plot(resultmedian$pres~resultmedian$TPR, type="s", ylim=c(0,1), xlim=c(min(resultmedian$TPR[!is.na(resultmedian$pres)], na.rm = T), 1),
     main = "pKSEA Precision-recall Curve", xlab= "Recall", ylab= "Precision")
abline(v=0.5)

View(scores$negreadout)
remap_neg[order(remap_neg$Condition),]
