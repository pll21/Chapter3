library(dplyr)

kinconpairs_df <- read.csv("Data/benchmark/kinconpairs.csv", stringsAsFactors = F)
pKSEAresults <- read.csv("Results/180125_Benchmark_all.csv", stringsAsFactors = F)

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
i <- 50
negset <- data.frame(Condition = character(),
                     Kinase = character(),
                     Regulation = character(),
                     stringsAsFactors = F)

for(i in 1:nrow(counts)){
  npos <- counts$no_rows[i]
  
  randomneg <- sample(row.names(kincononly_df[kincononly_df$Condition != counts$Condition[i],]), npos, replace = F)
  neg <- kincononly_df[row.names(kincononly_df) %in% randomneg,]
  neg$Condition <- kincononly_df$Condition[i]
  negset <- rbind(negset, neg)
}

#recode kinases
levels(factor(kincononly_df$Kinase))


#need recode- akt1 to PKBalpha, aurka/b to AuroraA/B, BRAF to RAF1, CAMk2a to CaMKIIalpha, LCK to Lck, PDGFRA/B 
  #to PDGFRalpha/beta, PRKACA to PKAalpha, PRKCA to PKCalpha, RPS6KB1 to p70S6Kb, use p70S6K as surrogate for MTOR,
  #ABL to Abl


dict <- list(AKT1 = "PKBalpha", AURKA = "AuroraA", AURKB = "AuroraB", BRAF = "RAF1", CAMK2A = "CaMKIIalpha", LCK= "Lck",
             PDGFRA= 'PDGFRalpha', PDGFRB = 'PDGFRbeta', PRKACA = "PKCalpha", PRKCA = "PKCalpha", RPS6KB1 = "p70S6Kb",
             MTOR = "p70S6K", ABL= "Abl")

remap_pos <- kincononly_df[!(kincononly_df$Kinase %in% removelist),]
remap_neg <- negset[!(negset$Kinase %in% removelist),]
for(i in 1:length(dict)){
  remap_pos <- replace(remap_pos, remap_pos == names(dict[i]), dict[i])
  remap_neg <- replace(remap_neg, remap_neg == names(dict[i]), dict[i])
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

function <
for(i in 1:nrow(remap_pos)){
  usecol <- paste("e",remap_pos$Condition[i], sep= "")
  if(!is.na(pKSEAresults[,usecol][pKSEAresults$Kinase == remap_pos$Kinase[i]])){
    if(remap_pos$Regulation[i] == "up"){
      if(pKSEAresults[,usecol][pKSEAresults$Kinase == remap_pos$Kinase[i]] >= upthresh){
        scoring_df$TP = scoring_df$TP + 1
      } else {
        scoring_df$FN = scoring_df$FN + 1
      }
    } else if(remap_pos$Regulation[i] == "down"){
      if(pKSEAresults[,usecol][pKSEAresults$Kinase == remap_pos$Kinase[i]] <= threshold){
        scoring_df$TP = scoring_df$TP + 1
      } else {
        scoring_df$FN = scoring_df$FN + 1
      }
    }
  }
  
  if(!is.na(pKSEAresults[,usecol][pKSEAresults$Kinase == remap_neg$Kinase[i]])){
    if(remap_neg$Regulation[i] == "up"){
      if(pKSEAresults[,usecol][pKSEAresults$Kinase == remap_neg$Kinase[i]] >= upthresh){
        scoring_df$FP = scoring_df$TP + 1
      } else {
        scoring_df$TN = scoring_df$FN + 1
      }
    } else if(remap_pos$Regulation[i] == "down"){
      if(pKSEAresults[,usecol][pKSEAresults$Kinase == remap_neg$Kinase[i]] <= threshold){
        scoring_df$FP = scoring_df$TP + 1
      } else {
        scoring_df$TN = scoring_df$FN + 1
      }
    }
  }
}
