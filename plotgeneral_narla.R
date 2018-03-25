rm(list=ls())
## Permutation Heatmap

#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}


#########################################################
### B) Reading in data and transform it into matrix format
#########################################################
# install.packages("reshape")
library(reshape)

datapath <- "Results/170119test/"
colguide <- c("character", "NULL", "numeric")

inputFiles <- list.files(datapath, pattern = "full")
inputFiles
inputFiles <- inputFiles[c(1,
                           4,
                           2,
                           3,
                           5)]
inputFiles
groupNames <- c("AZD",
                  "\nCombo/DT\n(AZD contribution)",
                  "Combination",
                  "\nCombo/AZD\n(DT contribution)",
                  "DT")

data <- list()

for (i in 1:length(inputFiles)){
  data[[i]] <- read.csv(file.path(datapath,inputFiles[i]), header=T, colClasses = colguide)
  colnames(data[[i]])[2] <- groupNames[i]
  # data[[i]]$genes <- rownames(data[[i]])
}

merge_data <- merge(data[[1]], data[[2]], by='X')

for (i in 3:length(inputFiles)){
  merge_data <- merge(merge_data, data[[i]], by='X')
}

names(merge_data)
# merge_data <- subset(merge_data, Das < 100| `Combo/Rap\n(Das)` < 100|Combo < 100 |`Combo/Das\n(Rap)` < 100 | Rap <100)

#Conversion to correct data format (genes as row names, matrix numeric)
Kinase <- merge_data$X
merge_data <- data.matrix(merge_data[,-1])
rownames(merge_data)<- Kinase

# merge_data<- t(merge_data)

#Filtering
merge_Filter <- merge_data
merge_Filter <- merge_Filter[rowSums((merge_Filter > 95 | merge_Filter < 5))>=3, ]
# merge_Filter <- merge_Filter[,colSums((merge_Filter > 95 | merge_Filter < 5))>1]

#########################################################
### C) Customizing and plotting the heat map
#########################################################
heatmapfolder <- file.path(datapath,"heatmap")
ifelse(!dir.exists(heatmapfolder), dir.create(heatmapfolder), F)

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("royalblue3", "white", "tomato3"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,20,length=100),  # for blue
               seq(21,80,length=100),           # for white
               seq(81,100,length=100))             # for red

# creates a 5 x 5 inch image
pdf(file.path(heatmapfolder, "pKSEA_narlacross.pdf"),    # create PNG for the heat map
    width = 10,        # 5 x 300 pixels
    height = 12,          # 300 pixels per inch
    pointsize = 9)        # smaller font size

heatmap.2(merge_data,
          #           cellnote = TC_Merge,  # same data set for cell labels
          main = "pKSEA Permutation Percentile", # heat map title
          cexCol = 2,
          cexRow = 0.5,
          keysize= 0.3,
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(15,10),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",
          Rowv=T, # F for transposed
          reorderfun = function(d, w) reorder(d, -w),
          Colv=F)            # turn off column clustering T for transposed

dev.off()               # close the PNG device


##PLOT FILTERED

my_palette <- colorRampPalette(c("royalblue3", "white", "tomato3"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,20,length=100),  # for blue
               seq(21,80,length=100),           # for white
               seq(81,100,length=100))             # for red

# creates a 5 x 5 inch image
pdf(file.path(heatmapfolder, "pKSEA_narlacross_filtered.pdf"),    # create PNG for the heat map
    width = 10,        # 5 x 300 pixels
    height = 12,          # 300 pixels per inch
    pointsize = 9)        # smaller font size


hm_statfilter <- heatmap.2(merge_Filter,
          #           cellnote = TC_Merge,  # same data set for cell labels
          main = "pKSEA Permutation Percentile (Filtered)", # heat map title
          cexCol = 2,
          cexRow = 1.5,
          notecex = 1.3,
          keysize= 0.30,
          notecol="black",      # change font color of cell labels to black
          cellnote= merge_Filter,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(15,10),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",
          Rowv=T,
          reorderfun = function(d, w) reorder(d, w, anglo.FUN = mean),
          Colv=F)            # turn off column clustering

dev.off()               # close the PNG device

#Dendrogram
library(dendextend)
dend <- hm_statfilter$rowDendrogram
cut_class <- cutree(dend, k=5)
dend <- color_branches(dend, k=5)
plot(dend)
dend2 <- sort(dend)
plot(dend2)

kinasegroups <- list()
for (i in 1:5){
  kinasegroups[[i]] <- names(cut_class[cut_class==i])
}

kinasegroups

cluster <-kmeans(merge_Filter, centers = 5)
cluster$cluster
kinasegroups_k <- list()
for (i in 1:5){
  kinasegroups_k[[i]] <- names(cut_class[cut_class==i])
}

#get significant kinase lists
kinase_id2ens <- read.csv("C:/Users/pll21/Dropbox/PhD/SystemPhos2/kinase_ensembl.csv")

exp_frame <- as.data.frame(merge_data)
com_up <- data.frame(kin = row.names(exp_frame[exp_frame$Combo >= 95,]),stringsAsFactors = F)
com_up$ens <- kinase_id2ens$ens[match(com_up$kin, kinase_id2ens$id)]
das_up <- data.frame(kin = row.names(exp_frame[exp_frame$Das >= 95,]),stringsAsFactors = F)
das_up$ens <- kinase_id2ens$ens[match(das_up$kin, kinase_id2ens$id)]
rap_up <- data.frame(kin = row.names(exp_frame[exp_frame$Rap >= 95,]),stringsAsFactors = F)
rap_up$ens <- kinase_id2ens$ens[match(rap_up$kin, kinase_id2ens$id)]

com_down <- data.frame(kin = row.names(exp_frame[exp_frame$Combo <=5,]),stringsAsFactors = F)
com_down$ens <- kinase_id2ens$ens[match(com_down$kin, kinase_id2ens$id)]
das_down <- data.frame(kin = row.names(exp_frame[exp_frame$Das <=5,]),stringsAsFactors = F)
das_down$ens <- kinase_id2ens$ens[match(das_down$kin, kinase_id2ens$id)]
rap_down <-data.frame(kin =  row.names(exp_frame[exp_frame$Rap <=5,]),stringsAsFactors = F)
rap_down$ens <- kinase_id2ens$ens[match(rap_down$kin, kinase_id2ens$id)]

dascom <- intersect(das_down$kin, com_down$kin)
nrow(das_down)
nrow(com_down)
nrow(rap_down)

length(intersect(das_down$kin, com_down$kin))
length(intersect(rap_down$kin, com_down$kin))
length(intersect(das_down$kin, rap_down$kin))
length(intersect(dascom, rap_down$kin))

ifelse(!dir.exists("kinaselists"), dir.create("kinaselists"), F)
write.table(com_up$ens, file = file.path("kinaselists", "comUens.txt"), quote = F, sep = "\t", row.names = F, col.names = F)
write.table(das_up$ens, file = file.path("kinaselists","dasUens.txt"), quote = F, sep = "\t", row.names = F, col.names = F)
write.table(rap_up$ens, file = file.path("kinaselists","rapUens.txt"), quote = F, sep = "\t", row.names = F, col.names = F)

write.table(com_down$ens, file = file.path("kinaselists","comDens.txt"), quote = F, sep = "\t", row.names = F, col.names = F)
write.table(das_down$ens, file = file.path("kinaselists","dasDens.txt"), quote = F, sep = "\t", row.names = F, col.names = F)
write.table(rap_down$ens, file = file.path("kinaselists","rapDens.txt"), quote = F, sep = "\t", row.names = F, col.names = F)
