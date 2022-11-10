library(edgeR)
library(ggplot2)
library(pheatmap)
library(matrixStats)
library(forcats)
library(RColorBrewer)
library(stats)
library(biomaRt)
library(DESeq2)
library(dplyr)
library(viridis)
library(reshape2)
library(tidyr)
library(tibble)
library(cluster)
library(factoextra)
library(ggVennDiagram)


setwd("C:/Users/am331/OneDrive/Documents/UCLA/Signaling_Systems_Lab/RNASeq_2022")

bmdm_clean <- read.table(file = "bmdm_clean.txt")

# Add 1.5 to column names
#(bmdm_clean)[4] <- "BGlu_1.5"
#colnames(bmdm_clean)[8] <- "IFNb_1.5"
#colnames(bmdm_clean)[12] <- "IFNg_1.5"
#colnames(bmdm_clean)[16] <- "LPS_1.5"
#colnames(bmdm_clean)[20] <- "PBS_1.5"
#write.table(bmdm_clean, 
#            file = "bmdm_clean.txt")




#rpkm Normalization

bmdm_rpkm <- as.data.frame(rpkm(bmdm_clean[,3:22], gene.length = bmdm_clean$Length))
#rpkm_summary<-as.data.frame(apply(bmdm_rpkm,2,summary))
summary(bmdm_rpkm)

#histogram for RPKM Data
#**melted dataframes have no rownames/gene symbols!!
bmdm_rpkm_melted <- melt(bmdm_rpkm, variable.name = "StimulusCondition", value.name = "GeneExpression", row.names = TRUE)

ggplot(bmdm_rpkm_melted, aes(x = log2(GeneExpression)))+
  geom_histogram() + 
  labs(x = "Gene Expression (log2(RPKM))", y = "Counts", title = "Gene Expression Distribution, RPKM values") +
  scale_x_continuous(breaks = seq(-8,15, 1), lim = c(-8, 15))
  #xlim(0,100) +ylim(0,45000)


#TPM Normalization

TPM_rpk <- as.data.frame(bmdm_clean[,3:22]/bmdm_clean[,2])
bmdm_tpm <- as.data.frame(t( t(TPM_rpk) * 1000000 / colSums(TPM_rpk) ))
summary(bmdm_tpm)

test <- as.data.frame(colSums(bmdm_tpm))
#histogram for TPM Data
#**melted dataframes have no rownames/gene symbols!!
bmdm_tpm_melted <- melt(bmdm_tpm, variable.name = "StimulusCondition", value.name = "GeneExpression", row.names = TRUE)

ggplot(bmdm_tpm_melted, aes(x = log2(GeneExpression)))+
  geom_histogram() + 
  labs(x = "Gene Expression (log2(TPM))", y = "Counts", title = "Gene Expression Distribution, TPM values") +
  scale_x_continuous(breaks = seq(-8,15, 1), lim = c(-8, 15))

  # xlim(0,150) +ylim(0,35000)


#add gene symbols to rownames
rownames(bmdm_rpkm) <- bmdm_clean$Symbol
rownames(bmdm_tpm) <- bmdm_clean$Symbol



#Add pseudocount to rpkm and TPM normalized datasets and take log2fc
rpkm_pseud <- bmdm_rpkm+0.1

tpm_pseud <- bmdm_tpm+0.25


#Create 4 dataframes for rpkm (PBS+BGlu, PBS+IFNb, PBS+IFNg, PBS+LPS) 
#with genes with at least one log2FC value >0.7 (Induced!)

#BGlu + PBS log2FC>0.7
rpkm_BGlu_FC_PBS_raw <- data.frame(row.names = rownames(rpkm_pseud),
                            PBS_1.5FC0 = log2(rpkm_pseud$PBS_1.5/rpkm_pseud$PBS_0),
                            PBS_3FC0 = log2(rpkm_pseud$PBS_3/rpkm_pseud$PBS_0),
                            PBS_8FC0 = log2(rpkm_pseud$PBS_8/rpkm_pseud$PBS_0),
                            BGlu_1.5FC0 = log2(rpkm_pseud$BGlu_1.5/rpkm_pseud$BGlu_0),
                            BGlu_3FC0 = log2(rpkm_pseud$BGlu_3/rpkm_pseud$BGlu_0),
                            BGlu_8FC0 = log2(rpkm_pseud$BGlu_8/rpkm_pseud$BGlu_0))

rpkm_BGlu_FC_PBS_0.7 <- subset(rpkm_BGlu_FC_PBS_raw, 
PBS_1.5FC0 > 0.7 | PBS_3FC0 > 0.7 | PBS_8FC0 > 0.7 | BGlu_1.5FC0 > 0.7 | BGlu_3FC0 > 0.7 | BGlu_8FC0 > 0.7)


#IFNb + PBS log2FC>0.7
rpkm_IFNb_FC_PBS_raw <- data.frame(row.names = rownames(rpkm_pseud),
                                   PBS_1.5FC0 = log2(rpkm_pseud$PBS_1.5/rpkm_pseud$PBS_0),
                                   PBS_3FC0 = log2(rpkm_pseud$PBS_3/rpkm_pseud$PBS_0),
                                   PBS_8FC0 = log2(rpkm_pseud$PBS_8/rpkm_pseud$PBS_0),
                                   IFNb_1.5FC0 = log2(rpkm_pseud$IFNb_1.5/rpkm_pseud$IFNb_0),
                                   IFNb_3FC0 = log2(rpkm_pseud$IFNb_3/rpkm_pseud$IFNb_0),
                                   IFNb_8FC0 = log2(rpkm_pseud$IFNb_8/rpkm_pseud$IFNb_0))
rpkm_IFNb_FC_PBS_0.7 <- subset(rpkm_IFNb_FC_PBS_raw, PBS_1.5FC0 > 0.7 | PBS_3FC0 > 0.7 | PBS_8FC0 > 0.7 | IFNb_1.5FC0 > 0.7 | IFNb_3FC0 > 0.7 | IFNb_8FC0 > 0.7)

#IFNg + PBS log2FC>0.7
rpkm_IFNg_FC_PBS_raw <- data.frame(row.names = rownames(rpkm_pseud),
                                   PBS_1.5FC0 = log2(rpkm_pseud$PBS_1.5/rpkm_pseud$PBS_0),
                                   PBS_3FC0 = log2(rpkm_pseud$PBS_3/rpkm_pseud$PBS_0),
                                   PBS_8FC0 = log2(rpkm_pseud$PBS_8/rpkm_pseud$PBS_0),
                                   IFNg_1.5FC0 = log2(rpkm_pseud$IFNg_1.5/rpkm_pseud$IFNg_0),
                                   IFNg_3FC0 = log2(rpkm_pseud$IFNg_3/rpkm_pseud$IFNg_0),
                                   IFNg_8FC0 = log2(rpkm_pseud$IFNg_8/rpkm_pseud$IFNg_0))
rpkm_IFNg_FC_PBS_0.7 <- subset(rpkm_IFNg_FC_PBS_raw, PBS_1.5FC0 > 0.7 | PBS_3FC0 > 0.7 | PBS_8FC0 > 0.7 | IFNg_1.5FC0 > 0.7 | IFNg_3FC0 > 0.7 | IFNg_8FC0 > 0.7)


#LPS + PBS log2FC>0.7
rpkm_LPS_FC_PBS_raw <- data.frame(row.names = rownames(rpkm_pseud),
                                   PBS_1.5FC0 = log2(rpkm_pseud$PBS_1.5/rpkm_pseud$PBS_0),
                                   PBS_3FC0 = log2(rpkm_pseud$PBS_3/rpkm_pseud$PBS_0),
                                   PBS_8FC0 = log2(rpkm_pseud$PBS_8/rpkm_pseud$PBS_0),
                                   LPS_1.5FC0 = log2(rpkm_pseud$LPS_1.5/rpkm_pseud$LPS_0),
                                   LPS_3FC0 = log2(rpkm_pseud$LPS_3/rpkm_pseud$LPS_0),
                                   LPS_8FC0 = log2(rpkm_pseud$LPS_8/rpkm_pseud$LPS_0))
rpkm_LPS_FC_PBS_0.7 <- subset(rpkm_LPS_FC_PBS_raw, PBS_1.5FC0 > 0.7 | PBS_3FC0 > 0.7 | PBS_8FC0 > 0.7 | LPS_1.5FC0 > 0.7 | LPS_3FC0 > 0.7 | LPS_8FC0 > 0.7)



#Create 4 dataframes for TPM (PBS+BGlu, PBS+IFNb, PBS+IFNg, PBS+LPS) with genes with at least one log2FC value >0.5
#Starting 10/25/2022, tested a few cutoff values from 0.5, 1, 2
#BGlu + PBS log2FC>0.5
tpm_BGlu_FC_PBS_raw <- data.frame(row.names = rownames(tpm_pseud),
                                   PBS_1.5FC0 = log2(tpm_pseud$PBS_1.5/tpm_pseud$PBS_0),
                                   PBS_3FC0 = log2(tpm_pseud$PBS_3/tpm_pseud$PBS_0),
                                   PBS_8FC0 = log2(tpm_pseud$PBS_8/tpm_pseud$PBS_0),
                                   BGlu_1.5FC0 = log2(tpm_pseud$BGlu_1.5/tpm_pseud$BGlu_0),
                                   BGlu_3FC0 = log2(tpm_pseud$BGlu_3/tpm_pseud$BGlu_0),
                                   BGlu_8FC0 = log2(tpm_pseud$BGlu_8/tpm_pseud$BGlu_0))
tpm_BGlu_FC_PBS_0.5 <- subset(tpm_BGlu_FC_PBS_raw, PBS_1.5FC0 > 0.5 | PBS_3FC0 > 0.5 | PBS_8FC0 > 0.5 | BGlu_1.5FC0 > 0.5 | BGlu_3FC0 > 0.5 | BGlu_8FC0 > 0.5)
tpm_BGlu_FC_PBS_1 <- subset(tpm_BGlu_FC_PBS_raw, PBS_1.5FC0 > 1 | PBS_3FC0 > 1 | PBS_8FC0 > 1 | BGlu_1.5FC0 > 1 | BGlu_3FC0 > 1 | BGlu_8FC0 > 1)
tpm_BGlu_FC_PBS_2 <- subset(tpm_BGlu_FC_PBS_raw, PBS_1.5FC0 > 2 | PBS_3FC0 > 2 | PBS_8FC0 > 2 | BGlu_1.5FC0 > 2 | BGlu_3FC0 > 2 | BGlu_8FC0 > 2)


#IFNb + PBS log2FC>0.5
tpm_IFNb_FC_PBS_raw <- data.frame(row.names = rownames(tpm_pseud),
                                   PBS_1.5FC0 = log2(tpm_pseud$PBS_1.5/tpm_pseud$PBS_0),
                                   PBS_3FC0 = log2(tpm_pseud$PBS_3/tpm_pseud$PBS_0),
                                   PBS_8FC0 = log2(tpm_pseud$PBS_8/tpm_pseud$PBS_0),
                                   IFNb_1.5FC0 = log2(tpm_pseud$IFNb_1.5/tpm_pseud$IFNb_0),
                                   IFNb_3FC0 = log2(tpm_pseud$IFNb_3/tpm_pseud$IFNb_0),
                                   IFNb_8FC0 = log2(tpm_pseud$IFNb_8/tpm_pseud$IFNb_0))
tpm_IFNb_FC_PBS_0.5 <- subset(tpm_IFNb_FC_PBS_raw, PBS_1.5FC0 > 0.5 | PBS_3FC0 > 0.5 | PBS_8FC0 > 0.5 | IFNb_1.5FC0 > 0.5 | IFNb_3FC0 > 0.5 | IFNb_8FC0 > 0.5)
tpm_IFNb_FC_PBS_1 <- subset(tpm_IFNb_FC_PBS_raw, PBS_1.5FC0 > 1 | PBS_3FC0 > 1 | PBS_8FC0 > 1 | IFNb_1.5FC0 > 1 | IFNb_3FC0 > 1 | IFNb_8FC0 > 1)
tpm_IFNb_FC_PBS_2 <- subset(tpm_IFNb_FC_PBS_raw, PBS_1.5FC0 > 2 | PBS_3FC0 > 2 | PBS_8FC0 > 2 | IFNb_1.5FC0 > 2 | IFNb_3FC0 > 2 | IFNb_8FC0 > 2)


#IFNg + PBS log2FC>0.5
tpm_IFNg_FC_PBS_raw <- data.frame(row.names = rownames(tpm_pseud),
                                   PBS_1.5FC0 = log2(tpm_pseud$PBS_1.5/tpm_pseud$PBS_0),
                                   PBS_3FC0 = log2(tpm_pseud$PBS_3/tpm_pseud$PBS_0),
                                   PBS_8FC0 = log2(tpm_pseud$PBS_8/tpm_pseud$PBS_0),
                                   IFNg_1.5FC0 = log2(tpm_pseud$IFNg_1.5/tpm_pseud$IFNg_0),
                                   IFNg_3FC0 = log2(tpm_pseud$IFNg_3/tpm_pseud$IFNg_0),
                                   IFNg_8FC0 = log2(tpm_pseud$IFNg_8/tpm_pseud$IFNg_0))
tpm_IFNg_FC_PBS_0.5 <- subset(tpm_IFNg_FC_PBS_raw, PBS_1.5FC0 > 0.5 | PBS_3FC0 > 0.5 | PBS_8FC0 > 0.5 | IFNg_1.5FC0 > 0.5 | IFNg_3FC0 > 0.5 | IFNg_8FC0 > 0.5)
tpm_IFNg_FC_PBS_1 <- subset(tpm_IFNg_FC_PBS_raw, PBS_1.5FC0 > 1 | PBS_3FC0 > 1 | PBS_8FC0 > 1 | IFNg_1.5FC0 > 1 | IFNg_3FC0 > 1 | IFNg_8FC0 > 1)
tpm_IFNg_FC_PBS_2 <- subset(tpm_IFNg_FC_PBS_raw, PBS_1.5FC0 > 2 | PBS_3FC0 > 2 | PBS_8FC0 > 2 | IFNg_1.5FC0 > 2 | IFNg_3FC0 > 2 | IFNg_8FC0 > 2)



#LPS + PBS log2FC>0.5
tpm_LPS_FC_PBS_raw <- data.frame(row.names = rownames(tpm_pseud),
                                  PBS_1.5FC0 = log2(tpm_pseud$PBS_1.5/tpm_pseud$PBS_0),
                                  PBS_3FC0 = log2(tpm_pseud$PBS_3/tpm_pseud$PBS_0),
                                  PBS_8FC0 = log2(tpm_pseud$PBS_8/tpm_pseud$PBS_0),
                                  LPS_1.5FC0 = log2(tpm_pseud$LPS_1.5/tpm_pseud$LPS_0),
                                  LPS_3FC0 = log2(tpm_pseud$LPS_3/tpm_pseud$LPS_0),
                                  LPS_8FC0 = log2(tpm_pseud$LPS_8/tpm_pseud$LPS_0))
tpm_LPS_FC_PBS_0.5 <- subset(tpm_LPS_FC_PBS_raw, PBS_1.5FC0 > 0.5 | PBS_3FC0 > 0.5 | PBS_8FC0 > 0.5 | LPS_1.5FC0 > 0.5 | LPS_3FC0 > 0.5 | LPS_8FC0 > 0.5)
tpm_LPS_FC_PBS_1 <- subset(tpm_LPS_FC_PBS_raw, PBS_1.5FC0 > 1 | PBS_3FC0 > 1 | PBS_8FC0 > 1 | LPS_1.5FC0 > 1 | LPS_3FC0 > 1 | LPS_8FC0 > 1)
tpm_LPS_FC_PBS_2 <- subset(tpm_LPS_FC_PBS_raw, PBS_1.5FC0 > 2 | PBS_3FC0 > 2 | PBS_8FC0 > 2 | LPS_1.5FC0 > 2 | LPS_3FC0 > 2 | LPS_8FC0 > 2)

#Generated FC dataframes with various induced gene cutoffs
rpkm_BGlu_FC_PBS_raw
rpkm_IFNb_FC_PBS_raw
rpkm_IFNg_FC_PBS_raw
rpkm_LPS_FC_PBS_raw
tpm_BGlu_FC_PBS_raw
tpm_IFNb_FC_PBS_raw
tpm_IFNg_FC_PBS_raw
tpm_LPS_FC_PBS_raw

#histogram for raw FC RPKM Data
#**melted dataframes have no rownames/gene symbols!!
FC_melted <- melt(rpkm_LPS_FC_PBS_raw, variable.name = "StimulusCondition", value.name = "GeneExpression", row.names = TRUE)

ggplot(FC_melted, aes(x = GeneExpression))+
  geom_histogram(bins=70) + 
  labs(x = "Gene Expression (log2FC, RPKM)", y = "Counts", title = "RPKM Log2FC (LPS+PBS) Distribution") +
  scale_x_continuous(breaks = seq(-2, 2, .5),lim = c(-2, 2))# +ylim(0,4500)

#IFN\U03b2

#histogram for raw FC TPM Data
#**melted dataframes have no rownames/gene symbols!!
FC_melted <- melt(tpm_BGlu_FC_PBS_raw, variable.name = "StimulusCondition", value.name = "GeneExpression", row.names = TRUE)

ggplot(FC_melted, aes(x = GeneExpression))+
  geom_histogram(bins=70) + 
  labs(x = "Gene Expression (log2FC, TPM)", y = "Counts", title = "TPM Log2FC (BGlu+PBS) Distribution") +
  scale_x_continuous(breaks = seq(-2, 2, .5),lim = c(-2, 2)) #ylim(0,4500)



#Using the the cutoff for rpkm and tpm FC values above, subset genes from rpkm and tpm
# and select for genes w at least one value>2

# generate dataframe with only induced and expressed genes for BGlu and PBS [RPKM]
#extract induced row names
induced_rpkm_BGlu_PBS_symbols <- rownames(rpkm_BGlu_FC_PBS_0.7)
#select for all rpkm induced genes
rpkm_induced_BGlu_PBS_all <- bmdm_rpkm[rownames(bmdm_rpkm) %in% induced_rpkm_BGlu_PBS_symbols,]
#create dataframe with only bglu/pbs induced genes
rpkm_induced_BGlu_PBS <- data.frame(row.names = rownames(rpkm_induced_BGlu_PBS_all),
                                    PBS_0 = rpkm_induced_BGlu_PBS_all$PBS_0,
                                    PBS_1.5 = rpkm_induced_BGlu_PBS_all$PBS_1.5,
                                    PBS_3 = rpkm_induced_BGlu_PBS_all$PBS_3,
                                    PBS_8 = rpkm_induced_BGlu_PBS_all$PBS_8,
                                    BGlu_0 = rpkm_induced_BGlu_PBS_all$BGlu_0,
                                    BGlu_1.5 = rpkm_induced_BGlu_PBS_all$BGlu_1.5,
                                    BGlu_3 = rpkm_induced_BGlu_PBS_all$BGlu_3,
                                    BGlu_8 = rpkm_induced_BGlu_PBS_all$BGlu_8)
#dataframe with only induced AND expressed genes for BGlu and PBS [RPKM]
rpkm_BGlu_PBS <- subset(rpkm_induced_BGlu_PBS, 
PBS_0 > 2 | PBS_1.5 > 2 | PBS_3 > 2 | PBS_8 > 2 | BGlu_0 > 2 | BGlu_1.5 > 2 | BGlu_3 > 2 | BGlu_8 >2)



# generate dataframe with only induced and expressed genes for IFNb and PBS [RPKM]
#extract induced row names
induced_rpkm_IFNb_PBS_symbols <- rownames(rpkm_IFNb_FC_PBS_0.7)
#select for all rpkm induced genes
rpkm_induced_IFNb_PBS_all <- bmdm_rpkm[rownames(bmdm_rpkm) %in% induced_rpkm_IFNb_PBS_symbols,]
#create dataframe with only IFNb/pbs induced genes
rpkm_induced_IFNb_PBS <- data.frame(row.names = rownames(rpkm_induced_IFNb_PBS_all),
                                    PBS_0 = rpkm_induced_IFNb_PBS_all$PBS_0,
                                    PBS_1.5 = rpkm_induced_IFNb_PBS_all$PBS_1.5,
                                    PBS_3 = rpkm_induced_IFNb_PBS_all$PBS_3,
                                    PBS_8 = rpkm_induced_IFNb_PBS_all$PBS_8,
                                    IFNb_0 = rpkm_induced_IFNb_PBS_all$IFNb_0,
                                    IFNb_1.5 = rpkm_induced_IFNb_PBS_all$IFNb_1.5,
                                    IFNb_3 = rpkm_induced_IFNb_PBS_all$IFNb_3,
                                    IFNb_8 = rpkm_induced_IFNb_PBS_all$IFNb_8)
#dataframe with only induced AND expressed genes for IFNb and PBS [RPKM]
rpkm_IFNb_PBS <- subset(rpkm_induced_IFNb_PBS, PBS_0 > 2 | PBS_1.5 > 2 | PBS_3 > 2 | PBS_8 > 2 | IFNb_0 > 2 | IFNb_1.5 > 2 | IFNb_3 > 2 | IFNb_8 >2)



# generate dataframe with only induced and expressed genes for IFNg and PBS [RPKM]
#extract induced row names
induced_rpkm_IFNg_PBS_symbols <- rownames(rpkm_IFNg_FC_PBS_0.7)
#select for all rpkm induced genes
rpkm_induced_IFNg_PBS_all <- bmdm_rpkm[rownames(bmdm_rpkm) %in% induced_rpkm_IFNg_PBS_symbols,]
#create dataframe with only IFNg/pbs induced genes
rpkm_induced_IFNg_PBS <- data.frame(row.names = rownames(rpkm_induced_IFNg_PBS_all),
                                    PBS_0 = rpkm_induced_IFNg_PBS_all$PBS_0,
                                    PBS_1.5 = rpkm_induced_IFNg_PBS_all$PBS_1.5,
                                    PBS_3 = rpkm_induced_IFNg_PBS_all$PBS_3,
                                    PBS_8 = rpkm_induced_IFNg_PBS_all$PBS_8,
                                    IFNg_0 = rpkm_induced_IFNg_PBS_all$IFNg_0,
                                    IFNg_1.5 = rpkm_induced_IFNg_PBS_all$IFNg_1.5,
                                    IFNg_3 = rpkm_induced_IFNg_PBS_all$IFNg_3,
                                    IFNg_8 = rpkm_induced_IFNg_PBS_all$IFNg_8)
#dataframe with only induced AND expressed genes for IFNg and PBS [RPKM]
rpkm_IFNg_PBS <- subset(rpkm_induced_IFNg_PBS, PBS_0 > 2 | PBS_1.5 > 2 | PBS_3 > 2 | PBS_8 > 2 | IFNg_0 > 2 | IFNg_1.5 > 2 | IFNg_3 > 2 | IFNg_8 >2)



# generate dataframe with only induced and expressed genes for LPS and PBS [RPKM]
#extract induced row names
induced_rpkm_LPS_PBS_symbols <- rownames(rpkm_LPS_FC_PBS_0.7)
#select for all rpkm induced genes
rpkm_induced_LPS_PBS_all <- bmdm_rpkm[rownames(bmdm_rpkm) %in% induced_rpkm_LPS_PBS_symbols,]
#create dataframe with only LPS/pbs induced genes
rpkm_induced_LPS_PBS <- data.frame(row.names = rownames(rpkm_induced_LPS_PBS_all),
                                    PBS_0 = rpkm_induced_LPS_PBS_all$PBS_0,
                                    PBS_1.5 = rpkm_induced_LPS_PBS_all$PBS_1.5,
                                    PBS_3 = rpkm_induced_LPS_PBS_all$PBS_3,
                                    PBS_8 = rpkm_induced_LPS_PBS_all$PBS_8,
                                    LPS_0 = rpkm_induced_LPS_PBS_all$LPS_0,
                                    LPS_1.5 = rpkm_induced_LPS_PBS_all$LPS_1.5,
                                    LPS_3 = rpkm_induced_LPS_PBS_all$LPS_3,
                                    LPS_8 = rpkm_induced_LPS_PBS_all$LPS_8)
#dataframe with only induced AND expressed genes for LPS and PBS [RPKM]
rpkm_LPS_PBS <- subset(rpkm_induced_LPS_PBS, PBS_0 > 2 | PBS_1.5 > 2 | PBS_3 > 2 | PBS_8 > 2 | LPS_0 > 2 | LPS_1.5 > 2 | LPS_3 > 2 | LPS_8 >2)


#10/18/2022 next steps, applied a cutoff of 10 rather than 2 :)
# generate dataframe with only induced and expressed genes for BGlu and PBS [TPM]
#extract induced row names
induced_tpm_BGlu_PBS_symbols <- rownames(tpm_BGlu_FC_PBS_0.5)
#select for all tpm induced genes
tpm_induced_BGlu_PBS_all <- bmdm_tpm[rownames(bmdm_tpm) %in% induced_tpm_BGlu_PBS_symbols,]
#create dataframe with only bglu/pbs induced genes
tpm_induced_BGlu_PBS <- data.frame(row.names = rownames(tpm_induced_BGlu_PBS_all),
                                    PBS_0 = tpm_induced_BGlu_PBS_all$PBS_0,
                                    PBS_1.5 = tpm_induced_BGlu_PBS_all$PBS_1.5,
                                    PBS_3 = tpm_induced_BGlu_PBS_all$PBS_3,
                                    PBS_8 = tpm_induced_BGlu_PBS_all$PBS_8,
                                    BGlu_0 = tpm_induced_BGlu_PBS_all$BGlu_0,
                                    BGlu_1.5 = tpm_induced_BGlu_PBS_all$BGlu_1.5,
                                    BGlu_3 = tpm_induced_BGlu_PBS_all$BGlu_3,
                                    BGlu_8 = tpm_induced_BGlu_PBS_all$BGlu_8)
#dataframe with only induced AND expressed genes for BGlu and PBS [TPM]
tpm_BGlu_PBS <- subset(tpm_induced_BGlu_PBS, PBS_0 > 10 | PBS_1.5 > 10 | PBS_3 > 10 | PBS_8 > 10 | BGlu_0 > 10 | BGlu_1.5 > 10 | BGlu_3 > 10 | BGlu_8 >10)



# generate dataframe with only induced and expressed genes for IFNb and PBS [TPM]
#extract induced row names
induced_tpm_IFNb_PBS_symbols <- rownames(tpm_IFNb_FC_PBS_0.5)
#select for all tpm induced genes
tpm_induced_IFNb_PBS_all <- bmdm_tpm[rownames(bmdm_tpm) %in% induced_tpm_IFNb_PBS_symbols,]
#create dataframe with only IFNb/pbs induced genes
tpm_induced_IFNb_PBS <- data.frame(row.names = rownames(tpm_induced_IFNb_PBS_all),
                                    PBS_0 = tpm_induced_IFNb_PBS_all$PBS_0,
                                    PBS_1.5 = tpm_induced_IFNb_PBS_all$PBS_1.5,
                                    PBS_3 = tpm_induced_IFNb_PBS_all$PBS_3,
                                    PBS_8 = tpm_induced_IFNb_PBS_all$PBS_8,
                                    IFNb_0 = tpm_induced_IFNb_PBS_all$IFNb_0,
                                    IFNb_1.5 = tpm_induced_IFNb_PBS_all$IFNb_1.5,
                                    IFNb_3 = tpm_induced_IFNb_PBS_all$IFNb_3,
                                    IFNb_8 = tpm_induced_IFNb_PBS_all$IFNb_8)
#dataframe with only induced AND expressed genes for IFNb and PBS [TPM]
tpm_IFNb_PBS <- subset(tpm_induced_IFNb_PBS, PBS_0 > 10 | PBS_1.5 > 10 | PBS_3 > 10 | PBS_8 > 10 | IFNb_0 > 10 | IFNb_1.5 > 10 | IFNb_3 > 10 | IFNb_8 >10)




# generate dataframe with only induced and expressed genes for IFNg and PBS [TPM]
#extract induced row names
induced_tpm_IFNg_PBS_symbols <- rownames(tpm_IFNg_FC_PBS_0.5)
#select for all tpm induced genes
tpm_induced_IFNg_PBS_all <- bmdm_tpm[rownames(bmdm_tpm) %in% induced_tpm_IFNg_PBS_symbols,]
#create dataframe with only IFNg/pbs induced genes
tpm_induced_IFNg_PBS <- data.frame(row.names = rownames(tpm_induced_IFNg_PBS_all),
                                    PBS_0 = tpm_induced_IFNg_PBS_all$PBS_0,
                                    PBS_1.5 = tpm_induced_IFNg_PBS_all$PBS_1.5,
                                    PBS_3 = tpm_induced_IFNg_PBS_all$PBS_3,
                                    PBS_8 = tpm_induced_IFNg_PBS_all$PBS_8,
                                    IFNg_0 = tpm_induced_IFNg_PBS_all$IFNg_0,
                                    IFNg_1.5 = tpm_induced_IFNg_PBS_all$IFNg_1.5,
                                    IFNg_3 = tpm_induced_IFNg_PBS_all$IFNg_3,
                                    IFNg_8 = tpm_induced_IFNg_PBS_all$IFNg_8)
#dataframe with only induced AND expressed genes for IFNg and PBS [TPM]
tpm_IFNg_PBS <- subset(tpm_induced_IFNg_PBS, PBS_0 > 10 | PBS_1.5 > 10 | PBS_3 > 10 | PBS_8 > 10 | IFNg_0 > 10 | IFNg_1.5 > 10 | IFNg_3 > 10 | IFNg_8 >10)




# generate dataframe with only induced and expressed genes for LPS and PBS [TPM]
#extract induced row names
induced_tpm_LPS_PBS_symbols <- rownames(tpm_LPS_FC_PBS_0.5)
#select for all tpm induced genes
tpm_induced_LPS_PBS_all <- bmdm_tpm[rownames(bmdm_tpm) %in% induced_tpm_LPS_PBS_symbols,]
#create dataframe with only LPS/pbs induced genes
tpm_induced_LPS_PBS <- data.frame(row.names = rownames(tpm_induced_LPS_PBS_all),
                                    PBS_0 = tpm_induced_LPS_PBS_all$PBS_0,
                                    PBS_1.5 = tpm_induced_LPS_PBS_all$PBS_1.5,
                                    PBS_3 = tpm_induced_LPS_PBS_all$PBS_3,
                                    PBS_8 = tpm_induced_LPS_PBS_all$PBS_8,
                                    LPS_0 = tpm_induced_LPS_PBS_all$LPS_0,
                                    LPS_1.5 = tpm_induced_LPS_PBS_all$LPS_1.5,
                                    LPS_3 = tpm_induced_LPS_PBS_all$LPS_3,
                                    LPS_8 = tpm_induced_LPS_PBS_all$LPS_8)
#dataframe with only induced AND expressed genes for LPS and PBS [TPM]
tpm_LPS_PBS <- subset(tpm_induced_LPS_PBS, PBS_0 > 10 | PBS_1.5 > 10 | PBS_3 > 10 | PBS_8 > 10 | LPS_0 > 10 | LPS_1.5 > 10 | LPS_3 > 10 | LPS_8 >10)

#Generated dataframes with cutoffs applied
rpkm_BGlu_PBS #1127 genes
rpkm_IFNb_PBS  #1019 genes
rpkm_IFNg_PBS  #875 genes
rpkm_LPS_PBS #1199 genes
tpm_BGlu_PBS #3580 genes  #2037 genes with cutoff of 10
tpm_IFNb_PBS #3459 genes  #1767 genes with cutoff of 10
tpm_IFNg_PBS  #2417 genes #1161 genes with cutoff of 10
tpm_LPS_PBS #3382 genes   #1788 genes with cutoff of 10

nrow(tpm_LPS_PBS)


#Remove mitochondrial genes from TPM datasets

#Extract rownames into a column of symbols
tpm_BGlu_PBS$Symbol <-rownames(tpm_BGlu_PBS)
tpm_IFNb_PBS$Symbol <-rownames(tpm_IFNb_PBS)
tpm_IFNg_PBS$Symbol <-rownames(tpm_IFNg_PBS)
tpm_LPS_PBS$Symbol <-rownames(tpm_LPS_PBS)

#Subset dataframes by removing gene rows with symbols of mitochrondrial genes "Mt-*"



tpm_BGlu_PBS <- subset(tpm_BGlu_PBS, !tpm_BGlu_PBS$Symbol %in% c("mt-Tw","mt-Tt","mt-Tq","mt-Tm","mt-Ti","mt-Nd5","mt-Nd4","mt-Nd2","mt-Nd1","mt-Cytb","mt-Co1"))
tpm_IFNb_PBS <- subset(tpm_IFNb_PBS, !tpm_IFNb_PBS$Symbol %in% c("mt-Tl1","mt-Ti","mt-Tq","mt-Tm","mt-Tw","mt-Tt","mt-Tp"))
tpm_IFNg_PBS <- subset(tpm_IFNg_PBS, !tpm_IFNg_PBS$Symbol %in% c("mt-Ti","mt-Tq","mt-Tm","mt-Tw","mt-Tn"))
tpm_LPS_PBS <- subset(tpm_LPS_PBS, !tpm_LPS_PBS$Symbol %in% c("mt-Tl1","mt-Nd1","mt-Ti","mt-Tq","mt-Nd2","mt-Tw","mt-Co1","mt-Nd4","mt-Nd5","mt-Nd6","mt-Cytb","mt-Tt","mt-Tp"))

nrow(tpm_LPS_PBS)

tpm_BGlu_PBS <- subset(tpm_BGlu_PBS, select = -c(Symbol))
tpm_IFNb_PBS <- subset(tpm_IFNb_PBS, select = -c(Symbol))
tpm_IFNg_PBS <- subset(tpm_IFNg_PBS, select = -c(Symbol))
tpm_LPS_PBS <- subset(tpm_LPS_PBS, select = -c(Symbol))


  

  ###new dataframes generated - now we want to look at heatmaps


#######################
#Generating heatmaps :)

##Heatmaps - just replace __________down here\/ and correct k = __ cluster value and heatmap title

##HEATMAP GENERATION
##Cluster data
scaled <- as.data.frame(t(scale(t(tpm_IFNg_PBS))))
set.seed(1) 
k = 4
kmcluster <- kmeans(scaled, 
                    iter.max = 1000, 
                    centers = k, 
                    nstart = 1)

##Extract the cluster information and add to the scaled dataframe
mat <- data.frame(cbind(scaled, cluster = kmcluster$cluster))
mat <- mat[order(mat$cluster), ]

##Create breaks on the heatmap and separate the clusters 
count <- k
for (i in 1:k) {
  count[i] <- length(kmcluster$cluster[kmcluster$cluster == i])
}
rowseps <- cumsum(count)

##Plotting the heatmap 
pheatmap(mat[,1:8], 
         cluster_rows = F, 
         cluster_cols = F, 
         show_colnames = T, 
         show_rownames = F, 
         scale = "none", 
         color = colorRampPalette(c("navy", "white", "red"))(50),
         gaps_row = rowseps,
         gaps_col = c(4),
         main = "IFN\U03b3 TPM 4 clusters, Version 4")

#IFN\U03b2



#Generating dataframe for clustered genes - symbol, tpm count, cluster

#TPM BGlu_PBS

scaled <- as.data.frame(t(scale(t(tpm_BGlu_PBS))))
set.seed(1) 
nrow(scaled)
k = 3
kmcluster <- kmeans(scaled, 
                    iter.max = 1000, 
                    centers = k, 
                    nstart = 1)

mat_tpm_BGlu <- data.frame(cbind(scaled, cluster = kmcluster$cluster))
mat_tpm_BGlu <- mat_tpm_BGlu[order(mat_tpm_BGlu$cluster), ]


#create data frame with gene cluster values
clustered_tpm_BGlu <- data.frame(Cluster = mat_tpm_BGlu$cluster,
                   tpm_BGlu_PBS[rownames(mat_tpm_BGlu),])  #taking the order from the mat_tpm_BGlu

match(rownames((clustered_tpm_BGlu)), rownames(mat_tpm_BGlu)) #check if rownames are in the same order



#TPM IFNb_PBS

scaled <- as.data.frame(t(scale(t(tpm_IFNb_PBS))))
set.seed(1) 
nrow(scaled)
k = 5
kmcluster <- kmeans(scaled, 
                    iter.max = 1000, 
                    centers = k, 
                    nstart = 1)

mat_tpm_IFNb <- data.frame(cbind(scaled, cluster = kmcluster$cluster))
mat_tpm_IFNb <- mat_tpm_IFNb[order(mat_tpm_IFNb$cluster), ]

#create data frame with gene cluster values
clustered_tpm_IFNb <- data.frame(Cluster = mat_tpm_IFNb$cluster,
                               tpm_IFNb_PBS[rownames(mat_tpm_IFNb),])  #taking the order from the mat_tpm_IFNb

match(rownames((clustered_tpm_IFNb)), rownames(mat_tpm_IFNb)) #check if rownames are in the same order


#TPM IFNg_PBS
scaled <- as.data.frame(t(scale(t(tpm_IFNg_PBS))))
set.seed(1) 
nrow(scaled)
k = 3
kmcluster <- kmeans(scaled, 
                    iter.max = 1000, 
                    centers = k, 
                    nstart = 1)

mat_tpm_IFNg <- data.frame(cbind(scaled, cluster = kmcluster$cluster))
mat_tpm_IFNg <- mat_tpm_IFNg[order(mat_tpm_IFNg$cluster), ]

#create data frame with gene cluster values
clustered_tpm_IFNg <- data.frame(Cluster = mat_tpm_IFNg$cluster,
                               tpm_IFNg_PBS[rownames(mat_tpm_IFNg),])  #taking the order from the mat_tpm_IFNg

match(rownames((clustered_tpm_IFNg)), rownames(mat_tpm_IFNg)) #check if rownames are in the same order


#TPM LPS_PBS
scaled <- as.data.frame(t(scale(t(tpm_LPS_PBS))))
set.seed(1) 
nrow(scaled)
k = 4
kmcluster <- kmeans(scaled, 
                    iter.max = 1000, 
                    centers = k, 
                    nstart = 1)

mat_tpm_LPS <- data.frame(cbind(scaled, cluster = kmcluster$cluster))
mat_tpm_LPS <- mat_tpm_LPS[order(mat_tpm_LPS$cluster), ]

#create data frame with gene cluster values
clustered_tpm_LPS <- data.frame(Cluster = mat_tpm_LPS$cluster,
                               tpm_LPS_PBS[rownames(mat_tpm_LPS),])  #taking the order from the mat_tpm_LPS

match(rownames((clustered_tpm_LPS)), rownames(mat_tpm_LPS)) #check if rownames are in the same order


#Generated dataframes of clustered genes, cutoffs applied already
clustered_tpm_BGlu
clustered_tpm_IFNb
clustered_tpm_IFNg
clustered_tpm_LPS

write.csv(clustered_tpm_BGlu,"C:/Users/am331/OneDrive/Documents/UCLA/Signaling_Systems_Lab/RNASeq_2022/Clustered_TPM_BGlu_V2.csv", row.names = TRUE)
write.csv(clustered_tpm_IFNb,"C:/Users/am331/OneDrive/Documents/UCLA/Signaling_Systems_Lab/RNASeq_2022/Clustered_TPM_IFNb_V2.csv", row.names = TRUE)
write.csv(clustered_tpm_IFNg,"C:/Users/am331/OneDrive/Documents/UCLA/Signaling_Systems_Lab/RNASeq_2022/Clustered_TPM_IFNg_V2.csv", row.names = TRUE)
write.csv(clustered_tpm_LPS,"C:/Users/am331/OneDrive/Documents/UCLA/Signaling_Systems_Lab/RNASeq_2022/Clustered_TPM_LPS_V2.csv", row.names = TRUE)


##subset LPS clustered dataset for only cluster 3 to plot in a heatmap

LPS_cluster_3 <- subset(clustered_tpm_LPS, Cluster == 3)
LPS_cluster_3 <- subset(LPS_cluster_3, select = -c(Cluster))


scaled <- as.data.frame(t(scale(t(LPS_cluster_3))))
set.seed(1) 
k = 5
kmcluster <- kmeans(scaled, 
                    iter.max = 1000, 
                    centers = k, 
                    nstart = 1)

##Extract the cluster information and add to the scaled dataframe
mat <- data.frame(cbind(scaled, cluster = kmcluster$cluster))
mat <- mat[order(mat$cluster), ]

##Create breaks on the heatmap and separate the clusters 
count <- k
for (i in 1:k) {
  count[i] <- length(kmcluster$cluster[kmcluster$cluster == i])
}
rowseps <- cumsum(count)

##Plotting the heatmap 
pheatmap(mat[,1:8], 
         cluster_rows = F, 
         cluster_cols = F, 
         show_colnames = T, 
         show_rownames = F, 
         scale = "none", 
         color = colorRampPalette(c("navy", "white", "red"))(50),
         gaps_row = rowseps,
         gaps_col = c(4),
         main = "LPS Version 2, Cluster 3; TPM 5 clusters")















#generate histograms to better determine cutoff values
#**melted dataframes have no rownames/gene symbols!!

#histogram for TPM_BGlu Clustered Data
tpm_BGlu_melted <- melt(clustered_tpm_BGlu, id = c("Cluster"), variable.name = "StimulusCondition", value.name = "GeneExpression", row.names = TRUE)

ggplot(tpm_BGlu_melted, aes(x = GeneExpression))+
  geom_histogram(bins=50) + 
  labs(x = "TPM", y = "Counts", title = "Gene Expression Distribution, \U0392Glu xlim(0,50)") +
  xlim(0,50) #+ylim(0,100)

#histogram for TPM_IFNb Clustered Data
tpm_IFNb_melted <- melt(clustered_tpm_IFNb, id = c("Cluster"), variable.name = "StimulusCondition", value.name = "GeneExpression")

ggplot(tpm_IFNb_melted, aes(x = GeneExpression))+
  geom_histogram(bins=50) + 
  labs(x = "TPM", y = "Counts", title = "Gene Expression Distribution, IFN\U03b2 xlim(0,50)") +
  xlim(0,50) #+ylim(0,100)

#histogram for TPM_IFNg Clustered Data
tpm_IFNg_melted <- melt(clustered_tpm_IFNg, id = c("Cluster"), variable.name = "StimulusCondition", value.name = "GeneExpression")

ggplot(tpm_IFNg_melted, aes(x = GeneExpression))+
  geom_histogram(bins=50) + 
  labs(x = "TPM", y = "Counts", title = "Gene Expression Distributions, IFN\U03b3 xlim(0,50)") +
  xlim(0,50) #+ylim(0,100)

#histogram for TPM_LPS Clustered Data
tpm_LPS_melted <- melt(clustered_tpm_LPS, id = c("Cluster"), variable.name = "StimulusCondition", value.name = "GeneExpression")

ggplot(tpm_LPS_melted, aes(x = GeneExpression))+
  geom_histogram(bins=50) + 
  labs(x = "TPM", y = "Counts", title = "Gene Expression Distributions, LPS xlim(0,50)") +
  xlim(0,50) #+ylim(0,100)

#Generated melted dataframes

tpm_BGlu_melted
tpm_IFNb_melted
tpm_IFNg_melted
tpm_LPS_melted


# Plot expression for one gene of interest, BGlu

gene <- clustered_tpm_BGlu["Ccl3",]

gene <- melt(gene, id = c("Cluster"), variable.name = "StimulusCondition", value.name = "GeneCount_TPM")
gene$Cluster <- NULL

gene <- gene %>% separate(StimulusCondition, c("StimulusCondition", "Timepoint"), sep = "_(?=[^_]+$)") 

ggplot(gene, aes(x=Timepoint, y = GeneCount_TPM, group = StimulusCondition, color = StimulusCondition))+
  geom_line() + geom_point() + labs(title = "Ccl3", y = "Gene Expression (TPM)", x = "Timepoint (hrs)", color = "Stimulus Condition") +
  scale_color_hue(labels = c("\U0392Glu","PBS")) 


# Plot expression for one gene of interest, IFNb

gene <- clustered_tpm_IFNb["Neat1",]

gene <- melt(gene, id = c("Cluster"), variable.name = "StimulusCondition", value.name = "GeneCount_TPM")
gene$Cluster <- NULL

gene <- gene %>% separate(StimulusCondition, c("StimulusCondition", "Timepoint"), sep = "_(?=[^_]+$)") 

ggplot(gene, aes(x=Timepoint, y = GeneCount_TPM, group = StimulusCondition, color = StimulusCondition))+
  geom_line() + geom_point() + labs(title = "Neat1", y = "Gene Expression (TPM)", x = "Timepoint (hrs)", color = "Stimulus Condition") +
  scale_color_hue(labels = c("IFN\U03b2","PBS")) 



# Plot expression for one gene of interest, IFNg

gene <- clustered_tpm_IFNg["Saa3",]

gene <- melt(gene, id = c("Cluster"), variable.name = "StimulusCondition", value.name = "GeneCount_TPM")
gene$Cluster <- NULL

gene <- gene %>% separate(StimulusCondition, c("StimulusCondition", "Timepoint"), sep = "_(?=[^_]+$)") 

ggplot(gene, aes(x=Timepoint, y = GeneCount_TPM, group = StimulusCondition, color = StimulusCondition))+
  geom_line() + geom_point() + labs(title = "Saa3", y = "Gene Expression (TPM)", x = "Timepoint (hrs)", color = "Stimulus Condition") +
  scale_color_hue(labels = c("IFN\U03b3","PBS")) 



# Plot expression for one gene of interest, LPS

gene <- clustered_tpm_LPS["S100a11",]

gene <- melt(gene, id = c("Cluster"), variable.name = "StimulusCondition", value.name = "GeneCount_TPM")
gene$Cluster <- NULL

gene <- gene %>% separate(StimulusCondition, c("StimulusCondition", "Timepoint"), sep = "_(?=[^_]+$)") 

ggplot(gene, aes(x=Timepoint, y = GeneCount_TPM, group = StimulusCondition, color = StimulusCondition))+
  geom_line() + geom_point() + labs(title = "S100a11", y = "Gene Expression (TPM)", x = "Timepoint (hrs)", color = "Stimulus Condition") +
  scale_color_hue(labels = c("LPS","PBS")) 


