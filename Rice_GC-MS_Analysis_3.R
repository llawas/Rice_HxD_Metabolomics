#load required packages and functions
library(VennDiagram)
library(grid)
library(RColorBrewer)
library(pcaMethods)
library(gplots)
library(outliers)
library(stats)
library(nortest)
library(plyr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(pheatmap)
source("H:/LMFLawas/PhD/Experiments/Protocols/R script/RemoveFactors_function.R")
source("H:/LMFLawas/PhD/Experiments/Protocols/R script/func_normalize.R")
source('H:/LMFLawas/PhD/Experiments/Protocols/R script/func_find_one_outlier.R')
source('H:/LMFLawas/PhD/Experiments/Protocols/R script/func_replace_outlier.R')
source('H:/LMFLawas/PhD/Experiments/Protocols/R script/func_hist_outlier.R')
source("H:/LMFLawas/PhD/Experiments/Protocols/R script/func_log_transform.R")
source("H:/LMFLawas/PhD/Experiments/Protocols/R script/func_log2_median_transform.R")


#NOTE: late stress means severe stress


##################################################
#DEVELOPING SEED#
##################################################


#load sample list - contains information on sample ID (Batch & Sequence of running sample in the GC-MS), Treatment, Timepoint, Cultivar, Description (the original sample name in raw data) - prepare in Excel
#order of cultivars: N22, Dular, Anjali
#order of timepoints (EGF): Control-Late stress, Control-12h rewatering (RW), Control-36h RW, Control-60h RW, Late stress, 12h RW, 36h RW, 60h RW
sample_list_developing_seed_2013 = read.table("sample_list_developing_seed_2013.txt", header = T, sep = "\t")
sample_list_developing_seed_2014 = read.table("sample_list_developing_seed_2014.txt", header = T, sep = "\t")
sample_list_developing_seed_2015 = read.table("sample_list_developing_seed_2015.txt", header = T, sep = "\t")

#edited sample list (for ease of use in some subsequent steps) - alternatively, generate a list like this 
sample_list_developing_seed_2013_2 = sample_list_developing_seed_2013
levels(sample_list_developing_seed_2013_2$Timepoint) = gsub("EGF ", "", levels(sample_list_developing_seed_2013_2$Timepoint))
sample_list_developing_seed_2013_2$Trt.Time = ifelse(sample_list_developing_seed_2013_2$Treatment == "Control", 
                                                     yes = as.character(sample_list_developing_seed_2013_2$Treatment), 
                                                     no = as.character(sample_list_developing_seed_2013_2$Timepoint))  #add Trt.Time column (for ease in assigning groups in PCA)
sample_list_developing_seed_2013_2$Trt.Time = as.factor(sample_list_developing_seed_2013_2$Trt.Time)  #convert to factor
sample_list_developing_seed_2013_2$Trt.Time2 = ifelse(sample_list_developing_seed_2013_2$Treatment == "Control", 
                                                      yes = as.character(interaction(sample_list_developing_seed_2013_2$Treatment,
                                                                                     sample_list_developing_seed_2013_2$Timepoint)), 
                                                      no = as.character(sample_list_developing_seed_2013_2$Timepoint))  #for ease in assigning groups in PCA
sample_list_developing_seed_2013_2$Trt.Time2 = as.factor(sample_list_developing_seed_2013_2$Trt.Time2)    #convert to factor

sample_list_developing_seed_2014_2 = sample_list_developing_seed_2014
levels(sample_list_developing_seed_2014_2$Timepoint) = gsub("EGF ", "", levels(sample_list_developing_seed_2014_2$Timepoint))
sample_list_developing_seed_2014_2$Trt.Time = ifelse(sample_list_developing_seed_2014_2$Treatment == "Control", 
                                                     yes = as.character(sample_list_developing_seed_2014_2$Treatment), 
                                                     no = as.character(sample_list_developing_seed_2014_2$Timepoint))  #add Trt.Time column (for ease in assigning groups in PCA)
sample_list_developing_seed_2014_2$Trt.Time = as.factor(sample_list_developing_seed_2014_2$Trt.Time)  #convert to factor
sample_list_developing_seed_2014_2$Trt.Time2 = ifelse(sample_list_developing_seed_2014_2$Treatment == "Control", 
                                                      yes = as.character(interaction(sample_list_developing_seed_2014_2$Treatment,
                                                                                     sample_list_developing_seed_2014_2$Timepoint)), 
                                                      no = as.character(sample_list_developing_seed_2014_2$Timepoint))  #for ease in assigning groups in PCA
sample_list_developing_seed_2014_2$Trt.Time2 = as.factor(sample_list_developing_seed_2014_2$Trt.Time2)    #convert to factor

sample_list_developing_seed_2015_2 = sample_list_developing_seed_2015
levels(sample_list_developing_seed_2015_2$Timepoint) = gsub("EGF ", "", levels(sample_list_developing_seed_2015_2$Timepoint))
sample_list_developing_seed_2015_2$Trt.Time = ifelse(sample_list_developing_seed_2015_2$Treatment == "Control", 
                                                     yes = as.character(sample_list_developing_seed_2015_2$Treatment), 
                                                     no = as.character(sample_list_developing_seed_2015_2$Timepoint))  #add Trt.Time column (for ease in assigning groups in PCA)
sample_list_developing_seed_2015_2$Trt.Time = as.factor(sample_list_developing_seed_2015_2$Trt.Time)  #convert to factor
sample_list_developing_seed_2015_2$Trt.Time2 = ifelse(sample_list_developing_seed_2015_2$Treatment == "Control", 
                                                      yes = as.character(interaction(sample_list_developing_seed_2015_2$Treatment,
                                                                                     sample_list_developing_seed_2015_2$Timepoint)), 
                                                      no = as.character(sample_list_developing_seed_2015_2$Timepoint))  #for ease in assigning groups in PCA
sample_list_developing_seed_2015_2$Trt.Time2 = as.factor(sample_list_developing_seed_2015_2$Trt.Time2)    #convert to factor


#load data
#input data: subset_profile with internal standards and reagents excluded and samples (in columns) aranged in preferred order - done in Excel
#exclude metabolites with > 66.67% NAs across all samples

#2013
subset_profile_developing_seed_HxD_2013 = read.table("subset_profile_developing_seed_2013.txt", sep = "\t", header = T, quote = "")   #load data - 222 metabolites, 120 samples
subset_profile_developing_seed_HxD_2013$Percent_NA = (rowSums(is.na(subset_profile_developing_seed_HxD_2013))/(ncol(subset_profile_developing_seed_HxD_2013)-3))*100    #calculate %NA per metabolite
subset_profile_developing_seed_HxD_2013_filtered = subset_profile_developing_seed_HxD_2013[which(subset_profile_developing_seed_HxD_2013$Percent_NA <= 67.5), ]    #only metabolites with < 66.67% NA but adjusted to 67.5% to include Hydroquinone which is also cultivar-specific but is slightly above the threshold due to two samples with NA - observations from manual inspection important!! - 164 metabolites
sample_list_data_developing_seed_2013 = subset_profile_developing_seed_HxD_2013_filtered[ , -c(1, 2, 3, ncol(subset_profile_developing_seed_HxD_2013_filtered))]   #filtered data set - data values only
metabolite_list_developing_seed_2013 = subset_profile_developing_seed_HxD_2013_filtered[ , c(1, 2, 3)]

#2014
subset_profile_developing_seed_HxD_2014 = read.table("subset_profile_developing_seed_2014.txt", sep = "\t", header = T, quote = "")   #load data - 144 metabolites, 112 samples
subset_profile_developing_seed_HxD_2014$Percent_NA = (rowSums(is.na(subset_profile_developing_seed_HxD_2014))/(ncol(subset_profile_developing_seed_HxD_2014)-3))*100    #calculate %NA per metabolite
subset_profile_developing_seed_HxD_2014_filtered = subset_profile_developing_seed_HxD_2014[which(subset_profile_developing_seed_HxD_2014$Percent_NA < 67.9), ]    #only metabolites with < 66.67% NA but adjusted to 67.9% to include 4-Hydroxyphenyl-beta-glucopyranoside and 4-hydroxy-benzoic acid which are also cultivar-specific but are slightly above the threshold due to unequal number of samples among cultivars - observations from manual inspection important!! - 92 metabolites  --> however, 5 other metabolites have been included but are not cultivar specific, hence, remove in next step

#remove A139005-101, A148006-101, Galactonic acid, Serine, O-acetyl-, and Cinnamic acid, 4-hydroxy-, trans- from list
subset_profile_developing_seed_HxD_2014_filtered_final = subset_profile_developing_seed_HxD_2014_filtered[-which(subset_profile_developing_seed_HxD_2014_filtered$Name == "A139005-101"), ]
subset_profile_developing_seed_HxD_2014_filtered_final = subset_profile_developing_seed_HxD_2014_filtered_final[-which(subset_profile_developing_seed_HxD_2014_filtered_final$Name == "A148006-101"), ]
subset_profile_developing_seed_HxD_2014_filtered_final = subset_profile_developing_seed_HxD_2014_filtered_final[-which(subset_profile_developing_seed_HxD_2014_filtered_final$Name == "Galactonic acid"), ]
subset_profile_developing_seed_HxD_2014_filtered_final = subset_profile_developing_seed_HxD_2014_filtered_final[-which(subset_profile_developing_seed_HxD_2014_filtered_final$Name == "Serine, O-acetyl-"), ]
subset_profile_developing_seed_HxD_2014_filtered_final = subset_profile_developing_seed_HxD_2014_filtered_final[-which(subset_profile_developing_seed_HxD_2014_filtered_final$Name == "Cinnamic acid, 4-hydroxy-, trans-"), ]
#final count: 87 metabolites

sample_list_data_developing_seed_2014 = subset_profile_developing_seed_HxD_2014_filtered_final[ , -c(1, 2, 3, ncol(subset_profile_developing_seed_HxD_2014_filtered_final))]   #filtered data set - data values only
metabolite_list_developing_seed_2014 = subset_profile_developing_seed_HxD_2014_filtered_final[ , c(1, 2, 3)]

#2015
subset_profile_developing_seed_HxD_2015 = read.table("subset_profile_developing_seed_2015.txt", sep = "\t", header = T, quote = "")   #load data - 230 metabolites, 120 samples
subset_profile_developing_seed_HxD_2015$Percent_NA = (rowSums(is.na(subset_profile_developing_seed_HxD_2015))/(ncol(subset_profile_developing_seed_HxD_2015)-3))*100    #calculate %NA per metabolite
subset_profile_developing_seed_HxD_2015_filtered = subset_profile_developing_seed_HxD_2015[which(subset_profile_developing_seed_HxD_2015$Percent_NA < 66.67), ]    #only metabolites with < 66.67% NA - 178 metabolites
sample_list_data_developing_seed_2015 = subset_profile_developing_seed_HxD_2015_filtered[ , -c(1, 2, 3, ncol(subset_profile_developing_seed_HxD_2015_filtered))]   #filtered data set - data values only
metabolite_list_developing_seed_2015 = subset_profile_developing_seed_HxD_2015_filtered[ , c(1, 2, 3)]


#metabolites common across the three experiments
overlap_developing_seed_2013_2014 = intersect(metabolite_list_developing_seed_2013$Name, metabolite_list_developing_seed_2014$Name)
overlap_developing_seed_20132014_2015 = intersect(overlap_developing_seed_2013_2014, metabolite_list_developing_seed_2015$Name)   #68 metabolites

#subset of NA-filtered data with overlapping metabolites only
reduced_metabolite_list_developing_seed_2013 = metabolite_list_developing_seed_2013[which(metabolite_list_developing_seed_2013$Name %in% overlap_developing_seed_20132014_2015), ]
reduced_sample_values_developing_seed_2013 = sample_list_data_developing_seed_2013[which(metabolite_list_developing_seed_2013$Name %in% overlap_developing_seed_20132014_2015), ]
reduced_metabolite_list_developing_seed_2014 = metabolite_list_developing_seed_2014[which(metabolite_list_developing_seed_2014$Name %in% overlap_developing_seed_20132014_2015), ]
reduced_sample_values_developing_seed_2014 = sample_list_data_developing_seed_2014[which(metabolite_list_developing_seed_2014$Name %in% overlap_developing_seed_20132014_2015), ]
reduced_metabolite_list_developing_seed_2015 = metabolite_list_developing_seed_2015[which(metabolite_list_developing_seed_2015$Name %in% overlap_developing_seed_20132014_2015), ]
reduced_sample_values_developing_seed_2015 = sample_list_data_developing_seed_2015[which(metabolite_list_developing_seed_2015$Name %in% overlap_developing_seed_20132014_2015), ]

#remove Ethanolamine (identified as additional contamination) - final number of overlapping metabolites: 67
reduced_metabolite_list_final_developing_seed_2013 = reduced_metabolite_list_developing_seed_2013[-which(reduced_metabolite_list_developing_seed_2013$Name == "Ethanolamine"), ]
reduced_metabolite_list_final_developing_seed_2014 = reduced_metabolite_list_developing_seed_2014[-which(reduced_metabolite_list_developing_seed_2014$Name == "Ethanolamine"), ] 
reduced_metabolite_list_final_developing_seed_2015 = reduced_metabolite_list_developing_seed_2015[-which(reduced_metabolite_list_developing_seed_2015$Name == "Ethanolamine"), ]

reduced_sample_values_final_developing_seed_2013 = reduced_sample_values_developing_seed_2013[-which(reduced_metabolite_list_developing_seed_2013$Name == "Ethanolamine"), ]
reduced_sample_values_final_developing_seed_2014 = reduced_sample_values_developing_seed_2014[-which(reduced_metabolite_list_developing_seed_2014$Name == "Ethanolamine"), ]
reduced_sample_values_final_developing_seed_2015 = reduced_sample_values_developing_seed_2015[-which(reduced_metabolite_list_developing_seed_2015$Name == "Ethanolamine"), ]


class(reduced_metabolite_list_final_developing_seed_2013$Name)    #factor
class(reduced_metabolite_list_final_developing_seed_2014$Name)    #factor
class(reduced_metabolite_list_final_developing_seed_2015$Name)    #factor

#convert metabolite names to character - to be able to check/compare with other lists
reduced_metabolite_list_final_developing_seed_2013$Name = as.character(reduced_metabolite_list_final_developing_seed_2013$Name)
reduced_metabolite_list_final_developing_seed_2014$Name = as.character(reduced_metabolite_list_final_developing_seed_2014$Name)
reduced_metabolite_list_final_developing_seed_2015$Name = as.character(reduced_metabolite_list_final_developing_seed_2015$Name)

class(reduced_metabolite_list_final_developing_seed_2013$Name)    #character
class(reduced_metabolite_list_final_developing_seed_2014$Name)    #character
class(reduced_metabolite_list_final_developing_seed_2015$Name)    #character


#final data structure to use - common metabolites across three experiments
data_overlap_developing_seed_2013 = reduced_sample_values_final_developing_seed_2013
rownames(data_overlap_developing_seed_2013) = reduced_metabolite_list_final_developing_seed_2013$Name
colnames(data_overlap_developing_seed_2013) = sample_list_developing_seed_2013$ID
data_overlap_developing_seed_2013 = t(data_overlap_developing_seed_2013)      #metabolites in columns, samples in rows

data_overlap_developing_seed_2014 = reduced_sample_values_final_developing_seed_2014
rownames(data_overlap_developing_seed_2014) = reduced_metabolite_list_final_developing_seed_2014$Name
colnames(data_overlap_developing_seed_2014) = sample_list_developing_seed_2014$ID
data_overlap_developing_seed_2014 = t(data_overlap_developing_seed_2014)      #metabolites in columns, samples in rows

data_overlap_developing_seed_2015 = reduced_sample_values_final_developing_seed_2015
rownames(data_overlap_developing_seed_2015) = reduced_metabolite_list_final_developing_seed_2015$Name
colnames(data_overlap_developing_seed_2015) = sample_list_developing_seed_2015$ID
data_overlap_developing_seed_2015 = t(data_overlap_developing_seed_2015)      #metabolites in columns, samples in rows


#edit names of MSTs (unknown metabolites) - remove -101
colnames(data_overlap_developing_seed_2013) = gsub("-101", "", colnames(data_overlap_developing_seed_2013))
colnames(data_overlap_developing_seed_2014) = gsub("-101", "", colnames(data_overlap_developing_seed_2014))
colnames(data_overlap_developing_seed_2015) = gsub("-101", "", colnames(data_overlap_developing_seed_2015))

reduced_metabolite_list_final_developing_seed_2013$Name = gsub("-101", "", reduced_metabolite_list_final_developing_seed_2013$Name)
reduced_metabolite_list_final_developing_seed_2014$Name = gsub("-101", "", reduced_metabolite_list_final_developing_seed_2014$Name)
reduced_metabolite_list_final_developing_seed_2015$Name = gsub("-101", "", reduced_metabolite_list_final_developing_seed_2015$Name)


#missing values

#NA count/percentage
(sum(is.na(data_overlap_developing_seed_2013))/(nrow(data_overlap_developing_seed_2013)*ncol(data_overlap_developing_seed_2013)))*100     #5.09%
(sum(is.na(data_overlap_developing_seed_2014))/(nrow(data_overlap_developing_seed_2014)*ncol(data_overlap_developing_seed_2014)))*100   
#14.66%
(sum(is.na(data_overlap_developing_seed_2015))/(nrow(data_overlap_developing_seed_2015)*ncol(data_overlap_developing_seed_2015)))*100   
#0.47%

#visualize NAs
heatmap.2(data_overlap_developing_seed_2013,
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          margins = c(15, 6), 
          main = "Developing seed - 2013 - Missing values", 
          lmat = rbind(c(0, 3, 4), c(2, 1, 0)), 
          lwid = c(1.5, 4, 2), 
          xlab = "Metabolite", 
          ylab = "Sample", 
          na.color = "grey")    

heatmap.2(data_overlap_developing_seed_2014,
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          margins = c(15, 6), 
          main = "Developing seed - 2014 - Missing values", 
          lmat = rbind(c(0, 3, 4), c(2, 1, 0)), 
          lwid = c(1.5, 4, 2), 
          xlab = "Metabolite", 
          ylab = "Sample", 
          na.color = "grey")   

heatmap.2(data_overlap_developing_seed_2015,
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          margins = c(15, 6), 
          main = "Developing seed - 2015 - Missing values", 
          lmat = rbind(c(0, 3, 4), c(2, 1, 0)), 
          lwid = c(1.5, 4, 2), 
          xlab = "Metabolite", 
          ylab = "Sample", 
          na.color = "grey")   


#replace NAs with half of the minimum value per metabolite
func_imp_minimum = function(samples) {
  min_samples = apply(samples,2, min, na.rm = TRUE)
  m_samples = as.matrix(samples)
  for(i in 1:nrow(samples)){
    for(j in 1:ncol(samples)){
      if(is.na(m_samples[i,j]) == TRUE){
        m_samples[i,j] = min_samples[j]/2
      }
    }
  }
  return(m_samples)
}

data_replaced_na_developing_seed_2013 = func_imp_minimum(data_overlap_developing_seed_2013)
data_replaced_na_developing_seed_2014 = func_imp_minimum(data_overlap_developing_seed_2014)
data_replaced_na_developing_seed_2015 = func_imp_minimum(data_overlap_developing_seed_2015)

#check for presence of NAs - should be 0!
sum(is.na(data_replaced_na_developing_seed_2013))
sum(is.na(data_replaced_na_developing_seed_2014))
sum(is.na(data_replaced_na_developing_seed_2015))

#histogram of NA-replaced data
pdf("histogram_NA_replaced_developing_seed_2013.pdf")
for(i in 1:ncol(data_replaced_na_developing_seed_2013)){
  hist(data_replaced_na_developing_seed_2013[ , i], 
       col = "grey", 
       breaks = 20, 
       main = colnames(data_replaced_na_developing_seed_2013)[i])
}
dev.off()

pdf("histogram_NA_replaced_developing_seed_2014.pdf")
for(i in 1:ncol(data_replaced_na_developing_seed_2014)){
  hist(data_replaced_na_developing_seed_2014[ , i], 
       col = "grey", 
       breaks = 20, 
       main = colnames(data_replaced_na_developing_seed_2014)[i])
}
dev.off()

pdf("histogram_NA_replaced_developing_seed_2015.pdf")
for(i in 1:ncol(data_replaced_na_developing_seed_2015)){
  hist(data_replaced_na_developing_seed_2015[ , i], 
       col = "grey", 
       breaks = 20, 
       main = colnames(data_replaced_na_developing_seed_2015)[i])
}
dev.off()


#combine all data sets/sample lists into one

#check first if data files to combine have the same column names
colnames(data_replaced_na_developing_seed_2013) == colnames(data_replaced_na_developing_seed_2014)    #OK
colnames(data_replaced_na_developing_seed_2013) == colnames(data_replaced_na_developing_seed_2015)    #OK

colnames(sample_list_developing_seed_2013) == colnames(sample_list_developing_seed_2014)    #OK
colnames(sample_list_developing_seed_2013) == colnames(sample_list_developing_seed_2015)    #OK


#combine data sets
data_replaced_na_developing_seed_131415 = rbind(data_replaced_na_developing_seed_2013, 
                                                data_replaced_na_developing_seed_2014, 
                                                data_replaced_na_developing_seed_2015)

#combine sample lists
sample_list_developing_seed_131415 = rbind(sample_list_developing_seed_2013, 
                                           sample_list_developing_seed_2014, 
                                           sample_list_developing_seed_2015)

#add additional columns on sample list - combined factors
sample_list_developing_seed_131415$Treatment.Timepoint = interaction(sample_list_developing_seed_131415$Treatment, sample_list_developing_seed_131415$Timepoint)
sample_list_developing_seed_131415$Treatment.Timepoint.Cultivar = interaction(sample_list_developing_seed_131415$Treatment, sample_list_developing_seed_131415$Timepoint, sample_list_developing_seed_131415$Cultivar)

#number of batches - GC-MS runs (combined data)
length(levels(sample_list_developing_seed_131415$Batch))    #33 batches
#since samples were ran on several batches, this might have an effect on the data --> apply ANOVA-based batch effect correction


#PCA before batch correction

#log10 transformation for PCA
data_log10_transformed_developing_seed_131415 = log10(data_replaced_na_developing_seed_131415)

#PCA - centered and scaled
pca_res_developing_seed_20131415_scaled_centered = pca(data_log10_transformed_developing_seed_131415,
                                                       method="ppca", 
                                                       nPcs=5, 
                                                       scale = "pareto", 
                                                       center = T)
pca_res_developing_seed_20131415_scaled_centered
pca_res_developing_seed_20131415_scaled_centered_scores = scores(pca_res_developing_seed_20131415_scaled_centered)

#bar plot of PC variance with values
par(mar = c(5.1, 6.1, 4.1, 2.1))
text(barplot(pca_res_developing_seed_20131415_scaled_centered@R2, 
             names.arg = c("PC1", "PC2", "PC3", "PC4", "PC5"), 
             main = "Probabilistic PCA - Developing seed - Centered, scaled", 
             ylab = "Variance", 
             ylim = c(0, 0.3), 
             cex.axis = 1.5, 
             cex.names = 1.5, 
             cex.lab = 1.5, 
             cex.main = 2), 
     0, round(pca_res_developing_seed_20131415_scaled_centered@R2,3), pos = 3, cex = 1.5)

#PCA score plot
palette(rainbow(length(levels(sample_list_developing_seed_131415$Batch))))
pairs(pca_res_developing_seed_20131415_scaled_centered_scores, 
      col = sample_list_developing_seed_131415$Batch, 
      main = "PC score plots - Developing seed - Centered, scaled")    #pairs

palette(rainbow(length(levels(sample_list_developing_seed_131415$Batch))))
par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = T)
plot(pca_res_developing_seed_20131415_scaled_centered_scores[ , 1], 
     pca_res_developing_seed_20131415_scaled_centered_scores[ , 2], 
     xlab = as.expression(paste("PC1 (", round(pca_res_developing_seed_20131415_scaled_centered@R2[1], 4)*100, "%)", sep = "")), 
     ylab = as.expression(paste("PC2 (", round(pca_res_developing_seed_20131415_scaled_centered@R2[2], 4)*100, "%)", sep = "")),
     main = "Scores of PC1 and PC2 - Developing seed - Centered, scaled", 
     col = sample_list_developing_seed_131415$Batch, pch = 19)
legend("topright", 
       inset = c(-0.13,0), 
       fill = rainbow(length(levels(sample_list_developing_seed_131415$Batch))), 
       levels(sample_list_developing_seed_131415$Batch))
text(pca_res_developing_seed_20131415_scaled_centered_scores[ , 1], pca_res_developing_seed_20131415_scaled_centered_scores[ , 2], 
     labels = sample_list_developing_seed_131415$Batch, col = "black", pos = 4, cex=0.6)


#batch effect correction

#functions for batch correction

#===============================================================================
# Name   : RemoveFactors by ANOVA
# Author : Jan Lisec
# Date   : 2013-07-12
# Version: 0.1
# Aim    : Provide an interface to remove the influence of technical/biiological factors from experimental variation using an ANOVA model
# Mail   : <<<lisec@mpimp-golm.mpg.de>>>
#===============================================================================
RemoveFactors <- function(y=NULL, sam=NULL, facs=NULL, keep=NULL, output=c("y_norm","y_lm","anova_y","anova_y_norm","boxplot")[1]) {
  # y : data to normalize (numeric + in same order as sam
  # sam : dataframe containing the factors/numerical vars for ANOVA model
  # facs : all factors to be incorporated in the model in the desired order
  # keep : all factors to be retained in the normalized data
  tdf <- data.frame(y, sam[,facs]) # set up dataframe for anova
  y.lm <- lm(y ~ ., data=tdf) # set up anova model
  ce <- coef(y.lm) # these are the coefficients of the individual factors
  if (any(is.na(ce))) {
    ce[is.na(ce)] <- 0
    warning("Some coefficients were NA and had to be set to 0.\nYou probably have nested factors. Please check if ANOVA is appropriate.")
  }
  tm <- rep(ce[1], length(y)) # this is the total mean
  re <- residuals(y.lm) # residuals
  y.norm <- y # this maintaines any names y might have had
  fi <- is.finite(y) # this preserves the NAs in y by restoring only the finite values
  y.norm[fi] <- tm[fi] + re
  for (i in 1:length(keep)) {
    if (!is.factor(sam[,keep[i]])) {
      warning(paste("Can't keep numeric factors.", keep[i], "was removed."))
    } else {
      fac_eff <- sapply(levels(sam[,keep[i]]), function(x) {ifelse(is.na(ce[paste(keep[i],x,sep="")]),0,ce[paste(keep[i],x,sep="")])})
      y.norm[fi] <- y.norm[fi] + fac_eff[as.numeric(sam[,keep[i]])][fi] # add group mean values
    }
  }
  y.norm[fi] <- y.norm[fi] + (mean(y[fi]) - mean(y.norm[fi])) # correct for the offset
  if (output=="y_norm") return(y.norm)
  else if (output=="y_lm") return(y.lm)
  else if (output=="anova_y") return(anova(y.lm))
  else if (output=="anova_y_norm") {
    tdf <- data.frame(y.norm, sam[,facs])
    return(anova(lm(y.norm ~ ., data=tdf)))
  }
  else if (output=="boxplot") {
    par(mfrow=c(1,2))
    plot(y ~ interaction(sam[,keep]), xlab=paste(facs, collapse=", "))
    plot(y.norm ~ interaction(sam[,keep]), xlab=paste(keep, collapse=", "))
    par(mfrow=c(1,1))
    invisible(NULL)
  }
  invisible(NA)
}


#===============================================================================
# Name   : ANOVA-Normalization: Apply RemoveFactors function to values table
# Author : Heike Sprenger
# Date   : 2014-07-15
# Version: 0.1
#===============================================================================


# TODO: please test, if col/row renaming also works with 
# overlapped_analytes_sample_subset_log10_values = jkitest1_values_cast_select_log10[,-65]

# this applies source("RemoveFactors_function.R") with the following parameters:
#
# overlapped_analytes_sample_subset_log10_values : data to normalize (numeric + in same order as trial_factors)
# rows: samples, columns: analytes
# trial_factors : dataframe containing the factors/numerical vars for ANOVA model
# facs : all factors to be incorporated in the model in the desired order
# keep : all factors to be retained in the normalized data

func_normalize <- function(overlapped_analytes_sample_subset_log10_values, trial_factors,
                           facs = c("cultivar", "treatment", "sample_time", "SequenceID", "BatchID", "log10_AvgAnnotated"),
                           keep = c("cultivar", "treatment", "sample_time")) {
  normalized_values <- apply(overlapped_analytes_sample_subset_log10_values, 2, RemoveFactors,
                             sam = trial_factors, facs=facs, keep=keep)
  colnames(normalized_values) <- colnames(overlapped_analytes_sample_subset_log10_values) 
  rownames(normalized_values) <- rownames(overlapped_analytes_sample_subset_log10_values)
  print("normalized values dim:")
  print(dim(normalized_values))
  print("first 3 rows/cols:")
  print(normalized_values[1:3,1:3])
  return(normalized_values)
}


#batch correction/normalization
data_normalized_developing_seed_131415 = func_normalize(data_replaced_na_developing_seed_131415, 
                                                        trial_factors = sample_list_developing_seed_131415, 
                                                        facs = c("Treatment", "Timepoint", "Cultivar", "Batch", "Sequence"),
                                                        keep = c("Treatment", "Timepoint", "Cultivar"))

#check if negative values obtained
table(data_normalized_developing_seed_131415 < 0)   #negative values obtained

##if normalization yields negative values##

#log10 transformation - already done in earlier step (line 251)

#batch correction on log10-transformed data
data_transformed_normalized_developing_seed_131415 = func_normalize(data_log10_transformed_developing_seed_131415, 
                                                                    trial_factors = sample_list_developing_seed_131415,
                                                                    facs = c("Treatment", "Timepoint", "Cultivar", "Batch", "Sequence"),
                                                                    keep = c("Treatment", "Timepoint", "Cultivar"))

#transform values back
data_transformed_normalized_backtransformed_developing_seed_131415 = 10^(data_transformed_normalized_developing_seed_131415)  #FINAL DATA to use!!!

#check if negative values obtained
table(is.na(data_transformed_normalized_backtransformed_developing_seed_131415))    #no negative values


#PCA after batch correction

#log10 transformation for PCA
data_logtransform_norm_backtransform_logtransform = log10(data_transformed_normalized_backtransformed_developing_seed_131415)

#PCA - centered and scaled
pca_res_norm_developing_seed_20131415_scaled_centered = pca(data_logtransform_norm_backtransform_logtransform, 
                                                            method="ppca", 
                                                            nPcs=5, 
                                                            scale = "pareto", 
                                                            center = T)
pca_res_norm_developing_seed_20131415_scaled_centered
pca_res_norm_developing_seed_20131415_scaled_centered_scores = scores(pca_res_norm_developing_seed_20131415_scaled_centered)

#bar plot of PC variance with values
par(mar = c(5.1, 6.1, 4.1, 2.1))
text(barplot(pca_res_norm_developing_seed_20131415_scaled_centered@R2, 
             names.arg = c("PC1", "PC2", "PC3", "PC4", "PC5"), 
             main = "Probabilistic PCA - Developing seed - Normalized - Centered, scaled", 
             ylab = "Variance", 
             ylim = c(0, 0.3), 
             cex.axis = 1.5, 
             cex.names = 1.5, 
             cex.lab = 1.5, cex.main = 2), 
     0, round(pca_res_norm_developing_seed_20131415_scaled_centered@R2, 3), pos = 3, cex = 1.5)

#PCA score plot
palette(rainbow(length(levels(sample_list_developing_seed_131415$Batch))))
pairs(pca_res_norm_developing_seed_20131415_scaled_centered_scores, 
      col = sample_list_developing_seed_131415$Batch, 
      main = "PC score plots - Developing seed - Normalized - Centered, scaled")    #pairs

palette(rainbow(length(levels(sample_list_developing_seed_131415$Batch))))
par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = T)
plot(pca_res_norm_developing_seed_20131415_scaled_centered_scores[ , 1], 
     pca_res_norm_developing_seed_20131415_scaled_centered_scores[ , 2], 
     xlab = as.expression(paste("PC1 (", round(pca_res_norm_developing_seed_20131415_scaled_centered@R2[1], 4)*100, "%)", sep = "")),
     ylab = as.expression(paste("PC2 (", round(pca_res_norm_developing_seed_20131415_scaled_centered@R2[2], 4)*100, "%)", sep = "")),
     main = "Scores of PC1 and PC2 - Developing seed - Normalized - Centered, scaled", 
     col = sample_list_developing_seed_131415$Batch, pch = 19)
legend("topright", 
       inset = c(-0.13,0), 
       fill = rainbow(length(levels(sample_list_developing_seed_131415$Batch))), 
       levels(sample_list_developing_seed_131415$Batch))
text(pca_res_norm_developing_seed_20131415_scaled_centered_scores[ , 1], 
     pca_res_norm_developing_seed_20131415_scaled_centered_scores[ , 2], 
     labels = sample_list_developing_seed_131415$Batch, col = "black", pos = 4, cex=0.6)


#separate batch-corrected data set per year
data_batchnorm_developing_seed_2013 = data_transformed_normalized_backtransformed_developing_seed_131415[which(rownames(data_transformed_normalized_backtransformed_developing_seed_131415) %in% sample_list_developing_seed_2013$ID), ]
data_batchnorm_developing_seed_2014 = data_transformed_normalized_backtransformed_developing_seed_131415[which(rownames(data_transformed_normalized_backtransformed_developing_seed_131415) %in% sample_list_developing_seed_2014$ID), ]
data_batchnorm_developing_seed_2015 = data_transformed_normalized_backtransformed_developing_seed_131415[which(rownames(data_transformed_normalized_backtransformed_developing_seed_131415) %in% sample_list_developing_seed_2015$ID), ]


#heatmaps - visualize batch-corrected data per year
heatmap.2(data_batchnorm_developing_seed_2013, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          margins = c(15, 6), 
          main = "Developing seed - Batch-corrected data - 2013", 
          lmat = rbind(c(0, 3, 4), c(2, 1, 0)), 
          lwid = c(1.5, 4, 2), 
          xlab = "Metabolite", 
          ylab = "Sample")  

heatmap.2(data_batchnorm_developing_seed_2014, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          margins = c(15, 6), 
          main = "Developing seed - Batch-corrected data - 2014", 
          lmat = rbind(c(0, 3, 4), c(2, 1, 0)), 
          lwid = c(1.5, 4, 2), 
          xlab = "Metabolite", 
          ylab = "Sample")  

heatmap.2(data_batchnorm_developing_seed_2015, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          margins = c(15, 6), 
          main = "Developing seed - Batch-corrected data - 2015", 
          lmat = rbind(c(0, 3, 4), c(2, 1, 0)), 
          lwid = c(1.5, 4, 2), 
          xlab = "Metabolite", 
          ylab = "Sample")  


#outlier detection

#functions for handling of outliers


# code adapted from Lukasz Komsta's grubbs.test

# outlier value ("o" in grubbs.test) is stored in $outlier_value
# row name of the outlier ("G" in grubbs.test) is stored in $outlier_rowname

library(outliers)

#func_find_one_outlier
find_one_outlier <- function (x, output="none", opposite = FALSE)
{
  DNAME <- deparse(substitute(x))
  x <- sort(x[complete.cases(x)])
  n <- length(x)
  
  if (xor(((x[n] - mean(x)) < (mean(x) - x[1])), opposite)) {
    alt = paste("lowest value", x[1], "is an outlier")
    o <- x[1]
    d <- x[2:n]
  }
  else {
    alt = paste("highest value", x[n], "is an outlier")
    o <- x[n]
    d <- x[1:(n - 1)]
  }
  g <- abs(o - mean(x))/sd(x)
  u <- var(d)/var(x) * (n - 2)/(n - 1)
  pval <- 1 - pgrubbs(g, n, type = 10)
  method <- "Grubbs test for one outlier"
  
  RVAL <- list(statistic = c(G = g, U = u), alternative = alt,
               p.value = pval, method = method, data.name = DNAME,
               outlier_value = o, outlier_rowname = g)
  class(RVAL) <- "htest"
  if (output=="name")
    return(names(RVAL$outlier_value))
  else if (output=="p.value")
    return(RVAL$p.value)
  else
    return(RVAL)
}


# set.seed(1234)
# x = rnorm(10)
# print(x)
# 
# names(x) <- c("A","B","C","D","E","F","G","H","I","J")
# 
# result <- find_one_outlier(x)
# result
# names(result$outlier_value)


#func_replace_outlier
func_replace_outlier <- function(values, threshold, original_values, output="original"){
  idx_outlier_col = 1
  
  while(length(idx_outlier_col)>0){
    
    outlier_pvalue <- rep(NA, ncol(values))
    outlier_name <- rep(NA, ncol(values))
    
    for (i in 1:ncol(values)){
      outlier_pvalue[i] <- find_one_outlier(values[,i], output="p.value")
      outlier_name[i] <- find_one_outlier(values[,i], output="name")
    }
    
    idx_outlier_col <- which(outlier_pvalue < threshold)
    
    # replace outliers by NA
    for (i in idx_outlier_col){
      values[outlier_name[i], i] <- NA
    }
    
    for (i in idx_outlier_col){
      original_values[outlier_name[i], i] <- NA
    }
    
  }
  
  if (output=="original")
    return(original_values)
  else 
    return(values)
}


#func_hist_outlier
func_hist_outlier <- function(values, threshold){
  
  outlier_pvalue <- rep(NA, ncol(values))
  outlier_name <- rep(NA, ncol(values))
  
  for (i in 1:ncol(values)){
    outlier_pvalue[i] <- find_one_outlier(values[,i], output="p.value")
    outlier_name[i] <- find_one_outlier(values[,i], output="name")
  }
  
  idx_outlier_col <- which(outlier_pvalue < threshold)
  
  
  for (i in idx_outlier_col){
    hist(values[,i], breaks=20, col="grey",
         main=paste(colnames(values)[i], 
                    "\n p-value grubbs test: ",
                    outlier_pvalue[i]))
  }
  
}


#outlier test - threshold: 1e-4 --> accounts for condition-specific outliers

#histogram of outliers
pdf("outlier_hist_1e-4_developing_seed_2013.pdf")      #remove samples!! do this tomorrow!
func_hist_outlier(data_batchnorm_developing_seed_2013, threshold = 1e-4)  
dev.off()

pdf("outlier_hist_1e-4_developing_seed_2014.pdf")
func_hist_outlier(data_batchnorm_developing_seed_2014, threshold = 1e-4)  
dev.off()

pdf("outlier_hist_1e-4_developing_seed_2015.pdf")
func_hist_outlier(data_batchnorm_developing_seed_2015, threshold = 1e-4)  
dev.off()


#replace outliers with NA (missing value)
data_replaced_outliers_1e_neg4_developing_seed_2013 = func_replace_outlier(data_batchnorm_developing_seed_2013, threshold = 1e-4, data_batchnorm_developing_seed_2013)

data_replaced_outliers_1e_neg4_developing_seed_2014 = func_replace_outlier(data_batchnorm_developing_seed_2014, threshold = 1e-4, data_batchnorm_developing_seed_2014)

data_replaced_outliers_1e_neg4_developing_seed_2015 = func_replace_outlier(data_batchnorm_developing_seed_2015, threshold = 1e-4, data_batchnorm_developing_seed_2015)


#number of outliers
(sum(is.na(data_replaced_outliers_1e_neg4_developing_seed_2013))/(nrow(data_replaced_outliers_1e_neg4_developing_seed_2013)*ncol(data_replaced_outliers_1e_neg4_developing_seed_2013)))*100   #0.3% NAs

(sum(is.na(data_replaced_outliers_1e_neg4_developing_seed_2014))/(nrow(data_replaced_outliers_1e_neg4_developing_seed_2014)*ncol(data_replaced_outliers_1e_neg4_developing_seed_2014)))*100   #1.68% NAs

(sum(is.na(data_replaced_outliers_1e_neg4_developing_seed_2015))/(nrow(data_replaced_outliers_1e_neg4_developing_seed_2015)*ncol(data_replaced_outliers_1e_neg4_developing_seed_2015)))*100   #0.52% NAs


#heatmap - outliers replaced with missing values
heatmap.2(data_replaced_outliers_1e_neg4_developing_seed_2013, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          margins = c(15, 6), 
          main = "Developing seed - Outliers replaced, threshold: 1e-4 - 2013", 
          lmat = rbind(c(0, 3, 4), c(2, 1, 0)), 
          lwid = c(1.5, 4, 2), 
          xlab = "Metabolite", 
          ylab = "Sample", 
          na.color = "grey")      

heatmap.2(data_replaced_outliers_1e_neg4_developing_seed_2014, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          margins = c(15, 6), 
          main = "Developing seed - Outliers replaced, threshold: 1e-4 - 2014", 
          lmat = rbind(c(0, 3, 4), c(2, 1, 0)), 
          lwid = c(1.5, 4, 2), 
          xlab = "Metabolite", 
          ylab = "Sample", 
          na.color = "grey")

heatmap.2(data_replaced_outliers_1e_neg4_developing_seed_2015, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          margins = c(15, 6), 
          main = "Developing seed - Outliers replaced, threshold: 1e-4 - 2015", 
          lmat = rbind(c(0, 3, 4), c(2, 1, 0)), 
          lwid = c(1.5, 4, 2), 
          xlab = "Metabolite", 
          ylab = "Sample", 
          na.color = "grey")


#extract NA-replaced outlier data and manually inspect in Excel
#for condition-specific outliers, do no treat as outlier and keep the original value (batch-corrected) --> keep at least 3 replicates!

#NA-replaced data (to be inspected and replaced if necessary) - 2013
data_samples_replaced_outliers_1e_neg4_developing_seed_2013 = cbind(data_replaced_outliers_1e_neg4_developing_seed_2013, 
                                                                    as.character(sample_list_developing_seed_2013$Code))   #add column for sample code
colnames(data_samples_replaced_outliers_1e_neg4_developing_seed_2013)[ncol(data_samples_replaced_outliers_1e_neg4_developing_seed_2013)] = "Code"
data_samples_replaced_outliers_1e_neg4_developing_seed_2013 = data_samples_replaced_outliers_1e_neg4_developing_seed_2013[ , c(ncol(data_samples_replaced_outliers_1e_neg4_developing_seed_2013), 1:(ncol(data_samples_replaced_outliers_1e_neg4_developing_seed_2013)-1))]   #rearrange columns - place Code column before data values
write.table(data_samples_replaced_outliers_1e_neg4_developing_seed_2013, file = "data_replaced_outliers_1e-4_developing_seed_2013.txt", quote = FALSE, sep = "\t")     #export as txt file
#batch-corrected data (values to replace outliers/NAs)
data_samples_batchnorm_developing_seed_2013 = cbind(data_batchnorm_developing_seed_2013, 
                                                    as.character(sample_list_developing_seed_2013$Code))   #add column for sample code
colnames(data_samples_batchnorm_developing_seed_2013)[ncol(data_samples_batchnorm_developing_seed_2013)] = "Code"
data_samples_batchnorm_developing_seed_2013 = data_samples_batchnorm_developing_seed_2013[ , c(ncol(data_samples_batchnorm_developing_seed_2013), 1:(ncol(data_samples_batchnorm_developing_seed_2013)-1))]   #rearrange columns - place Code column before data values
write.table(data_samples_batchnorm_developing_seed_2013, file = "data_batchnorm_developing_seed_2013.txt", quote = FALSE, sep = "\t")

#NA-replaced data (to be inspected and replaced if necessary) - 2014
data_samples_replaced_outliers_1e_neg4_developing_seed_2014 = cbind(data_replaced_outliers_1e_neg4_developing_seed_2014, 
                                                                    as.character(sample_list_developing_seed_2014$Code))   #add column for sample code
colnames(data_samples_replaced_outliers_1e_neg4_developing_seed_2014)[ncol(data_samples_replaced_outliers_1e_neg4_developing_seed_2014)] = "Code"
data_samples_replaced_outliers_1e_neg4_developing_seed_2014 = data_samples_replaced_outliers_1e_neg4_developing_seed_2014[ , c(ncol(data_samples_replaced_outliers_1e_neg4_developing_seed_2014), 1:(ncol(data_samples_replaced_outliers_1e_neg4_developing_seed_2014)-1))]   #rearrange columns - place Code column before data values
write.table(data_samples_replaced_outliers_1e_neg4_developing_seed_2014, file = "data_replaced_outliers_1e-4_developing_seed_2014.txt", quote = FALSE, sep = "\t")     #export as txt file
#batch-corrected data (values to replace outliers/NAs)
data_samples_batchnorm_developing_seed_2014 = cbind(data_batchnorm_developing_seed_2014, 
                                                    as.character(sample_list_developing_seed_2014$Code))   #add column for sample code
colnames(data_samples_batchnorm_developing_seed_2014)[ncol(data_samples_batchnorm_developing_seed_2014)] = "Code"
data_samples_batchnorm_developing_seed_2014 = data_samples_batchnorm_developing_seed_2014[ , c(ncol(data_samples_batchnorm_developing_seed_2014), 1:(ncol(data_samples_batchnorm_developing_seed_2014)-1))]   #rearrange columns - place Code column before data values
write.table(data_samples_batchnorm_developing_seed_2014, file = "data_batchnorm_developing_seed_2014.txt", quote = FALSE, sep = "\t")

#NA-replaced data (to be inspected and replaced if necessary) - 2015
data_samples_replaced_outliers_1e_neg4_developing_seed_2015 = cbind(data_replaced_outliers_1e_neg4_developing_seed_2015, 
                                                                    as.character(sample_list_developing_seed_2015$Code))   #add column for sample code
colnames(data_samples_replaced_outliers_1e_neg4_developing_seed_2015)[ncol(data_samples_replaced_outliers_1e_neg4_developing_seed_2015)] = "Code"
data_samples_replaced_outliers_1e_neg4_developing_seed_2015 = data_samples_replaced_outliers_1e_neg4_developing_seed_2015[ , c(ncol(data_samples_replaced_outliers_1e_neg4_developing_seed_2015), 1:(ncol(data_samples_replaced_outliers_1e_neg4_developing_seed_2015)-1))]   #rearrange columns - place Code column before data values
write.table(data_samples_replaced_outliers_1e_neg4_developing_seed_2015, file = "data_replaced_outliers_1e-4_developing_seed_2015.txt", quote = FALSE, sep = "\t")     #export as txt file
#batch-corrected data (values to replace outliers/NAs)
data_samples_batchnorm_developing_seed_2015 = cbind(data_batchnorm_developing_seed_2015, 
                                                    as.character(sample_list_developing_seed_2015$Code))   #add column for sample code
colnames(data_samples_batchnorm_developing_seed_2015)[ncol(data_samples_batchnorm_developing_seed_2015)] = "Code"
data_samples_batchnorm_developing_seed_2015 = data_samples_batchnorm_developing_seed_2015[ , c(ncol(data_samples_batchnorm_developing_seed_2015), 1:(ncol(data_samples_batchnorm_developing_seed_2015)-1))]   #rearrange columns - place Code column before data values
write.table(data_samples_batchnorm_developing_seed_2015, file = "data_batchnorm_developing_seed_2015.txt", quote = FALSE, sep = "\t")


#load manually-curated outlier-replaced data

#2013 - no changes
data_outliers_replaced_manual_developing_seed_2013 = read.table("data_replaced_outliers_manual_1e-4_developing_seed_2013.txt", header = TRUE, sep = "\t", row.names = 1)
data_outliers_replaced_manual_developing_seed_2013 = data_outliers_replaced_manual_developing_seed_2013[ , -which(colnames(data_outliers_replaced_manual_developing_seed_2013) == "Code")]   #remove Code column
data_outliers_replaced_manual_developing_seed_2013 = as.matrix(data_outliers_replaced_manual_developing_seed_2013)  #convert data frame to matrix

#2014
data_outliers_replaced_manual_developing_seed_2014 = read.table("data_replaced_outliers_manual_1e-4_developing_seed_2014.txt", header = TRUE, sep = "\t", row.names = 1)
data_outliers_replaced_manual_developing_seed_2014 = data_outliers_replaced_manual_developing_seed_2014[ , -which(colnames(data_outliers_replaced_manual_developing_seed_2014) == "Code")]   #remove Code column
data_outliers_replaced_manual_developing_seed_2014 = as.matrix(data_outliers_replaced_manual_developing_seed_2014)  #convert data frame to matrix

#2015
data_outliers_replaced_manual_developing_seed_2015 = read.table("data_replaced_outliers_manual_1e-4_developing_seed_2015.txt", header = TRUE, sep = "\t", row.names = 1)
data_outliers_replaced_manual_developing_seed_2015 = data_outliers_replaced_manual_developing_seed_2015[ , -which(colnames(data_outliers_replaced_manual_developing_seed_2015) == "Code")]   #remove Code column
data_outliers_replaced_manual_developing_seed_2015 = as.matrix(data_outliers_replaced_manual_developing_seed_2015)  #convert data frame to matrix


#heatmap visualization - outliers replaced by R function and manual inspection

heatmap.2(data_outliers_replaced_manual_developing_seed_2013, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          margins = c(15, 6),
          main = "Developing seed - Outliers replaced_Final - 2013", 
          lmat = rbind(c(0, 3, 4), c(2, 1, 0)), 
          lwid = c(1.5, 4, 2), 
          xlab = "Metabolite", 
          ylab = "Sample", 
          labCol = reduced_metabolite_list_final_developing_seed_2013$Name, 
          na.color = "grey")

heatmap.2(data_outliers_replaced_manual_developing_seed_2014, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          margins = c(15, 6),
          main = "Developing seed - Outliers replaced_Final - 2014", 
          lmat = rbind(c(0, 3, 4), c(2, 1, 0)), 
          lwid = c(1.5, 4, 2), 
          xlab = "Metabolite", 
          ylab = "Sample",    
          labCol = reduced_metabolite_list_final_developing_seed_2014$Name,  
          na.color = "grey")

heatmap.2(data_outliers_replaced_manual_developing_seed_2015, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          margins = c(15, 6),
          main = "Developing seed - Outliers replaced_Final - 2015", 
          lmat = rbind(c(0, 3, 4), c(2, 1, 0)), 
          lwid = c(1.5, 4, 2), 
          xlab = "Metabolite", 
          ylab = "Sample",  
          labCol = reduced_metabolite_list_final_developing_seed_2015$Name,
          na.color = "grey")


##################

#combined analysis

#check if column names of data and sample lists to combine are the same
colnames(data_outliers_replaced_manual_developing_seed_2013) == colnames(data_outliers_replaced_manual_developing_seed_2014)  #OK
colnames(data_outliers_replaced_manual_developing_seed_2013) == colnames(data_outliers_replaced_manual_developing_seed_2015)  #OK
colnames(sample_list_developing_seed_2013_2) == colnames(sample_list_developing_seed_2014_2)    #OK
colnames(sample_list_developing_seed_2013_2) == colnames(sample_list_developing_seed_2015_2)    #OK

#combined 3-year data and sample list
data_combined_developing_seed = rbind(data_outliers_replaced_manual_developing_seed_2013, 
                                      data_outliers_replaced_manual_developing_seed_2014, 
                                      data_outliers_replaced_manual_developing_seed_2015) 
#combine 3-year sample list
sample_list_combined_developing_seed = rbind(sample_list_developing_seed_2013_2, 
                                             sample_list_developing_seed_2014_2, 
                                             sample_list_developing_seed_2015_2) 


#log2-median transformation

#func_log2_median_transform

#===============================================================================
# Name        : Log2-Median-Transformation
# Description : Calculate the median metabolite level for metabolite (over all samples) and
#               calculate log2-ratio this value from the respective metabolite levels in all samples
# Author      : adapted from Heike Sprenger (original: log10-median transformation)
# Date        : 2017-11-15
# Version     : 0.1
#===============================================================================


# input matrix with samples in columns and analytes in rows

func_log2_median_transform <- function(samples) {
  
  median_samples <- apply(samples, 1, median, na.rm=TRUE) # calculate median rowwise --> per analyte over all samples
  
  transformed <- data.frame(matrix(rep(NA, nrow(samples)*ncol(samples)), nrow=nrow(samples)))
  
  for (i  in 1:length(samples[1,])) {
    for (j in 1:length(samples[ ,1])) {
      transformed[j,i] <- log2(samples[j,i]/median_samples[j]) # calculate log2 of ratio of value/median
    }
  }
  colnames(transformed) <- colnames(samples)
  rownames(transformed) <- rownames(samples)
  return(transformed)
}


# input matrix with samples in rows and analytes in columns

func_log2_median_transform_2 <- function(samples) {
  
  samples <- t(samples)
  
  median_samples <- apply(samples, 1, median, na.rm=TRUE) # calculate median rowwise --> per analyte over all samples
  
  transformed <- data.frame(matrix(rep(NA, nrow(samples)*ncol(samples)), nrow=nrow(samples)))
  
  for (i  in 1:length(samples[1,])) {
    for (j in 1:length(samples[ ,1])) {
      transformed[j,i] <- log2(samples[j,i]/median_samples[j]) # calculate log2 of ratio of value/median
    }
  }
  colnames(transformed) <- colnames(samples)
  rownames(transformed) <- rownames(samples)
  
  transformed <- t(transformed)
  return(transformed)
}


data_combined_log2_median_transformed_developing_seed = func_log2_median_transform_2(data_combined_developing_seed)

#calculate mean per cultivar per timepoint
data_combined_log2_mean_developing_seed = aggregate(data_combined_log2_median_transformed_developing_seed, 
                                                    by = list(sample_list_combined_developing_seed$Treatment, 
                                                              sample_list_combined_developing_seed$Timepoint, 
                                                              sample_list_combined_developing_seed$Cultivar, 
                                                              sample_list_combined_developing_seed$Code), 
                                                    mean, na.rm = TRUE)      
colnames(data_combined_log2_mean_developing_seed)[1:4] = c("Treatment", "Timepoint", "Cultivar", "Code")   #rename column names

#mean data with code as rownames - for use in next steps
data_combined_log2_mean_developing_seed_2 = data_combined_log2_mean_developing_seed[ , 4:ncol(data_combined_log2_mean_developing_seed)]    #sample code and data only
rownames(data_combined_log2_mean_developing_seed_2) = data_combined_log2_mean_developing_seed_2$Code  #sample code as rowname
data_combined_log2_mean_developing_seed_2 = data_combined_log2_mean_developing_seed_2[ , -1]    #remove Code column
colnames(data_combined_log2_mean_developing_seed_2) = reduced_metabolite_list_final_developing_seed_2013$Name


#data structure to use for further analysis
data_combined_log2median_developing_seed = cbind(as.character(sample_list_combined_developing_seed$Cultivar), 
                                                 as.character(sample_list_combined_developing_seed$Trt.Time2), 
                                                 data_combined_log2_median_transformed_developing_seed)  #code, cultivar, timepoint, data
colnames(data_combined_log2median_developing_seed)[1:2] = c("Cultivar", "Timepoint")
data_combined_log2median_developing_seed = as.data.frame(data_combined_log2median_developing_seed)
data_combined_log2median_developing_seed[ , 3:ncol(data_combined_log2median_developing_seed)] = sapply(data_combined_log2median_developing_seed[ , 3:ncol(data_combined_log2median_developing_seed)], function(x) as.numeric(as.character(x))) #convert from factor to numeric


##comparison between control and stress##

#control vs. late stress

#log2-fold change - late stress vs. control
control_latestress_order_EGF = c("GNCLS", "GNLS", "GDCLS", "GDLS", "GACLS", "GALS")
data_combined_control_latestress_log2_mean_developing_seed = data_combined_log2_mean_developing_seed_2[control_latestress_order_EGF, ]  #control and late stress data in preferred order
data_combined_control_latestress_log2_mean_developing_seed = t(data_combined_control_latestress_log2_mean_developing_seed)
log2fc_combined_N22_control_latestress_developing_seed = data_combined_control_latestress_log2_mean_developing_seed[ , 2] - data_combined_control_latestress_log2_mean_developing_seed[ , 1]
log2fc_combined_Dular_control_latestress_developing_seed = data_combined_control_latestress_log2_mean_developing_seed[ , 4] - data_combined_control_latestress_log2_mean_developing_seed[ , 3]
log2fc_combined_Anjali_control_latestress_developing_seed = data_combined_control_latestress_log2_mean_developing_seed[ , 6] - data_combined_control_latestress_log2_mean_developing_seed[ , 5]
log2fc_combined_control_latestress_developing_seed = cbind(log2fc_combined_N22_control_latestress_developing_seed, 
                                                           log2fc_combined_Dular_control_latestress_developing_seed, 
                                                           log2fc_combined_Anjali_control_latestress_developing_seed)
colnames(log2fc_combined_control_latestress_developing_seed) = c(paste(colnames(data_combined_control_latestress_log2_mean_developing_seed)[2], "-", colnames(data_combined_control_latestress_log2_mean_developing_seed)[1]), 
                                                                 paste(colnames(data_combined_control_latestress_log2_mean_developing_seed)[4], "-", colnames(data_combined_control_latestress_log2_mean_developing_seed)[3]), 
                                                                 paste(colnames(data_combined_control_latestress_log2_mean_developing_seed)[6], "-", colnames(data_combined_control_latestress_log2_mean_developing_seed)[5]))

#heatmap - overview
heatmap.2(log2fc_combined_control_latestress_developing_seed, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none",
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5, 20),  
          main = "Developing seed - Combined - Log2 FC, Late Stress/Control", 
          cexCol = 1.5, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          labCol = c("N22", "Dular", "Anjali"), 
          srtCol = 0, 
          adjCol = c(0.5, 0.5),
          colsep = seq(0, ncol(log2fc_combined_control_latestress_developing_seed)), 
          rowsep = seq(0, nrow(log2fc_combined_control_latestress_developing_seed)), 
          sepcolor = "black", 
          sepwidth = c(0.01, 0.01))


#N22
data_combined_N22_log2median_developing_seed = data_combined_log2median_developing_seed[which(data_combined_log2median_developing_seed$Cultivar == "N22"), ]

#data distribution
shapiro_res_combined_N22_control_latestress_developing_seed = t(data.frame(lapply(data_combined_N22_log2median_developing_seed[which(data_combined_N22_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Late stress")), 3:ncol(data_combined_N22_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_control_latestress_developing_seed) = c("p.value")
shapiro_res_combined_N22_control_latestress_developing_seed = as.data.frame(shapiro_res_combined_N22_control_latestress_developing_seed)
shapiro_res_combined_N22_control_latestress_developing_seed$Distribution = ifelse(shapiro_res_combined_N22_control_latestress_developing_seed$p.value < 0.05, 
                                                                                  yes = "Non-normal", 
                                                                                  no = "Normal")
sum(shapiro_res_combined_N22_control_latestress_developing_seed$Distribution == "Normal")   #52 variables
sum(shapiro_res_combined_N22_control_latestress_developing_seed$Distribution == "Non-normal")   #15 variables

#Wilcoxon test
wilcox_res_N22_CLS_LS_developing_seed_combined = t(data.frame(lapply(data_combined_N22_log2median_developing_seed[ , 3:ncol(data_combined_N22_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_N22_log2median_developing_seed$Timepoint, subset = data_combined_N22_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Late stress"))$p.value)))
colnames(wilcox_res_N22_CLS_LS_developing_seed_combined) = c("p.value")
rownames(wilcox_res_N22_CLS_LS_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_N22_CLS_LS_developing_seed_combined = as.data.frame(wilcox_res_N22_CLS_LS_developing_seed_combined)
wilcox_res_N22_CLS_LS_developing_seed_combined$Significance = with(wilcox_res_N22_CLS_LS_developing_seed_combined,
                                                                   ifelse(p.value <= 0.001, "***", ifelse(
                                                                     p.value <= 0.01, "**", ifelse(
                                                                       p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_CLS_LS_developing_seed_combined$Significance == "*") #7 metabolites
sum(wilcox_res_N22_CLS_LS_developing_seed_combined$Significance == "**") #5 metabolites
sum(wilcox_res_N22_CLS_LS_developing_seed_combined$Significance == "***") #4 metabolites
sum(wilcox_res_N22_CLS_LS_developing_seed_combined$Significance != "ns") #16 metabolites
sum(wilcox_res_N22_CLS_LS_developing_seed_combined$Significance == "ns") #51 metabolites 


#Dular
data_combined_Dular_log2median_developing_seed = data_combined_log2median_developing_seed[which(data_combined_log2median_developing_seed$Cultivar == "Dular"), ]

#data distribution
shapiro_res_combined_Dular_control_latestress_developing_seed = t(data.frame(lapply(data_combined_Dular_log2median_developing_seed[which(data_combined_Dular_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Late stress")), 3:ncol(data_combined_Dular_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_control_latestress_developing_seed) = c("p.value")
shapiro_res_combined_Dular_control_latestress_developing_seed = as.data.frame(shapiro_res_combined_Dular_control_latestress_developing_seed)
shapiro_res_combined_Dular_control_latestress_developing_seed$Distribution = ifelse(shapiro_res_combined_Dular_control_latestress_developing_seed$p.value < 0.05, 
                                                                                    yes = "Non-normal", 
                                                                                    no = "Normal")
sum(shapiro_res_combined_Dular_control_latestress_developing_seed$Distribution == "Normal")   #44 variables
sum(shapiro_res_combined_Dular_control_latestress_developing_seed$Distribution == "Non-normal")   #23 variables

#Wilcoxon test
wilcox_res_Dular_CLS_LS_developing_seed_combined = t(data.frame(lapply(data_combined_Dular_log2median_developing_seed[ , 3:ncol(data_combined_Dular_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_developing_seed$Timepoint, subset = data_combined_Dular_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Late stress"))$p.value)))
colnames(wilcox_res_Dular_CLS_LS_developing_seed_combined) = c("p.value")
rownames(wilcox_res_Dular_CLS_LS_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_Dular_CLS_LS_developing_seed_combined = as.data.frame(wilcox_res_Dular_CLS_LS_developing_seed_combined)
wilcox_res_Dular_CLS_LS_developing_seed_combined$Significance = with(wilcox_res_Dular_CLS_LS_developing_seed_combined,
                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                       p.value <= 0.01, "**", ifelse(
                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_CLS_LS_developing_seed_combined$Significance == "*") #8 metabolites
sum(wilcox_res_Dular_CLS_LS_developing_seed_combined$Significance == "**") #4 metabolites
sum(wilcox_res_Dular_CLS_LS_developing_seed_combined$Significance == "***") #1 metabolites
sum(wilcox_res_Dular_CLS_LS_developing_seed_combined$Significance != "ns") #13 metabolites
sum(wilcox_res_Dular_CLS_LS_developing_seed_combined$Significance == "ns") #54 metabolites 


#Anjali
data_combined_Anjali_log2median_developing_seed = data_combined_log2median_developing_seed[which(data_combined_log2median_developing_seed$Cultivar == "Anjali"), ]

#data distribution
shapiro_res_combined_Anjali_control_latestress_developing_seed = t(data.frame(lapply(data_combined_Anjali_log2median_developing_seed[which(data_combined_Anjali_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Late stress")), 3:ncol(data_combined_Anjali_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_control_latestress_developing_seed) = c("p.value")
shapiro_res_combined_Anjali_control_latestress_developing_seed = as.data.frame(shapiro_res_combined_Anjali_control_latestress_developing_seed)
shapiro_res_combined_Anjali_control_latestress_developing_seed$Distribution = ifelse(shapiro_res_combined_Anjali_control_latestress_developing_seed$p.value < 0.05, 
                                                                                     yes = "Non-normal", 
                                                                                     no = "Normal")
sum(shapiro_res_combined_Anjali_control_latestress_developing_seed$Distribution == "Normal")   #51 variables
sum(shapiro_res_combined_Anjali_control_latestress_developing_seed$Distribution == "Non-normal")   #16 variables

#Wilcoxon test
wilcox_res_Anjali_CLS_LS_developing_seed_combined = t(data.frame(lapply(data_combined_Anjali_log2median_developing_seed[ , 3:ncol(data_combined_Anjali_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_developing_seed$Timepoint, subset = data_combined_Anjali_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Late stress"))$p.value)))
colnames(wilcox_res_Anjali_CLS_LS_developing_seed_combined) = c("p.value")
rownames(wilcox_res_Anjali_CLS_LS_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_Anjali_CLS_LS_developing_seed_combined = as.data.frame(wilcox_res_Anjali_CLS_LS_developing_seed_combined)
wilcox_res_Anjali_CLS_LS_developing_seed_combined$Significance = with(wilcox_res_Anjali_CLS_LS_developing_seed_combined,
                                                                      ifelse(p.value <= 0.001, "***", ifelse(
                                                                        p.value <= 0.01, "**", ifelse(
                                                                          p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_CLS_LS_developing_seed_combined$Significance == "*") #5 metabolites
sum(wilcox_res_Anjali_CLS_LS_developing_seed_combined$Significance == "**") #7 metabolites
sum(wilcox_res_Anjali_CLS_LS_developing_seed_combined$Significance == "***") #9 metabolites
sum(wilcox_res_Anjali_CLS_LS_developing_seed_combined$Significance != "ns") #21 metabolites
sum(wilcox_res_Anjali_CLS_LS_developing_seed_combined$Significance == "ns") #46 metabolites 


#metabolites with Wilcoxon test significant results

wilcox_sig_CLS_LS_developing_seed_combined = cbind(wilcox_res_N22_CLS_LS_developing_seed_combined, 
                                                   wilcox_res_Dular_CLS_LS_developing_seed_combined$Significance, 
                                                   wilcox_res_Anjali_CLS_LS_developing_seed_combined$Significance)
wilcox_sig_CLS_LS_developing_seed_combined = wilcox_sig_CLS_LS_developing_seed_combined[ , -1]
colnames(wilcox_sig_CLS_LS_developing_seed_combined)[1] = "wilcox_res_N22_CLS_LS_developing_seed_combined$Significance"
wilcox_sig_CLS_LS_developing_seed_combined = data.frame(lapply(wilcox_sig_CLS_LS_developing_seed_combined, 
                                                               function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_CLS_LS_developing_seed_combined) = rownames(wilcox_res_N22_CLS_LS_developing_seed_combined)
wilcox_sig_CLS_LS_developing_seed_combined$Nonsig = rowSums(wilcox_sig_CLS_LS_developing_seed_combined == " ")  #count number of nonsignificant
wilcox_sig_final_CLS_LS_developing_seed_combined = wilcox_sig_CLS_LS_developing_seed_combined[-which(wilcox_sig_CLS_LS_developing_seed_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars

wilcox_sig_log2fc_CLS_LS_developing_seed_combined = log2fc_combined_control_latestress_developing_seed[which(rownames(log2fc_combined_control_latestress_developing_seed) %in% rownames(wilcox_sig_final_CLS_LS_developing_seed_combined)), ]  #log2FC values of metabolites significant in at least one of the cultivars

#visualization
heatmap.2(wilcox_sig_log2fc_CLS_LS_developing_seed_combined, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)),        #1 value outside the scale
          margins = c(5, 20),  
          main = "Developing seed-Combined-Wilcoxon \nLog2 FC, Late stress/Control-Late stress", 
          cexCol = 1.5, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          labCol = c("N22", "Dular", "Anjali"), 
          srtCol = 0, 
          adjCol = c(0.5, 0.5), 
          cellnote = wilcox_sig_final_CLS_LS_developing_seed_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, ncol(wilcox_sig_log2fc_CLS_LS_developing_seed_combined))), 
          rowsep = c(seq(0, nrow(wilcox_sig_log2fc_CLS_LS_developing_seed_combined))), 
          sepwidth = c(0.01, 0.01))


#venn diagram - common among cultivars

#N22
wilcox_sig_N22_CLS_LS_developing_seed_combined = wilcox_res_N22_CLS_LS_developing_seed_combined[-which(wilcox_res_N22_CLS_LS_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_N22_CLS_LS_developing_seed_combined = as.data.frame(wilcox_sig_log2fc_CLS_LS_developing_seed_combined[match(rownames(wilcox_sig_N22_CLS_LS_developing_seed_combined), rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_N22_CLS_LS_developing_seed_combined) = colnames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined)[1]
wilcox_sig_log2fc_N22_CLS_LS_developing_seed_combined$Up.Down = ifelse(wilcox_sig_log2fc_N22_CLS_LS_developing_seed_combined$`GNLS - GNCLS` > 0, 
                                                                       yes = "Up", 
                                                                       no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_N22_CLS_LS_developing_seed_combined$Up.Down == "Up")    #8 metabolites
sum(wilcox_sig_log2fc_N22_CLS_LS_developing_seed_combined$Up.Down == "Down")  #8 metabolites

#Dular
wilcox_sig_Dular_CLS_LS_developing_seed_combined = wilcox_res_Dular_CLS_LS_developing_seed_combined[-which(wilcox_res_Dular_CLS_LS_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Dular_CLS_LS_developing_seed_combined = as.data.frame(wilcox_sig_log2fc_CLS_LS_developing_seed_combined[match(rownames(wilcox_sig_Dular_CLS_LS_developing_seed_combined), rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Dular_CLS_LS_developing_seed_combined) = colnames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined)[2]
wilcox_sig_log2fc_Dular_CLS_LS_developing_seed_combined$Up.Down = ifelse(wilcox_sig_log2fc_Dular_CLS_LS_developing_seed_combined$`GDLS - GDCLS` > 0, 
                                                                         yes = "Up", 
                                                                         no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Dular_CLS_LS_developing_seed_combined$Up.Down == "Up")    #9 metabolites
sum(wilcox_sig_log2fc_Dular_CLS_LS_developing_seed_combined$Up.Down == "Down")  #4 metabolites

#Anjali
wilcox_sig_Anjali_CLS_LS_developing_seed_combined = wilcox_res_Anjali_CLS_LS_developing_seed_combined[-which(wilcox_res_Anjali_CLS_LS_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Anjali_CLS_LS_developing_seed_combined = as.data.frame(wilcox_sig_log2fc_CLS_LS_developing_seed_combined[match(rownames(wilcox_sig_Anjali_CLS_LS_developing_seed_combined), rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Anjali_CLS_LS_developing_seed_combined) = colnames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined)[3]
wilcox_sig_log2fc_Anjali_CLS_LS_developing_seed_combined$Up.Down = ifelse(wilcox_sig_log2fc_Anjali_CLS_LS_developing_seed_combined$`GALS - GACLS` > 0, 
                                                                          yes = "Up", 
                                                                          no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Anjali_CLS_LS_developing_seed_combined$Up.Down == "Up")    #17 metabolites
sum(wilcox_sig_log2fc_Anjali_CLS_LS_developing_seed_combined$Up.Down == "Down")  #4 metabolites


#venn - increased - late stress/control-late stress
vennlist_up_wilcox_CLS_LS_developing_seed_combined = list(rownames(wilcox_sig_log2fc_N22_CLS_LS_developing_seed_combined)[wilcox_sig_log2fc_N22_CLS_LS_developing_seed_combined$Up.Down == "Up"],
                                                          rownames(wilcox_sig_log2fc_Dular_CLS_LS_developing_seed_combined)[wilcox_sig_log2fc_Dular_CLS_LS_developing_seed_combined$Up.Down == "Up"], 
                                                          rownames(wilcox_sig_log2fc_Anjali_CLS_LS_developing_seed_combined)[wilcox_sig_log2fc_Anjali_CLS_LS_developing_seed_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_CLS_LS_developing_seed_combined) = c("N22", "Dular", "Anjali")
venn_up_wilcox_CLS_LS_developing_seed_combined = venn(vennlist_up_wilcox_CLS_LS_developing_seed_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_late_stress_developing_seed_combined = venn.diagram(vennlist_up_wilcox_CLS_LS_developing_seed_combined, 
                                                            category.names = c(paste(names(vennlist_up_wilcox_CLS_LS_developing_seed_combined[1]), " (", length(vennlist_up_wilcox_CLS_LS_developing_seed_combined$N22), ")", sep = ""), 
                                                                               paste(names(vennlist_up_wilcox_CLS_LS_developing_seed_combined[2]), " (", length(vennlist_up_wilcox_CLS_LS_developing_seed_combined$Dular), ")", sep = ""), 
                                                                               paste(names(vennlist_up_wilcox_CLS_LS_developing_seed_combined[3]), " (", length(vennlist_up_wilcox_CLS_LS_developing_seed_combined$Anjali), ")", sep = "")),
                                                            filename = NULL, 
                                                            fill = c ("orchid", "yellow", "light blue"), 
                                                            cex = 4.5, 
                                                            fontfamily = "Arial", 
                                                            cat.cex = 2.5, 
                                                            cat.default.pos = c("text"),
                                                            cat.pos = c(12, -14, 0), 
                                                            cat.dist = c(0.1, 0.1, -0.07), 
                                                            cat.fontfamily = "Arial",
                                                            cat.fontface = "bold", 
                                                            euler.d = FALSE, 
                                                            scaled = FALSE, 
                                                            main = "A", 
                                                            main.pos = c(0.05, 0.95),
                                                            main.fontface = "bold", 
                                                            main.fontfamily = "Arial", 
                                                            main.cex = 3)   #with total number of metabolites per cultivar
#fig label A - as part of the final figure in the paper
intersect_wilcox_CLS_LS_venn_up_developing_seed_combined = attr(venn_up_wilcox_CLS_LS_developing_seed_combined, "intersections")  #list of unique and common metabolites
venn_up_wilcox_metabolites_CLS_LS_developing_seed_combined = t(ldply(intersect_wilcox_CLS_LS_venn_up_developing_seed_combined, rbind))  #convert from list to matrix
colnames(venn_up_wilcox_metabolites_CLS_LS_developing_seed_combined) = venn_up_wilcox_metabolites_CLS_LS_developing_seed_combined[1, ]
venn_up_wilcox_metabolites_CLS_LS_developing_seed_combined = venn_up_wilcox_metabolites_CLS_LS_developing_seed_combined[-1, ]


#venn - decreased - late stress/control-late stress
vennlist_down_wilcox_CLS_LS_developing_seed_combined = list(rownames(wilcox_sig_log2fc_N22_CLS_LS_developing_seed_combined)[wilcox_sig_log2fc_N22_CLS_LS_developing_seed_combined$Up.Down == "Down"],
                                                            rownames(wilcox_sig_log2fc_Dular_CLS_LS_developing_seed_combined)[wilcox_sig_log2fc_Dular_CLS_LS_developing_seed_combined$Up.Down == "Down"], 
                                                            rownames(wilcox_sig_log2fc_Anjali_CLS_LS_developing_seed_combined)[wilcox_sig_log2fc_Anjali_CLS_LS_developing_seed_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_CLS_LS_developing_seed_combined) = c("N22", "Dular", "Anjali")
venn_down_wilcox_CLS_LS_developing_seed_combined = venn(vennlist_down_wilcox_CLS_LS_developing_seed_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_late_stress_developing_seed_combined = venn.diagram(vennlist_down_wilcox_CLS_LS_developing_seed_combined, 
                                                              category.names = c(paste(names(vennlist_down_wilcox_CLS_LS_developing_seed_combined[1]), " (", length(vennlist_down_wilcox_CLS_LS_developing_seed_combined$N22), ")", sep = ""),
                                                                                 paste(names(vennlist_down_wilcox_CLS_LS_developing_seed_combined[2]), " (", length(vennlist_down_wilcox_CLS_LS_developing_seed_combined$Dular), ")", sep = ""), 
                                                                                 paste(names(vennlist_down_wilcox_CLS_LS_developing_seed_combined[3]), " (", length(vennlist_down_wilcox_CLS_LS_developing_seed_combined$Anjali), ")", sep = "")), 
                                                              filename = NULL, 
                                                              fill = c ("orchid", "yellow", "light blue"),
                                                              cex = 4.5, 
                                                              fontfamily = "Arial", 
                                                              cat.cex = 2.5, 
                                                              cat.default.pos = c("text"),
                                                              cat.pos = c(12, -14, 0), 
                                                              cat.dist = c(0.1, 0.1, -0.07), 
                                                              cat.fontfamily = "Arial",
                                                              cat.fontface = "bold", 
                                                              euler.d = FALSE, 
                                                              scaled = FALSE,
                                                              main = "B", 
                                                              main.pos = c(0.05, 0.95),
                                                              main.fontface = "bold", 
                                                              main.fontfamily = "Arial", 
                                                              main.cex = 3)   #with total number of metabolites per cultivar
#fig label B - as part of the final figure in the paper
intersect_wilcox_CLS_LS_venn_down_developing_seed_combined = attr(venn_down_wilcox_CLS_LS_developing_seed_combined, "intersections")  #list of unique and common metabolites
venn_down_wilcox_metabolites_CLS_LS_developing_seed_combined = t(ldply(intersect_wilcox_CLS_LS_venn_down_developing_seed_combined, rbind))  #convert from list to matrix
colnames(venn_down_wilcox_metabolites_CLS_LS_developing_seed_combined) = venn_down_wilcox_metabolites_CLS_LS_developing_seed_combined[1, ]
venn_down_wilcox_metabolites_CLS_LS_developing_seed_combined = venn_down_wilcox_metabolites_CLS_LS_developing_seed_combined[-1, ]


#export tables
write.table(venn_up_wilcox_metabolites_CLS_LS_developing_seed_combined, file = "wilcox_sig_up_metabolites_developing_seed_late_stress_combined.txt", sep = "\t", quote = F)
write.table(venn_down_wilcox_metabolites_CLS_LS_developing_seed_combined, file = "wilcox_sig_down_metabolites_developing_seed_late_stress_combined.txt", sep = "\t", quote = F)


##################

#final venn diagram figure for paper

png("venn_up_down_developing_seed_combined.png", width = 16*300, height = 8*300, res = 300)
grid.arrange(gTree(children = venn_up_late_stress_developing_seed_combined), 
             gTree(children = venn_down_late_stress_developing_seed_combined),
             ncol = 2)
dev.off()

##################

##comparison among controls##

#subset data and sample list for controls only - control corresponding to late stress
data_combined_control_developing_seed = data_combined_developing_seed[which(sample_list_combined_developing_seed$Trt.Time2 == "Control.Late stress"), ]
sample_list_combined_control_developing_seed = sample_list_combined_developing_seed[which(sample_list_combined_developing_seed$Trt.Time2 == "Control.Late stress"), ]

#log2-median transformation
data_combined_control_log2_median_transformed_developing_seed = func_log2_median_transform_2(data_combined_control_developing_seed)

#mean per condition
data_combined_control_log2_median_mean_developing_seed = aggregate(data_combined_control_log2_median_transformed_developing_seed, 
                                                                   by = list(sample_list_combined_control_developing_seed$Treatment, 
                                                                             sample_list_combined_control_developing_seed$Timepoint, 
                                                                             sample_list_combined_control_developing_seed$Cultivar, 
                                                                             sample_list_combined_control_developing_seed$Code),
                                                                   mean, na.rm = TRUE)
colnames(data_combined_control_log2_median_mean_developing_seed)[1:4] = c("Treatment", "Timepoint", "Cultivar", "Code")   #rename columns

#mean data with code as rownames
data_combined_control_log2_median_mean_developing_seed_2 = data_combined_control_log2_median_mean_developing_seed[ , -c(1, 2, 3)]
rownames(data_combined_control_log2_median_mean_developing_seed_2) = data_combined_control_log2_median_mean_developing_seed_2$Code
data_combined_control_log2_median_mean_developing_seed_2 = data_combined_control_log2_median_mean_developing_seed_2[ , -1]

#transform data - metabolites in rows, samples in columns
data_combined_control_log2_median_mean_developing_seed_2 = t(data_combined_control_log2_median_mean_developing_seed_2)

#order based on tolerance - N22, Dular, Anjali
data_combined_control_log2_median_mean_developing_seed_2 = data_combined_control_log2_median_mean_developing_seed_2[, c(3, 2, 1)]

#heatmap - overview
heatmap.2(data_combined_control_log2_median_mean_developing_seed_2, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(4, 20),  main = "Developing seed - Combined_Controls-Late stress", 
          cexCol = 1.5, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5),
          labCol = c("N22", "Dular", "Anjali"), 
          srtCol = 0, 
          adjCol = c(0.5, 0.5),
          colsep = seq(0, ncol(data_combined_control_log2_median_mean_developing_seed_2)), 
          rowsep = seq(0, nrow(data_combined_control_log2_median_mean_developing_seed_2)), 
          sepcolor = "black", 
          sepwidth = c(0.01, 0.01))


#log2-fold difference
log2fc_combined_N22_Dular_control_developing_seed = data_combined_control_log2_median_mean_developing_seed_2[ , 1] - data_combined_control_log2_median_mean_developing_seed_2[ , 2]   #N22/Dular
log2fc_combined_N22_Anjali_control_developing_seed = data_combined_control_log2_median_mean_developing_seed_2[ , 1] - data_combined_control_log2_median_mean_developing_seed_2[ , 3]    #N22/Anjali
log2fc_combined_Dular_Anjali_control_developing_seed = data_combined_control_log2_median_mean_developing_seed_2[ , 2] - data_combined_control_log2_median_mean_developing_seed_2[ , 3]    #Dular/Anjali
log2fc_combined_control_developing_seed = cbind(log2fc_combined_N22_Dular_control_developing_seed, 
                                                log2fc_combined_N22_Anjali_control_developing_seed, 
                                                log2fc_combined_Dular_Anjali_control_developing_seed)   #combine in one data set
colnames(log2fc_combined_control_developing_seed) = c("N22/Dular", "N22/Anjali", "Dular/Anjali")
rownames(log2fc_combined_control_developing_seed) = reduced_metabolite_list_final_developing_seed_2013$Name

#heatmap - overview
heatmap.2(log2fc_combined_control_developing_seed,
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(4, 20), 
          main = "Developing seed - Combined_Controls - Log2 FD", 
          cexCol = 1.5, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          srtCol = 0, 
          adjCol = c(0.5, 0.5),
          colsep = seq(0, ncol(log2fc_combined_control_developing_seed)), 
          rowsep = seq(0, nrow(log2fc_combined_control_developing_seed)), 
          sepcolor = "black", 
          sepwidth = c(0.01, 0.01))


#data structure to use for further analysis
data_combined_log2median_control_developing_seed = cbind(as.character(sample_list_combined_control_developing_seed$Trt.Time2), 
                                                         as.character(sample_list_combined_control_developing_seed$Cultivar), 
                                                         data_combined_control_log2_median_transformed_developing_seed)
colnames(data_combined_log2median_control_developing_seed)[1:2] = c("Trt.Timepoint", "Cultivar")
data_combined_log2median_control_developing_seed = as.data.frame(data_combined_log2median_control_developing_seed)
data_combined_log2median_control_developing_seed[ , 3:ncol(data_combined_log2median_control_developing_seed)] = sapply(data_combined_log2median_control_developing_seed[3:ncol(data_combined_log2median_control_developing_seed)], function(x) as.numeric(as.character(x)))    #convert from factor to numeric


#N22/Dular

#data distribution
shapiro_res_combined_control_N22_Dular_developing_seed = t(as.data.frame(lapply(data_combined_log2median_control_developing_seed[which(data_combined_log2median_control_developing_seed$Cultivar %in% c("N22", "Dular")), 3:ncol(data_combined_log2median_control_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_control_N22_Dular_developing_seed) = "p.value"
shapiro_res_combined_control_N22_Dular_developing_seed = as.data.frame(shapiro_res_combined_control_N22_Dular_developing_seed)
shapiro_res_combined_control_N22_Dular_developing_seed$Distribution = ifelse(shapiro_res_combined_control_N22_Dular_developing_seed$p.value < 0.05, 
                                                                             yes = "Non-normal", 
                                                                             no = "Normal")
sum(shapiro_res_combined_control_N22_Dular_developing_seed$Distribution == "Non-normal")   #20 metabolites
sum(shapiro_res_combined_control_N22_Dular_developing_seed$Distribution == "Normal")       #47 metabolites

#Wilcoxon test
wilcox_res_combined_control_N22_Dular_developing_seed = t(data.frame(lapply(data_combined_log2median_control_developing_seed[ , 3:ncol(data_combined_log2median_control_developing_seed)], function(x) wilcox.test(x ~ data_combined_log2median_control_developing_seed$Cultivar, subset = data_combined_log2median_control_developing_seed$Cultivar %in% c("N22", "Dular"))$p.value)))
colnames(wilcox_res_combined_control_N22_Dular_developing_seed) = "p.value"
wilcox_res_combined_control_N22_Dular_developing_seed = as.data.frame(wilcox_res_combined_control_N22_Dular_developing_seed)
rownames(wilcox_res_combined_control_N22_Dular_developing_seed) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_combined_control_N22_Dular_developing_seed$Significance = with(wilcox_res_combined_control_N22_Dular_developing_seed,
                                                                          ifelse(p.value <= 0.001, "***",
                                                                                 ifelse(p.value <= 0.01, "**",
                                                                                        ifelse(p.value < 0.05, "*", "ns"))))
sum(wilcox_res_combined_control_N22_Dular_developing_seed$Significance == "*") #4 metabolites
sum(wilcox_res_combined_control_N22_Dular_developing_seed$Significance == "**") #8 metabolites
sum(wilcox_res_combined_control_N22_Dular_developing_seed$Significance == "***") #7 metabolites
sum(wilcox_res_combined_control_N22_Dular_developing_seed$Significance != "ns") #19 metabolites
sum(wilcox_res_combined_control_N22_Dular_developing_seed$Significance == "ns") #48 metabolites


#N22/Anjali

#data distribution
shapiro_res_combined_control_N22_Anjali_developing_seed = t(as.data.frame(lapply(data_combined_log2median_control_developing_seed[which(data_combined_log2median_control_developing_seed$Cultivar %in% c("N22", "Anjali")), 3:ncol(data_combined_log2median_control_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_control_N22_Anjali_developing_seed) = "p.value"
shapiro_res_combined_control_N22_Anjali_developing_seed = as.data.frame(shapiro_res_combined_control_N22_Anjali_developing_seed)
shapiro_res_combined_control_N22_Anjali_developing_seed$Distribution = ifelse(shapiro_res_combined_control_N22_Anjali_developing_seed$p.value < 0.05, 
                                                                              yes = "Non-normal",
                                                                              no = "Normal")
sum(shapiro_res_combined_control_N22_Anjali_developing_seed$Distribution == "Non-normal")   #25 metabolites
sum(shapiro_res_combined_control_N22_Anjali_developing_seed$Distribution == "Normal")       #42 metabolites

#Wilcoxon test
wilcox_res_combined_control_N22_Anjali_developing_seed = t(data.frame(lapply(data_combined_log2median_control_developing_seed[ , 3:ncol(data_combined_log2median_control_developing_seed)], function(x) wilcox.test(x ~ data_combined_log2median_control_developing_seed$Cultivar, subset = data_combined_log2median_control_developing_seed$Cultivar %in% c("N22", "Anjali"))$p.value)))
colnames(wilcox_res_combined_control_N22_Anjali_developing_seed) = "p.value"
wilcox_res_combined_control_N22_Anjali_developing_seed = as.data.frame(wilcox_res_combined_control_N22_Anjali_developing_seed)
rownames(wilcox_res_combined_control_N22_Anjali_developing_seed) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_combined_control_N22_Anjali_developing_seed$Significance = with(wilcox_res_combined_control_N22_Anjali_developing_seed,
                                                                           ifelse(p.value <= 0.001, "***",
                                                                                  ifelse(p.value <= 0.01, "**",
                                                                                         ifelse(p.value < 0.05, "*", "ns"))))
sum(wilcox_res_combined_control_N22_Anjali_developing_seed$Significance == "*") #3 metabolites
sum(wilcox_res_combined_control_N22_Anjali_developing_seed$Significance == "**") #9 metabolites
sum(wilcox_res_combined_control_N22_Anjali_developing_seed$Significance == "***") #14 metabolites
sum(wilcox_res_combined_control_N22_Anjali_developing_seed$Significance != "ns") #26 metabolites
sum(wilcox_res_combined_control_N22_Anjali_developing_seed$Significance == "ns") #41 metabolites


#Dular/Anjali

#data distribution
shapiro_res_combined_control_Dular_Anjali_developing_seed = t(as.data.frame(lapply(data_combined_log2median_control_developing_seed[which(data_combined_log2median_control_developing_seed$Cultivar %in% c("Dular", "Anjali")), 3:ncol(data_combined_log2median_control_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_control_Dular_Anjali_developing_seed) = "p.value"
shapiro_res_combined_control_Dular_Anjali_developing_seed = as.data.frame(shapiro_res_combined_control_Dular_Anjali_developing_seed)
shapiro_res_combined_control_Dular_Anjali_developing_seed$Distribution = ifelse(shapiro_res_combined_control_Dular_Anjali_developing_seed$p.value < 0.05, 
                                                                                yes = "Non-normal",
                                                                                no = "Normal")
sum(shapiro_res_combined_control_Dular_Anjali_developing_seed$Distribution == "Non-normal")   #17 metabolites
sum(shapiro_res_combined_control_Dular_Anjali_developing_seed$Distribution == "Normal")       #50 metabolites

#Wilcoxon test
wilcox_res_combined_control_Dular_Anjali_developing_seed = t(data.frame(lapply(data_combined_log2median_control_developing_seed[ , 3:ncol(data_combined_log2median_control_developing_seed)], function(x) wilcox.test(x ~ data_combined_log2median_control_developing_seed$Cultivar, subset = data_combined_log2median_control_developing_seed$Cultivar %in% c("Dular", "Anjali"))$p.value)))
colnames(wilcox_res_combined_control_Dular_Anjali_developing_seed) = "p.value"
wilcox_res_combined_control_Dular_Anjali_developing_seed = as.data.frame(wilcox_res_combined_control_Dular_Anjali_developing_seed)
rownames(wilcox_res_combined_control_Dular_Anjali_developing_seed) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_combined_control_Dular_Anjali_developing_seed$Significance = with(wilcox_res_combined_control_Dular_Anjali_developing_seed,
                                                                             ifelse(p.value <= 0.001, "***",
                                                                                    ifelse(p.value <= 0.01, "**",
                                                                                           ifelse(p.value < 0.05, "*", "ns"))))
sum(wilcox_res_combined_control_Dular_Anjali_developing_seed$Significance == "*") #12 metabolites
sum(wilcox_res_combined_control_Dular_Anjali_developing_seed$Significance == "**") #5 metabolites
sum(wilcox_res_combined_control_Dular_Anjali_developing_seed$Significance == "***") #11 metabolites
sum(wilcox_res_combined_control_Dular_Anjali_developing_seed$Significance != "ns") #28 metabolites
sum(wilcox_res_combined_control_Dular_Anjali_developing_seed$Significance == "ns") #39 metabolites


#metabolites with significant Wilcoxon test results
wilcox_sig_combined_control_developing_seed = cbind(wilcox_res_combined_control_N22_Dular_developing_seed,
                                                    wilcox_res_combined_control_N22_Anjali_developing_seed$Significance, 
                                                    wilcox_res_combined_control_Dular_Anjali_developing_seed$Significance)
wilcox_sig_combined_control_developing_seed = wilcox_sig_combined_control_developing_seed[ , -1]
colnames(wilcox_sig_combined_control_developing_seed)[1] = "wilcox_res_combined_control_N22_Dular_developing_seed$Significance"
wilcox_sig_combined_control_developing_seed = data.frame(lapply(wilcox_sig_combined_control_developing_seed, 
                                                                function(x) gsub("ns", " ", x)))
rownames(wilcox_sig_combined_control_developing_seed) = rownames(wilcox_res_combined_control_N22_Dular_developing_seed)
wilcox_sig_combined_control_developing_seed$Nonsig = rowSums(wilcox_sig_combined_control_developing_seed == " ")    #number of nonsignificant results
wilcox_sig_final_combined_control_developing_seed = wilcox_sig_combined_control_developing_seed[-which(wilcox_sig_combined_control_developing_seed$Nonsig == "3"), ]   #metabolites significant in at least one of the comparisons

#log2FD values of metabolites significant in at least one of the comparisons
wilcox_sig_log2fc_combined_control_developing_seed = log2fc_combined_control_developing_seed[which(rownames(log2fc_combined_control_developing_seed) %in% rownames(wilcox_sig_final_combined_control_developing_seed)), ]


#visualization - heatmap with row colors to categorize metabolites

#make another file as data frame - to add column for metabolite class for heatmap category
wilcox_sig_log2fc_combined_control_developing_seed_2 = wilcox_sig_log2fc_combined_control_developing_seed

#convert to data frame to be able to add another column
wilcox_sig_log2fc_combined_control_developing_seed_2 = as.data.frame(wilcox_sig_log2fc_combined_control_developing_seed_2)

#add column for metabolite class to clasify in heatmap
wilcox_sig_log2fc_combined_control_developing_seed_2$Metabolite.Class = reduced_metabolite_list_final_developing_seed_2013[match(rownames(wilcox_sig_log2fc_combined_control_developing_seed_2), reduced_metabolite_list_final_developing_seed_2013$Name), "Class"]

#convert/save as factor
metabolite_class_developing_seed = factor(wilcox_sig_log2fc_combined_control_developing_seed_2$Metabolite.Class)

#color palette for metabolite class category - from flowering spikelet (11 categories)- to have the same color per metabolite class - except color for Aromatic & Phenylpropanoid groups - not in developing seed data set
#same for FL and EGF
color_heatmap_developing_seed = c("#8DD3C7", "#FFFFB3", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5")

#heatmap - Developing seed-Combined-Wilcoxon-Log2 fold difference, Controls (Control-Late stress) - with metabolite class categories
png("wilcox_sig_log2fc_control_developing_seed_combined.png", width = 10*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_combined_control_developing_seed, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5.8, 51.5), 
          cexCol = 1.2, 
          cexRow = 1.2, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 11), 
          lwid = c(1.1, 5), 
          srtCol = 90, 
          adjCol = c(0.9, 0.5), 
          cellnote = wilcox_sig_final_combined_control_developing_seed, 
          notecol = "black", 
          notecex = 2, 
          sepcolor = "black", 
          colsep = c(seq(0, ncol(wilcox_sig_log2fc_combined_control_developing_seed))), 
          rowsep = c(seq(0, nrow(wilcox_sig_log2fc_combined_control_developing_seed))), 
          sepwidth = c(0.01, 0.01),
          RowSideColors = color_heatmap_developing_seed[metabolite_class_developing_seed])
legend(x = -0.06, y = 1.03, xpd = TRUE, legend = expression(paste("Log" ["2"], "-fold difference")), bty = "n", cex = 0.8)
legend(x = 0.5, y = 1, legend = levels(metabolite_class_developing_seed), fill = color_heatmap_developing_seed[factor(levels(metabolite_class_developing_seed))], bty = "n", cex = 0.8)
legend(x = 0.145, y = 1.125, xpd = TRUE, legend = c("B"), bty = "n", cex = 1.5, text.font = 2)
dev.off()


#venn diagram - common and specific metabolites among comparisons

#N22/Dular
wilcox_sig_combined_control_N22_Dular_developing_seed = wilcox_res_combined_control_N22_Dular_developing_seed[-which(wilcox_res_combined_control_N22_Dular_developing_seed$Significance == "ns"), ]    #significantly different metabolites
wilcox_sig_log2fc_combined_control_N22_Dular_developing_seed = as.data.frame(wilcox_sig_log2fc_combined_control_developing_seed[match(rownames(wilcox_sig_combined_control_N22_Dular_developing_seed), rownames(wilcox_sig_log2fc_combined_control_developing_seed)), 1])   #log2fc values of significantly different metabolites
colnames(wilcox_sig_log2fc_combined_control_N22_Dular_developing_seed) = colnames(wilcox_sig_log2fc_combined_control_developing_seed)[1]
wilcox_sig_log2fc_combined_control_N22_Dular_developing_seed$Up.Down = ifelse(wilcox_sig_log2fc_combined_control_N22_Dular_developing_seed$`N22/Dular` > 0, 
                                                                              yes = "Up", 
                                                                              no = "Down") #categorize if FD is up/down
sum(wilcox_sig_log2fc_combined_control_N22_Dular_developing_seed$Up.Down == "Up")    #9 metabolites
sum(wilcox_sig_log2fc_combined_control_N22_Dular_developing_seed$Up.Down == "Down")  #10 metabolites

#N22/Anjali
wilcox_sig_combined_control_N22_Anjali_developing_seed = wilcox_res_combined_control_N22_Anjali_developing_seed[-which(wilcox_res_combined_control_N22_Anjali_developing_seed$Significance == "ns"), ]    #significantly different metabolites
wilcox_sig_log2fc_combined_control_N22_Anjali_developing_seed = as.data.frame(wilcox_sig_log2fc_combined_control_developing_seed[match(rownames(wilcox_sig_combined_control_N22_Anjali_developing_seed), rownames(wilcox_sig_log2fc_combined_control_developing_seed)), 2])   #log2fc values of significantly different metabolites
colnames(wilcox_sig_log2fc_combined_control_N22_Anjali_developing_seed) = colnames(wilcox_sig_log2fc_combined_control_developing_seed)[2]
wilcox_sig_log2fc_combined_control_N22_Anjali_developing_seed$Up.Down = ifelse(wilcox_sig_log2fc_combined_control_N22_Anjali_developing_seed$`N22/Anjali` > 0, 
                                                                               yes = "Up", 
                                                                               no = "Down") #categorize if FD is up/down
sum(wilcox_sig_log2fc_combined_control_N22_Anjali_developing_seed$Up.Down == "Up")    #25 metabolites
sum(wilcox_sig_log2fc_combined_control_N22_Anjali_developing_seed$Up.Down == "Down")  #1 metabolite

#Dular/Anjali
wilcox_sig_combined_control_Dular_Anjali_developing_seed = wilcox_res_combined_control_Dular_Anjali_developing_seed[-which(wilcox_res_combined_control_Dular_Anjali_developing_seed$Significance == "ns"), ]    #significantly different metabolites
wilcox_sig_log2fc_combined_control_Dular_Anjali_developing_seed = as.data.frame(wilcox_sig_log2fc_combined_control_developing_seed[match(rownames(wilcox_sig_combined_control_Dular_Anjali_developing_seed), rownames(wilcox_sig_log2fc_combined_control_developing_seed)), 3])   #log2fc values of significantly different metabolites
colnames(wilcox_sig_log2fc_combined_control_Dular_Anjali_developing_seed) = colnames(wilcox_sig_log2fc_combined_control_developing_seed)[3]
wilcox_sig_log2fc_combined_control_Dular_Anjali_developing_seed$Up.Down = ifelse(wilcox_sig_log2fc_combined_control_Dular_Anjali_developing_seed$`Dular/Anjali` > 0, 
                                                                                 yes = "Up", 
                                                                                 no = "Down") #categorize if FD is up/down
sum(wilcox_sig_log2fc_combined_control_Dular_Anjali_developing_seed$Up.Down == "Up")    #27 metabolites
sum(wilcox_sig_log2fc_combined_control_Dular_Anjali_developing_seed$Up.Down == "Down")  #1 metabolite


#three-way venn - increased & decreased combined in one venn
vennlist_up_down_wilcox_combined_control_developing_seed = list(rownames(wilcox_sig_log2fc_combined_control_N22_Dular_developing_seed),
                                                                rownames(wilcox_sig_log2fc_combined_control_N22_Anjali_developing_seed),
                                                                rownames(wilcox_sig_log2fc_combined_control_Dular_Anjali_developing_seed))
names(vennlist_up_down_wilcox_combined_control_developing_seed) = colnames(wilcox_sig_log2fc_combined_control_developing_seed)
venn_up_down_wilcox_combined_control_developing_seed = venn(vennlist_up_down_wilcox_combined_control_developing_seed)    #plots venn diagram and lists results
grid.newpage()
venn_up_down_combined_control_developing_seed = venn.diagram(vennlist_up_down_wilcox_combined_control_developing_seed, 
                                                             category.names = c(paste(names(vennlist_up_down_wilcox_combined_control_developing_seed[1]), " (", length(vennlist_up_down_wilcox_combined_control_developing_seed$`N22/Dular`), ")", sep = ""), 
                                                                                paste(names(vennlist_up_down_wilcox_combined_control_developing_seed[2]), " (", length(vennlist_up_down_wilcox_combined_control_developing_seed$`N22/Anjali`), ")", sep = ""), 
                                                                                paste(names(vennlist_up_down_wilcox_combined_control_developing_seed[3]), " (", length(vennlist_up_down_wilcox_combined_control_developing_seed$`Dular/Anjali`), ")", sep = "")),
                                                             filename = NULL, 
                                                             fill = c ("orchid", "yellow", "light blue"),
                                                             cex = 4,
                                                             fontfamily = "Arial", 
                                                             cat.cex = 2, 
                                                             cat.default.pos = c("text"), 
                                                             cat.pos = c(21, -21, 0),
                                                             cat.dist = c(0.1, 0.1, -0.07), 
                                                             cat.fontfamily = "Arial", 
                                                             cat.fontface = "bold",
                                                             euler.d = FALSE, 
                                                             scaled = FALSE, 
                                                             main = "A", 
                                                             main.pos = c(0.05, 0.95),
                                                             main.fontface = "bold", 
                                                             main.fontfamily = "Arial", 
                                                             main.cex = 3)   #with total number of metabolites per comparison
#fig label A - as part of the final figure in the paper
intersect_venn_up_down_wilcox_combined_control_developing_seed = attr(venn_up_down_wilcox_combined_control_developing_seed, "intersections")  #list of unique and common metabolites
venn_up_down_wilcox_metabolites_combined_control_developing_seed = t(ldply(intersect_venn_up_down_wilcox_combined_control_developing_seed, rbind))  #convert from list to matrix
colnames(venn_up_down_wilcox_metabolites_combined_control_developing_seed) = venn_up_down_wilcox_metabolites_combined_control_developing_seed[1, ]
venn_up_down_wilcox_metabolites_combined_control_developing_seed = venn_up_down_wilcox_metabolites_combined_control_developing_seed[-1, ]


##################

#final venn diagram figure for paper - 2-column figure so that the size will be the same as in flag leaf figure
png("venn_up&down_combined_control_developing_seed.png", width = 16*300, height = 8*300, res = 300)
grid.arrange(gTree(children = venn_up_down_combined_control_developing_seed), 
             gTree(children = NULL), 
             ncol = 2)
dev.off()


##################

##correlation analysis##

#correlation between change in yield and metabolite changes
#correlation between change in proportion of chalky grains (>50% chalk content) and metabolite changes

#yield
yield_HxD_data_combined = read.table("yield_HxD.txt", header = TRUE, sep = "\t")

#chalkiness
chalk_50_75_HxD_data_combined = read.table("chalkiness_50_75_HxD.txt", header = TRUE, sep = "\t")


#EGF stage

yield_HxD_EGF = yield_HxD_data_combined[which(yield_HxD_data_combined$Stage == "EGF"), ]    #yield under EGF only
chalk_50_75_HxD_EGF = chalk_50_75_HxD_data_combined[which(chalk_50_75_HxD_data_combined$Stage == "EGF"), ]    #grains with 50-75% chalk under EGF only
chalk_50_75_HxD_EGF$Chalky_grains_50_75 = chalk_50_75_HxD_EGF$X50.75..chalk + chalk_50_75_HxD_EGF$X.75..chalk     #sum of grains with 50-75% chalk


#yield and chalky grains from the three experiments in one data frame 

#make sure to have the same ID per row first
yield_HxD_EGF$Year == chalk_50_75_HxD_EGF$Year
yield_HxD_EGF$Cultivar == chalk_50_75_HxD_EGF$Cultivar
yield_HxD_EGF$Treatment == chalk_50_75_HxD_EGF$Treatment

#combine in one data set
yield_chalk_HxD_EGF = cbind(yield_HxD_EGF, chalk_50_75_HxD_EGF$Chalky_grains_50_75)
colnames(yield_chalk_HxD_EGF)[ncol(yield_chalk_HxD_EGF)] = "Chalky_50_75"   #edit column name

yield_chalk_HxD_EGF_mean = aggregate(yield_chalk_HxD_EGF[ , c("Yield", "Chalky_50_75")], 
                                     by = list(yield_chalk_HxD_EGF$Year, 
                                               yield_chalk_HxD_EGF$Cultivar, 
                                               yield_chalk_HxD_EGF$Treatment), 
                                     mean)     #mean per year, cultivar, treatment
colnames(yield_chalk_HxD_EGF_mean)[1:3] = c("Year", "Cultivar", "Treatment")
rownames(yield_chalk_HxD_EGF_mean) = interaction(yield_chalk_HxD_EGF_mean$Year, yield_chalk_HxD_EGF_mean$Cultivar, yield_chalk_HxD_EGF_mean$Treatment)    #assign rownames as Year.Cultivar.Treatment
yield_chalk_HxD_EGF_mean = yield_chalk_HxD_EGF_mean[ , 4:ncol(yield_chalk_HxD_EGF_mean)]     #data and sample ID as rownames
yield_chalk_HxD_EGF_mean = t(yield_chalk_HxD_EGF_mean)        #transform - sample ID in columns
yield_chalk_change_HxD_EGF = ((yield_chalk_HxD_EGF_mean[ , 10:ncol(yield_chalk_HxD_EGF_mean)] - yield_chalk_HxD_EGF_mean[ , 1:9])/yield_chalk_HxD_EGF_mean[ , 1:9])*100     #change in yield and chalkiness (stress-control)
colnames(yield_chalk_change_HxD_EGF) = gsub("\\.", " ", colnames(yield_chalk_change_HxD_EGF))
colnames(yield_chalk_change_HxD_EGF) = gsub(" Stress", "", colnames(yield_chalk_change_HxD_EGF))
rownames(yield_chalk_change_HxD_EGF) = c("Change in yield", "Change in chalky_50-75%")


#metabolite log2-fold change from the three experiments in one data frame

#calculate first log2-fold change from the means of median-normalized data per experiment

#mean of log2-median transformed values

#2013
data_log2_median_transformed_developing_seed_2013 = func_log2_median_transform_2(data_outliers_replaced_manual_developing_seed_2013)
data_log2_mean_developing_seed_2013 = aggregate(data_log2_median_transformed_developing_seed_2013, 
                                                by = list(sample_list_developing_seed_2013$Treatment, 
                                                          sample_list_developing_seed_2013$Timepoint, 
                                                          sample_list_developing_seed_2013$Cultivar, 
                                                          sample_list_developing_seed_2013$Code), 
                                                mean, na.rm = TRUE)
colnames(data_log2_mean_developing_seed_2013)[1:4] = c("Treatment", "Timepoint", "Cultivar", "Code")

data_log2_mean_developing_seed_2013_2 = data_log2_mean_developing_seed_2013[ , 4:ncol(data_log2_mean_developing_seed_2013)]   #data and code as sample ID
rownames(data_log2_mean_developing_seed_2013_2) = data_log2_mean_developing_seed_2013_2$Code
data_log2_mean_developing_seed_2013_2 = data_log2_mean_developing_seed_2013_2[ , -1]
colnames(data_log2_mean_developing_seed_2013_2) = reduced_metabolite_list_final_developing_seed_2013$Name

#2014
data_log2_median_transformed_developing_seed_2014 = func_log2_median_transform_2(data_outliers_replaced_manual_developing_seed_2014)
data_log2_mean_developing_seed_2014 = aggregate(data_log2_median_transformed_developing_seed_2014, 
                                                by = list(sample_list_developing_seed_2014$Treatment, 
                                                          sample_list_developing_seed_2014$Timepoint, 
                                                          sample_list_developing_seed_2014$Cultivar, 
                                                          sample_list_developing_seed_2014$Code), 
                                                mean, na.rm = TRUE)
colnames(data_log2_mean_developing_seed_2014)[1:4] = c("Treatment", "Timepoint", "Cultivar", "Code")

data_log2_mean_developing_seed_2014_2 = data_log2_mean_developing_seed_2014[ , 4:ncol(data_log2_mean_developing_seed_2014)]   #data and code as sample ID
rownames(data_log2_mean_developing_seed_2014_2) = data_log2_mean_developing_seed_2014_2$Code
data_log2_mean_developing_seed_2014_2 = data_log2_mean_developing_seed_2014_2[ , -1]
colnames(data_log2_mean_developing_seed_2014_2) = reduced_metabolite_list_final_developing_seed_2014$Name

#2015
data_log2_median_transformed_developing_seed_2015 = func_log2_median_transform_2(data_outliers_replaced_manual_developing_seed_2015)
data_log2_mean_developing_seed_2015 = aggregate(data_log2_median_transformed_developing_seed_2015, 
                                                by = list(sample_list_developing_seed_2015$Treatment, 
                                                          sample_list_developing_seed_2015$Timepoint, 
                                                          sample_list_developing_seed_2015$Cultivar, 
                                                          sample_list_developing_seed_2015$Code), 
                                                mean, na.rm = TRUE)
colnames(data_log2_mean_developing_seed_2015)[1:4] = c("Treatment", "Timepoint", "Cultivar", "Code")

data_log2_mean_developing_seed_2015_2 = data_log2_mean_developing_seed_2015[ , 4:ncol(data_log2_mean_developing_seed_2015)]   #data and code as sample ID
rownames(data_log2_mean_developing_seed_2015_2) = data_log2_mean_developing_seed_2015_2$Code
data_log2_mean_developing_seed_2015_2 = data_log2_mean_developing_seed_2015_2[ , -1]
colnames(data_log2_mean_developing_seed_2015_2) = reduced_metabolite_list_final_developing_seed_2015$Name


#calculate log2-fold change (stress/control)

control_stress_order_EGF = c("GNCLS", "GNLS", "GDCLS", "GDLS", "GACLS", "GALS")  #preferred order


#2013
data_control_stress_log2_mean_developing_seed_2013 = data_log2_mean_developing_seed_2013_2[control_stress_order_EGF, ]  #data with control and stress values only
data_control_stress_log2_mean_developing_seed_2013 = t(data_control_stress_log2_mean_developing_seed_2013)  #samples in columns

log2fc_N22_control_stress_developing_seed_2013 = data_control_stress_log2_mean_developing_seed_2013[ , 2] - data_control_stress_log2_mean_developing_seed_2013[ , 1]    #late stress/control-late stress
log2fc_Dular_control_stress_developing_seed_2013 = data_control_stress_log2_mean_developing_seed_2013[ , 4] - data_control_stress_log2_mean_developing_seed_2013[ , 3]     #late stress/control-late stress
log2fc_Anjali_control_stress_developing_seed_2013 = data_control_stress_log2_mean_developing_seed_2013[ , 6] - data_control_stress_log2_mean_developing_seed_2013[ , 5]     #late stress/control-late stress
log2fc_control_stress_developing_seed_2013 = cbind(log2fc_N22_control_stress_developing_seed_2013, 
                                                   log2fc_Dular_control_stress_developing_seed_2013, 
                                                   log2fc_Anjali_control_stress_developing_seed_2013)   #combine in one data set
#rename columns
colnames(log2fc_control_stress_developing_seed_2013) = c(paste(colnames(data_control_stress_log2_mean_developing_seed_2013)[2], "-", colnames(data_control_stress_log2_mean_developing_seed_2013)[1]),
                                                         paste(colnames(data_control_stress_log2_mean_developing_seed_2013)[4], "-", colnames(data_control_stress_log2_mean_developing_seed_2013)[3]),
                                                         paste(colnames(data_control_stress_log2_mean_developing_seed_2013)[6], "-", colnames(data_control_stress_log2_mean_developing_seed_2013)[5]))

#2014
data_control_stress_log2_mean_developing_seed_2014 = data_log2_mean_developing_seed_2014_2[control_stress_order_EGF, ]  #data with control and stress values only
data_control_stress_log2_mean_developing_seed_2014 = t(data_control_stress_log2_mean_developing_seed_2014)  #samples in columns

log2fc_N22_control_stress_developing_seed_2014 = data_control_stress_log2_mean_developing_seed_2014[ , 2] - data_control_stress_log2_mean_developing_seed_2014[ , 1]    #late stress/control-late stress
log2fc_Dular_control_stress_developing_seed_2014 = data_control_stress_log2_mean_developing_seed_2014[ , 4] - data_control_stress_log2_mean_developing_seed_2014[ , 3]     #late stress/control-late stress
log2fc_Anjali_control_stress_developing_seed_2014 = data_control_stress_log2_mean_developing_seed_2014[ , 6] - data_control_stress_log2_mean_developing_seed_2014[ , 5]     #late stress/control-late stress
log2fc_control_stress_developing_seed_2014 = cbind(log2fc_N22_control_stress_developing_seed_2014, 
                                                   log2fc_Dular_control_stress_developing_seed_2014, 
                                                   log2fc_Anjali_control_stress_developing_seed_2014)   #combine in one data set
#rename columns
colnames(log2fc_control_stress_developing_seed_2014) = c(paste(colnames(data_control_stress_log2_mean_developing_seed_2014)[2], "-", colnames(data_control_stress_log2_mean_developing_seed_2014)[1]),
                                                         paste(colnames(data_control_stress_log2_mean_developing_seed_2014)[4], "-", colnames(data_control_stress_log2_mean_developing_seed_2014)[3]),
                                                         paste(colnames(data_control_stress_log2_mean_developing_seed_2014)[6], "-", colnames(data_control_stress_log2_mean_developing_seed_2014)[5]))

#2015
data_control_stress_log2_mean_developing_seed_2015 = data_log2_mean_developing_seed_2015_2[control_stress_order_EGF, ]  #data with control and stress values only
data_control_stress_log2_mean_developing_seed_2015 = t(data_control_stress_log2_mean_developing_seed_2015)  #samples in columns

log2fc_N22_control_stress_developing_seed_2015 = data_control_stress_log2_mean_developing_seed_2015[ , 2] - data_control_stress_log2_mean_developing_seed_2015[ , 1]    #late stress/control-late stress
log2fc_Dular_control_stress_developing_seed_2015 = data_control_stress_log2_mean_developing_seed_2015[ , 4] - data_control_stress_log2_mean_developing_seed_2015[ , 3]     #late stress/control-late stress
log2fc_Anjali_control_stress_developing_seed_2015 = data_control_stress_log2_mean_developing_seed_2015[ , 6] - data_control_stress_log2_mean_developing_seed_2015[ , 5]     #late stress/control-late stress
log2fc_control_stress_developing_seed_2015 = cbind(log2fc_N22_control_stress_developing_seed_2015, 
                                                   log2fc_Dular_control_stress_developing_seed_2015, 
                                                   log2fc_Anjali_control_stress_developing_seed_2015)   #combine in one data set
#rename columns
colnames(log2fc_control_stress_developing_seed_2015) = c(paste(colnames(data_control_stress_log2_mean_developing_seed_2015)[2], "-", colnames(data_control_stress_log2_mean_developing_seed_2015)[1]),
                                                         paste(colnames(data_control_stress_log2_mean_developing_seed_2015)[4], "-", colnames(data_control_stress_log2_mean_developing_seed_2015)[3]),
                                                         paste(colnames(data_control_stress_log2_mean_developing_seed_2015)[6], "-", colnames(data_control_stress_log2_mean_developing_seed_2015)[5]))


#check first if rownames are the same before combining data
rownames(log2fc_control_stress_developing_seed_2013) == rownames(log2fc_control_stress_developing_seed_2014)
rownames(log2fc_control_stress_developing_seed_2013) == rownames(log2fc_control_stress_developing_seed_2015)

#log2-fold change of metabolites
log2fc_developing_seed_control_latestress_20131415 = cbind(log2fc_control_stress_developing_seed_2013[ , c(3, 2, 1)], 
                                                           log2fc_control_stress_developing_seed_2014[ , c(3, 2, 1)], 
                                                           log2fc_control_stress_developing_seed_2015[ , c(3, 2, 1)])    #log2fc of late stress relative to control - 3-yr data in one data table - order: Anjali, Dular, N22
colnames(log2fc_developing_seed_control_latestress_20131415) = c("2013 Anjali", "2013 Dular", "2013 N22", 
                                                                 "2014 Anjali", "2014 Dular", "2014 N22",
                                                                 "2015 Anjali", "2015 Dular", "2015 N22")
log2fc_developing_seed_control_latestress_20131415 = log2fc_developing_seed_control_latestress_20131415[ , match(colnames(yield_chalk_change_HxD_EGF), colnames(log2fc_developing_seed_control_latestress_20131415))]     #reorder columns to match column order of yield and chalkiness data


#check first if rownames are the same before combining data
colnames(yield_chalk_change_HxD_EGF) == colnames(log2fc_developing_seed_control_latestress_20131415)

#yield, chalkiness, and log2-fold change in one data table
yield_chalk_log2fc_control_latestress_developing_seed_HxD = rbind(yield_chalk_change_HxD_EGF, 
                                                                  log2fc_developing_seed_control_latestress_20131415)
yield_chalk_log2fc_control_latestress_developing_seed_HxD = t(yield_chalk_log2fc_control_latestress_developing_seed_HxD)
yield_chalk_log2fc_control_latestress_developing_seed_HxD = as.data.frame(yield_chalk_log2fc_control_latestress_developing_seed_HxD)    #convert to data frame


#data distribution

#func_shapiro_test
func_shapiro_test = function(data){
  shapiro_test = apply(data, 2, shapiro.test)
  shapiro_res = sapply(shapiro_test, `[`, c("statistic","p.value"))
  return(data.frame(t(shapiro_res)))
}

shapiro_yield_chalk_log2fc_developing_seed_latestress_HxD = func_shapiro_test(yield_chalk_log2fc_control_latestress_developing_seed_HxD)
shapiro_yield_chalk_log2fc_developing_seed_latestress_HxD$Distribution = ifelse(shapiro_yield_chalk_log2fc_developing_seed_latestress_HxD$p.value < 0.05, 
                                                                                yes = "Non-normal", 
                                                                                no = "Normal")
sum(shapiro_yield_chalk_log2fc_developing_seed_latestress_HxD$Distribution == "Normal")   #67 variables
sum(shapiro_yield_chalk_log2fc_developing_seed_latestress_HxD$Distribution == "Non-normal")   #2 variables


#correlation test - since not all variables have normally distributed data, use spearman method

#correlation coefficient
cor_res_spearman_rho_developing_seed_latestress_HxD = cor(yield_chalk_log2fc_control_latestress_developing_seed_HxD, 
                                                          method = "spearman")

#test for significance

#yield
yield_cor_res_spearman_pval_developing_seed_latestress_HXD = t(data.frame(lapply(yield_chalk_log2fc_control_latestress_developing_seed_HxD[3:ncol(yield_chalk_log2fc_control_latestress_developing_seed_HxD)], function(x) cor.test(yield_chalk_log2fc_control_latestress_developing_seed_HxD$`Change in yield`, x, method = "spearman")$p.value)))    #significance of correlation
yield_cor_res_spearman_pval_developing_seed_latestress_HXD = as.data.frame(yield_cor_res_spearman_pval_developing_seed_latestress_HXD)
colnames(yield_cor_res_spearman_pval_developing_seed_latestress_HXD) = "p.value"
rownames(yield_cor_res_spearman_pval_developing_seed_latestress_HXD) = reduced_metabolite_list_final_developing_seed_2013$Name
yield_cor_res_spearman_pval_developing_seed_latestress_HXD$Significance = with(yield_cor_res_spearman_pval_developing_seed_latestress_HXD,
                                                                               ifelse(p.value <= 0.001, "***", ifelse(
                                                                                 p.value <= 0.01, "**", ifelse(
                                                                                   p.value < 0.05, "*", "ns"))))
length(which(yield_cor_res_spearman_pval_developing_seed_latestress_HXD$p.value < 0.05))   #5 metabolites
yield_cor_sig_metabolites_developing_seed_latestress_HxD = rownames(yield_cor_res_spearman_pval_developing_seed_latestress_HXD)[yield_cor_res_spearman_pval_developing_seed_latestress_HXD$Significance != "ns"]   #metabolites with significant correlation with change in yield
yield_sig_cor_rho_spearman_developing_seed_latestress_HxD = data.frame(cor_res_spearman_rho_developing_seed_latestress_HxD[1, which(colnames(cor_res_spearman_rho_developing_seed_latestress_HxD) %in% yield_cor_sig_metabolites_developing_seed_latestress_HxD)])     #metabolites with significant correlation with change in yield and rho values
colnames(yield_sig_cor_rho_spearman_developing_seed_latestress_HxD) = "rho value"
yield_sig_cor_rho_spearman_developing_seed_latestress_HxD = as.matrix(yield_sig_cor_rho_spearman_developing_seed_latestress_HxD)    #convert to matrix


#proportion of chalky grains
chalk_50_75_cor_res_spearman_pval_developing_seed_latestress_HxD = t(data.frame(lapply(yield_chalk_log2fc_control_latestress_developing_seed_HxD[3:ncol(yield_chalk_log2fc_control_latestress_developing_seed_HxD)], function(x) cor.test(yield_chalk_log2fc_control_latestress_developing_seed_HxD$`Change in chalky_50-75%`, x, method = "spearman")$p.value)))    #significance of correlation
chalk_50_75_cor_res_spearman_pval_developing_seed_latestress_HxD = as.data.frame(chalk_50_75_cor_res_spearman_pval_developing_seed_latestress_HxD)
colnames(chalk_50_75_cor_res_spearman_pval_developing_seed_latestress_HxD) = "p.value"
rownames(chalk_50_75_cor_res_spearman_pval_developing_seed_latestress_HxD) = reduced_metabolite_list_final_developing_seed_2013$Name
chalk_50_75_cor_res_spearman_pval_developing_seed_latestress_HxD$Significance = with(chalk_50_75_cor_res_spearman_pval_developing_seed_latestress_HxD,
                                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                                       p.value <= 0.01, "**", ifelse(
                                                                                         p.value < 0.05, "*", "ns"))))
length(which(chalk_50_75_cor_res_spearman_pval_developing_seed_latestress_HxD$p.value < 0.05))   #0 metabolite


#correlation heatmap - yield

#matrix with 2 columns of the same data - dummy - to be able to make a 1-column heatmap
yield_sig_cor_rho_spearman_developing_seed_latestress_HxD_2 = cbind(yield_sig_cor_rho_spearman_developing_seed_latestress_HxD,
                                                                    yield_sig_cor_rho_spearman_developing_seed_latestress_HxD)
yield_sig_cor_rho_spearman_developing_seed_latestress_HxD_3 = yield_sig_cor_rho_spearman_developing_seed_latestress_HxD_2[order(rownames(yield_sig_cor_rho_spearman_developing_seed_latestress_HxD_2)), ] #metabolites sorted alphabetically

#create a new matrix with dummy data to make equal-sized heatmap with the others - 14 metabolites in total
dummy_data_1 = matrix(c(0, 0), nrow = 9, ncol = 2, dimnames = list(c(rep("Z", 9))))
yield_sig_cor_rho_spearman_developing_seed_latestress_HxD_4 = rbind(yield_sig_cor_rho_spearman_developing_seed_latestress_HxD_3,
                                                                    dummy_data_1)

#heatmap 
png("cor_yield_late_stress_developing_seed_HxD.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(yield_sig_cor_rho_spearman_developing_seed_latestress_HxD_4, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          breaks = c(seq(-1, -0.01, length.out = 300), seq(0.01, 1, length.out = 300)), 
          symbreaks = TRUE, 
          col = bluered(599), 
          colsep = c(0, ncol(yield_sig_cor_rho_spearman_developing_seed_latestress_HxD_4)), 
          rowsep = c(seq(0, nrow(yield_sig_cor_rho_spearman_developing_seed_latestress_HxD_4))), 
          sepcolor = "black", 
          sepwidth = c(0.001, 0.005), 
          margins = c(20, 35), 
          cexRow = 1, 
          labCol = "", 
          lhei = c(0.45, 5), 
          lwid = c(1.5, 5), 
          density.info = "none",
          key.title = NA, 
          key.xlab = NA,
          key.xtickfun = function() {
            cex <- par("cex")*par("cex.axis")
            side <- 1
            line <- 0
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("-1", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("0", side=side, at=0.5, adj=0.5,
                  line=line, cex=cex, col=col, font=font)
            mtext("1", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          })
legend(x = -0.083, y = 1.06, xpd = TRUE, legend = expression(paste("Correlation coefficient")), bty = "n", cex = 0.8)
legend(0.15, 1.12, xpd = TRUE, legend = c("F"), bty = "n", cex = 1.5)
dev.off()
#crop dummy data at the end of the figure


##################

##correlation analysis##

#correlation between metabolite levels under control and change in yield and proportion of chalky grains

#control data only
control_order_EGF = c("GNCLS", "GDCLS", "GACLS")

#2013
data_control_log2_mean_developing_seed_2013 = data_log2_mean_developing_seed_2013_2[control_order_EGF, ]    #data with controls (corresponding to late stress) only
data_control_log2_mean_developing_seed_2013 = t(data_control_log2_mean_developing_seed_2013)  #samples in columns

#2014
data_control_log2_mean_developing_seed_2014 = data_log2_mean_developing_seed_2014_2[control_order_EGF, ]    #data with controls (corresponding to late stress) only
data_control_log2_mean_developing_seed_2014 = t(data_control_log2_mean_developing_seed_2014)  #samples in columns

#2015
data_control_log2_mean_developing_seed_2015 = data_log2_mean_developing_seed_2015_2[control_order_EGF, ]    #data with controls (corresponding to late stress) only
data_control_log2_mean_developing_seed_2015 = t(data_control_log2_mean_developing_seed_2015)  #samples in columns


#check first if rownames are the same before combining data
rownames(data_control_log2_mean_developing_seed_2013) == rownames(data_control_log2_mean_developing_seed_2014)
rownames(data_control_log2_mean_developing_seed_2013) == rownames(data_control_log2_mean_developing_seed_2015)

#combine into one data set
control_developing_seed_20131415 = cbind(data_control_log2_mean_developing_seed_2013[, c(3, 2, 1)], 
                                         data_control_log2_mean_developing_seed_2014[, c(3, 2, 1)], 
                                         data_control_log2_mean_developing_seed_2015[, c(3, 2, 1)])    #control values - 3-yr data in one data table - order: Anjali, Dular, N22
colnames(control_developing_seed_20131415) = c("2013 Anjali", "2013 Dular", "2013 N22", 
                                               "2014 Anjali", "2014 Dular", "2014 N22", 
                                               "2015 Anjali", "2015 Dular", "2015 N22")
control_developing_seed_20131415 = control_developing_seed_20131415[ , match(colnames(yield_chalk_change_HxD_EGF), colnames(control_developing_seed_20131415))]   #reorder columns to match column order of yield and chalkiness data


#combine agronomic & metabolomic data

#check first if rownames are the same before combining data
colnames(yield_chalk_change_HxD_EGF) == colnames(control_developing_seed_20131415)

#yield, chalkiness, and control values in one data table
yield_chalk_control_developing_seed_HxD = rbind(yield_chalk_change_HxD_EGF, control_developing_seed_20131415)
yield_chalk_control_developing_seed_HxD = t(yield_chalk_control_developing_seed_HxD)
yield_chalk_control_developing_seed_HxD = as.data.frame(yield_chalk_control_developing_seed_HxD)


#data distribution
shapiro_yield_chalk_control_developing_seed_HxD = func_shapiro_test(yield_chalk_control_developing_seed_HxD)
shapiro_yield_chalk_control_developing_seed_HxD$Distribution = ifelse(shapiro_yield_chalk_control_developing_seed_HxD$p.value < 0.05, 
                                                                      yes = "Non-normal",
                                                                      no = "Normal")
sum(shapiro_yield_chalk_control_developing_seed_HxD$Distribution == "Normal")   #63 variables
sum(shapiro_yield_chalk_control_developing_seed_HxD$Distribution == "Non-normal")   #6 variables


#correlation test - since not all variables have normally distributed data, use spearman method

#correlation coefficient
cor_res_spearman_rho_developing_seed_control_HxD = cor(yield_chalk_control_developing_seed_HxD, method = "spearman")

#test for significance

#yield
yield_cor_res_spearman_pval_developing_seed_control_HxD = t(data.frame(lapply(yield_chalk_control_developing_seed_HxD[3:ncol(yield_chalk_control_developing_seed_HxD)], function(x) cor.test(yield_chalk_control_developing_seed_HxD$`Change in yield`, x, method = "spearman")$p.value)))    #significance of correlation
yield_cor_res_spearman_pval_developing_seed_control_HxD = as.data.frame(yield_cor_res_spearman_pval_developing_seed_control_HxD)
colnames(yield_cor_res_spearman_pval_developing_seed_control_HxD) = "p.value"
rownames(yield_cor_res_spearman_pval_developing_seed_control_HxD) = reduced_metabolite_list_final_developing_seed_2013$Name
yield_cor_res_spearman_pval_developing_seed_control_HxD$Significance = with(yield_cor_res_spearman_pval_developing_seed_control_HxD,
                                                                            ifelse(p.value <= 0.001, "***", ifelse(
                                                                              p.value <= 0.01, "**", ifelse(
                                                                                p.value < 0.05, "*", "ns"))))
length(which(yield_cor_res_spearman_pval_developing_seed_control_HxD$p.value < 0.05))   #2 metabolites
yield_cor_sig_metabolites_developing_seed_control_HxD = rownames(yield_cor_res_spearman_pval_developing_seed_control_HxD)[yield_cor_res_spearman_pval_developing_seed_control_HxD$Significance != "ns"]   #metabolites with significant correlation with change in yield
yield_sig_cor_rho_spearman_developing_seed_control_HxD = data.frame(cor_res_spearman_rho_developing_seed_control_HxD[1, which(colnames(cor_res_spearman_rho_developing_seed_control_HxD) %in% yield_cor_sig_metabolites_developing_seed_control_HxD)])     #metabolites with significant correlation with change in yield and rho values
colnames(yield_sig_cor_rho_spearman_developing_seed_control_HxD) = "rho value"
yield_sig_cor_rho_spearman_developing_seed_control_HxD = as.matrix(yield_sig_cor_rho_spearman_developing_seed_control_HxD)


#proportion of chalky grains
chalk_50_75_cor_res_spearman_pval_developing_seed_control_HxD = t(data.frame(lapply(yield_chalk_control_developing_seed_HxD[3:ncol(yield_chalk_control_developing_seed_HxD)], function(x) cor.test(yield_chalk_control_developing_seed_HxD$`Change in chalky_50-75%`, x, method = "spearman")$p.value)))    #significance of correlation
chalk_50_75_cor_res_spearman_pval_developing_seed_control_HxD = as.data.frame(chalk_50_75_cor_res_spearman_pval_developing_seed_control_HxD)
colnames(chalk_50_75_cor_res_spearman_pval_developing_seed_control_HxD) = "p.value"
rownames(chalk_50_75_cor_res_spearman_pval_developing_seed_control_HxD) = reduced_metabolite_list_final_developing_seed_2013$Name
chalk_50_75_cor_res_spearman_pval_developing_seed_control_HxD$Significance = with(chalk_50_75_cor_res_spearman_pval_developing_seed_control_HxD,
                                                                                  ifelse(p.value <= 0.001, "***", ifelse(
                                                                                    p.value <= 0.01, "**", ifelse(
                                                                                      p.value < 0.05, "*", "ns"))))
length(which(chalk_50_75_cor_res_spearman_pval_developing_seed_control_HxD$p.value < 0.05))    #3 metabolites
chalk_50_75_cor_sig_metabolites_developing_seed_control_HxD = rownames(chalk_50_75_cor_res_spearman_pval_developing_seed_control_HxD)[chalk_50_75_cor_res_spearman_pval_developing_seed_control_HxD$Significance != "ns"]   #metabolites with significant correlation with change in proportion of grains with >50% chalk
chalk_50_75_sig_cor_rho_spearman_developing_seed_control_HxD = data.frame(cor_res_spearman_rho_developing_seed_control_HxD[2, which(colnames(cor_res_spearman_rho_developing_seed_control_HxD) %in% chalk_50_75_cor_sig_metabolites_developing_seed_control_HxD)])     #metabolites with significant correlation with change in chalky grain proportion and rho values
colnames(chalk_50_75_sig_cor_rho_spearman_developing_seed_control_HxD) = "rho value"
chalk_50_75_sig_cor_rho_spearman_developing_seed_control_HxD = as.matrix(chalk_50_75_sig_cor_rho_spearman_developing_seed_control_HxD)


#correlation heatmap - yield

#matrix with 2 columns of the same data - dummy - to be able to make a 1-column heatmap
yield_sig_cor_rho_spearman_developing_seed_control_HxD_2 = cbind(yield_sig_cor_rho_spearman_developing_seed_control_HxD,
                                                                 yield_sig_cor_rho_spearman_developing_seed_control_HxD)

#create a new matrix with dummy data to make equal-sized heatmap with the others - 14 metabolites in total
dummy_data_2 = matrix(c(0, 0), nrow = 12, ncol = 2, dimnames = list(c(rep("Z", 12))))
yield_sig_cor_rho_spearman_developing_seed_control_HxD_3 = rbind(yield_sig_cor_rho_spearman_developing_seed_control_HxD_2, dummy_data_2)

png("cor_yield_control_developing_seed_HxD.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(yield_sig_cor_rho_spearman_developing_seed_control_HxD_3, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          breaks = c(seq(-1, -0.01, length.out = 300), seq(0.01, 1, length.out = 300)), 
          symbreaks = TRUE, 
          col = bluered(599), 
          colsep = c(0, ncol(yield_sig_cor_rho_spearman_developing_seed_control_HxD_3)), 
          rowsep = c(seq(0, nrow(yield_sig_cor_rho_spearman_developing_seed_control_HxD_3))), 
          sepcolor = "black", 
          sepwidth = c(0.001, 0.005),
          margins = c(20, 35), 
          cexRow = 1, 
          labCol = "", 
          lhei = c(0.45, 5), 
          lwid = c(1.5, 5), 
          density.info = "none",
          key.title = NA, 
          key.xlab = NA,
          key.xtickfun = function() {
            cex <- par("cex")*par("cex.axis")
            side <- 1
            line <- 0
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("-1", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("0", side=side, at=0.5, adj=0.5,
                  line=line, cex=cex, col=col, font=font)
            mtext("1", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          })
legend(x = -0.083, y = 1.06, xpd = TRUE, legend = expression(paste("Correlation coefficient")), bty = "n", cex = 0.8)
legend(0.15, 1.12, xpd = TRUE, legend = c("J"), bty = "n", cex = 1.5)
dev.off()
#crop dummy data at the end of the figure


#correlation heatmap - chalky grain proportion

#matrix with 2 columns of the same data - dummy - to be able to make a 1-column heatmap
chalk_50_75_sig_cor_rho_spearman_developing_seed_control_HxD_2 = cbind(chalk_50_75_sig_cor_rho_spearman_developing_seed_control_HxD,
                                                                       chalk_50_75_sig_cor_rho_spearman_developing_seed_control_HxD)


#create a new matrix with dummy data to make equal-sized heatmap with the others - 11 metabolites in total
dummy_data_3 = matrix(c(0, 0), nrow = 8, ncol = 2, dimnames = list(c(rep("Z", 8))))
chalk_50_75_sig_cor_rho_spearman_developing_seed_control_HxD_3 = rbind(chalk_50_75_sig_cor_rho_spearman_developing_seed_control_HxD_2,
                                                                       dummy_data_3)

png("cor_chalk_50_75_control_developing_seed_HxD.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(chalk_50_75_sig_cor_rho_spearman_developing_seed_control_HxD_3, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          breaks = c(seq(-1, -0.01, length.out = 300), seq(0.01, 1, length.out = 300)), 
          symbreaks = TRUE, 
          col = bluered(599), 
          colsep = c(0, ncol(chalk_50_75_sig_cor_rho_spearman_developing_seed_control_HxD_3)), 
          rowsep = c(seq(0, nrow(chalk_50_75_sig_cor_rho_spearman_developing_seed_control_HxD_3))), 
          sepcolor = "black", 
          sepwidth = c(0.001, 0.005),
          margins = c(25.2, 35), 
          cexRow = 1, 
          labCol = "", 
          lhei = c(0.45, 5), 
          lwid = c(1.5, 5), 
          density.info = "none",
          key.title = NA, 
          key.xlab = NA, 
          key.xtickfun = function() {
            cex <- par("cex")*par("cex.axis")
            side <- 1
            line <- 0
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("-1", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("0", side=side, at=0.5, adj=0.5,
                  line=line, cex=cex, col=col, font=font)
            mtext("1", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          })
legend(x = -0.083, y = 1.06, xpd = TRUE, legend = expression(paste("Correlation coefficient")), bty = "n", cex = 0.8)
legend(0.15, 1.12, xpd = TRUE, legend = c("I"), bty = "n", cex = 1.5)
dev.off()
#crop dummy data at the end of the figure