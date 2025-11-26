#This script was made to compare the different datasets in terms of how the variants are distributed across different regions of titin
#the csv was made by retrieving number of variants in each genomic range in a dataset and compiling them 
#for example: bcftools view -r 2:178618584-178777049 TTN-ABraOM-mSNVs.vcf.gz | grep -v "^#" | wc -l

#directory
setwd("/Users/stella/Desktop/Master_Projekt/vcf_stats")
getwd()

#load data
data <- read.csv("Dataset-Distribution_comparison.csv")
data$mSNV_count <- as.numeric(gsub("[^0-9]", "", data$mSNV_count))
str(data)
summary(data)
print(data)


#make table to make chi-square test to see if the distribution is signficantly different
table_data <- xtabs(mSNV_count ~ dataset + region, data = data)
print(table_data)
chisq_result <- chisq.test(table_data)
print(chisq_result)
expected_table <- chisq_result$expected

#comapre with expected if it were not signficantly different
table_data
expected_table

#check residuals 
residuals <- chisq_result$residuals
print(residuals)

#Graphic Description of Results

#install.packages("graphics")
library("graphics")
mosaicplot(data, shade = TRUE, las=2,
           main = "housetasks")

ordered_regions <- c("Z-Disk", "I-band", "A-band", "M-band")
ordered_datasets<- c("Clinvar", "gnomAD", "ABraOM")
residuals_reord <- residuals[ordered_datasets, ordered_regions]

#install.packages("corrplot")
library(corrplot)
corrplot(residuals_reord,
         tl.col = "black",
         is.cor = FALSE)

data

#compare datasets pairwise

Cvar_vs_gAD <- data[5:12,]
Cvar_vs_gAD


table_data2 <- xtabs(mSNV_count ~ dataset + region, data = Cvar_vs_gAD)
print(table_data2)
chisq_result2 <- chisq.test(table_data2)
print(chisq_result2)
expected_table2 <- chisq_result2$expected

table_data2
expected_table2

residuals2 <- chisq_result2$residuals
print(residuals2)
corrplot(residuals2,
         tl.col = "black",
         is.cor = FALSE)

Cvar_vs_ABr <- data[1:8,]
Cvar_vs_ABr


table_data3 <- xtabs(mSNV_count ~ dataset + region, data = Cvar_vs_ABr)
print(table_data3)
chisq_result3 <- chisq.test(table_data3)
print(chisq_result3)
expected_table3 <- chisq_result3$expected

table_data3
expected_table3

residuals3 <- chisq_result3$residuals
print(residuals3)
corrplot(residuals3,
         tl.col = "black",
         is.cor = FALSE)

gAD_vs_ABr <- data[c(1:4, 9:12), ]
gAD_vs_ABr


table_data4 <- xtabs(mSNV_count ~ dataset + region, data = gAD_vs_ABr)
print(table_data4)
chisq_result4 <- chisq.test(table_data4)
print(chisq_result4)
expected_table4 <- chisq_result4$expected

table_data4
expected_table4

residuals4 <- chisq_result4$residuals
print(residuals4)
corrplot(residuals4,
         tl.col = "black",
         is.cor = FALSE)

