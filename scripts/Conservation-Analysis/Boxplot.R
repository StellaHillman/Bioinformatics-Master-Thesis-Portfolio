#this script was made to create boxplot of homoglus titin a band ig domains from human, medaka A and medaka B
#only works after running "BLAST-results-comparison-and-plotting.R"
#this script was also used to check the datasets and make corrections
#for a few cases a better match was outside the domain, so the equivalent domain within the repeat pattern was used instead

p2 <- ggplot(combinedHi4, aes(x = comparison, y = per_identity, fill = comparison)) + #combinedHi4 were the final datapoints, see further down for how this was built
  geom_boxplot(outlier.shape = NA, alpha = 1) +
  geom_jitter(width = 0.2,  alpha = 0.3) +
  scale_fill_manual(
    name = "Type",
    values = c(
      "AB" = "#ffd898",
      "AHu" = "#bcc483",
      "BHu" = "#7bac77"
    ),
    labels = c(
      "AB" = "Homologous A-band Igs A-medaka vs B-medaka",
      "AHu" = "Homologous A-band Igs A-medaka vs Human",
      "BHu" = "Homologous A-band Igs B-medaka vs Human"
    )
  ) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_x_discrete(labels = str_wrap(c("Ig: medaka A vs B", "Ig: medaka A vs Human", 
                                       "Ig: medaka B vs Human", "FnIII: medaka B vs Human"), 
                                     width = 8)) +
  theme_minimal() +
  labs(title = "", y = "% Identity") +
  theme(
    axis.text.x = element_text(angle = 0, size = 12, hjust = 0.5, vjust = 1), 
    legend.text = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, hjust = 0.4),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 8)
  )
p2

#this plot was combined with the other one in BLAST-results-comparison-and-plotting.R (line 631)


#data merging and checking
hiMerge <- merge(highestAB, highestAHu, by= "subject_id") #merge best hits between AB and AHu by subject_id
summary(hiMerge)

Himatch <- data.frame(hiMerge$subject_id,hiMerge$query_id.x, hiMerge$query_id.y) #just checking what matches were there

#check for duplicates
query_counts <- table(combined_data$query_id)
any(query_counts >= 3) #check if any query shows up 3 times or more

query_triplicates <- names(query_counts[query_counts >= 3])
print(query_triplicates)

subject_counts <- table(combined_data$subject_id)
any(subject_counts >= 3) #same for subject ids

subject_triplicates <- names(subject_counts[subject_counts >= 3])
print(subject_triplicates)

#triplicates checked manually
triplicate_subjects <- c("olaTTNAIg109/193", "olaTTNAIg64/193", "olaTTNAIg89/196", "olaTTNAIg95/199")
triplicate_data <- combinedHi[combined_data$subject_id %in% triplicate_subjects, ]
print(triplicate_data)

#summary stats
nr_comparisons <- combinedHi %>%
  group_by(comparison) %>%
  summarise(count = n())
print(nr_comparisons)

#how many unique queries per comparison
unique_query <- combinedHi %>%
  group_by(comparison) %>%
  summarise(unique_query_ids = n_distinct(query_id))
print(unique_query)

#same for subject
unique_subject <- combinedHi3 %>%
  group_by(comparison) %>%
  summarise(unique_subject_ids = n_distinct(subject_id))
print(unique_subject)

#filter some data for manual checking
combined2 <-combined_data %>% filter(subject_id %in% 
                                       c("olaTTNAIg109/193", "olaTTNAIg64/193", "olaTTNAIg89/196", "olaTTNAIg95/199"))
combined2

#assign comparison labels for merging later
highestAB$comparison <- "AB"
highestAHu$comparison <- "AHu"
highestBHu$comparison <- "BHu"

#manually add hits that were better aligned than the original blast match
combined3 <- combined_data %>% filter(query_id == "olaTTNnewIg_1104/1104")
olaTTNnewIg_1104 <- combined3 %>% filter(subject_id == "TTNA15")
olaTTNnewIg_1104$comparison <- "BHu"

combined3 <- combined_data %>% filter(query_id == "olaTTNIg86")
olaTTNIg86 <- combined3 %>% filter(subject_id == "olaTTNAIg86/197")
olaTTNIg86$comparison <- "AB"

combined3 <- combined_data %>% filter(query_id == "TTNA83")
TTNA83 <- combined3 %>% filter(subject_id == "olaTTNAIg86/197")
TTNA83$comparison <- "AHu"

#TTNA138 was removed at some point
str(TTNA138)
str(olaTTNnewIg_1104)
str(olaTTNIg86)
str(TTNA83)

#final version of combinedHi3 before adding manually corrected hits
str(combinedHi3)
str(combined_data)
str(medAB)
str(medAHu)
str(medBHu)
str(fnIIBHu)

#add manually corrected alignments
newdomains <- rbind(olaTTNnewIg_1104, olaTTNIg86, TTNA83)
combinedHi4 <- rbind(combinedHi3, newdomains)
summary(combinedHi4)

#check for duplicates again
query_counts <- table(combinedHi4$query_id)
any(query_counts >= 3)
query_triplicates <- names(query_counts[query_counts >= 3])
print(query_triplicates)

subject_counts <- table(combinedHi4$subject_id)
any(subject_counts >= 3)
subject_triplicates <- names(subject_counts[subject_counts >= 3])
print(subject_triplicates)

#this one was identified as a duplicate
triplicate_subjects <- c("olaTTNAIg101/196")
triplicate_data <- combinedHi4[combinedHi4$subject_id %in% triplicate_subjects, ]
triplicate_data

#stats
nr_comparisons <- combinedHi %>%
  group_by(comparison) %>%
  summarise(count = n())
print(nr_comparisons)

unique_query <- combinedHi4 %>%
  group_by(comparison) %>%
  summarise(unique_query_ids = n_distinct(query_id))
print(unique_query)

unique_subject <- combinedHi4 %>%
  group_by(comparison) %>%
  summarise(unique_subject_ids = n_distinct(subject_id))
print(unique_subject)

#check distribution
hist(combinedHi4$per_identity)
qqnorm(combinedHi4$per_identity)
qqline(combinedHi4$per_identity)

#test for normality
shapiro.test(combinedHi4$per_identity[combinedHi4$comparison == "AB"])
shapiro.test(combinedHi4$per_identity[combinedHi4$comparison == "AHu"])
shapiro.test(combinedHi4$per_identity[combinedHi4$comparison == "BHu"])

#test for equal variance
library(car)
leveneTest(per_identity ~ comparison, data = combinedHi4)

#anova if normal
anova_result <- aov(per_identity ~ comparison, data = combinedHi4)
summary(anova_result)

#post-hoc
TukeyHSD(anova_result)

#pairwise t-tests
pairwise_ttest <- pairwise.t.test(combinedHi4$per_identity, combinedHi4$comparison, p.adjust.method = "bonferroni")
print(pairwise_ttest)

#welsch test in case of unequal variance
t_test_result <- t.test(per_identity ~ comparison, data = combinedHi4, var.equal = FALSE)
print(t_test_result)


#identify outliers
Q1 <- quantile(combinedHi3$per_identity, 0.25)
Q3 <- quantile(combinedHi3$per_identity, 0.75)
outliers <- combinedHi4[combinedHi4$per_identity < Q1 | combinedHi4$per_identity > Q3, ]
summary(outliers)

#10% and 90% percentile cutoff
twenty <- quantile(combinedHi4$per_identity, 0.1)
eighty <- quantile(combinedHi4$per_identity, 0.9)
outliers <- combinedHi4[combinedHi4$per_identity < twenty | combinedHi4$per_identity > eighty, ]
outliers1 <- combinedHi4[combinedHi4$per_identity < twenty, ]
outliers2 <- combinedHi4[combinedHi4$per_identity > eighty, ]
print(outliers1, n=30)
print(outliers2, n=30)

#quick check of structure and format
colnames(combined_data)
class(combined_data)
class(combined_data$dataset)
print(head(combined_data))

combined_data <- as_tibble(combined_data)

#how many unique query-subject pairs per dataset
unique_counts <- combined_data %>%
  group_by(dataset = .data$dataset) %>%
  summarize(unique_pairs = n_distinct(paste0(.data$query_id, "_", .data$subject_id)))
print(unique_counts)

#fix comparison column (if needed)
unique_counts$comparison <- recode(unique_counts$comparison,
                                   "medAB" = "AB",
                                   "medAHu" = "AHu",
                                   "medBHu" = "BHu",
                                   "fnIIBHu" = "FnIII")

#total number of unique queries
unique_query_count <- combined_data %>%
  summarize(unique_queries = n_distinct(query_id))
print(unique_query_count)

#unique queries by dataset
unique_query_counts_by_dataset <- combined_data %>%
  group_by(dataset) %>%
  summarize(unique_queries = n_distinct(query_id))
print(unique_query_counts_by_dataset)

#msa stuff for ggmsa plot
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("ggmsa")
library(ggmsa)
library(Biostrings)

setwd("/Users/stella/Desktop/Master_Projekt/MSA_Final")
IgMSA <- readAAMultipleAlignment("MSA_mafft_all_igs_trimmed_species-type_manually-adjusted_NEW-4-6-25.fasta")
IgMSA1 <- readAAMultipleAlignment("MSA-for-figure_first-third.fasta")
ggmsa(IgMSA)

#white color msa plot
my_cutstom <- data.frame(names = c(LETTERS[1:26],"-"), color = "#FFFFFF", stringsAsFactors = FALSE)

ggmsa(IgMSA1, 
      custom_color = my_cutstom, 
      char_width = 0.5, 
      show.legend = FALSE,
      none_bg = TRUE,
      seq_name = TRUE)

#rename sequences from external file
names1 <- readLines("names1.txt")
seqs <- as.character(unmasked(IgMSA1))
IgMSA1 <- AAStringSet(IgMSA1)
names(IgMSA1) <- names1

#plot with new names
ggmsa(IgMSA1,
      custom_color = my_cutstom, 
      char_width = 0.5, 
      show.legend = FALSE,
      none_bg = TRUE,
      seq_name = TRUE)

#cleaner version
msa <- readAAMultipleAlignment("MSA-for-figure_first-third.fasta")
seqs <- as(msa, "AAStringSet")
if(length(names1) != length(seqs)) stop("Mismatch between number of sequences and number of names.")
names(seqs) <- names1
ggmsa(seqs, 1, 127, custom_color = my_cutstom, char_width = 0.5, show.legend = FALSE, none_bg = TRUE, seq_name = TRUE)

#check unique queries from cleaned data
cleaned_data %>%
  filter(dataset == "medAHu") %>%
  summarize(unique_queries = n_distinct(query_id))

cleaned_data %>%
  filter(dataset == "medAHu") %>%
  filter(query_id == "TTNA4" | subject_id == "TTNA4")





