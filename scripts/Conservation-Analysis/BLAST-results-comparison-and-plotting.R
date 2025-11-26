#This scrip was made to process BLAST output files and compare sequence identity between domains (titin from human and medaka)
#The input tsv file consists of pairwise BLASTP comparison of all-vs-all domain alignments (every domain in protein A compared to every domain in protein B).
#The script identifies the highest-identity matches (homologous domains), summarizes domain-level similarity, 
#performs statistical tests (normality, ANOVA/Kruskal–Wallis, pairwise comparisons), 
#and generates plots to visualize sequence identity across different titin comparisons.
#2025 Stella Hillman R version 4.4.2

#packages
packages <- c(
  "ggplot2",
  "psych",
  "dplyr",
  "car",
  "stringr",
  "scales",
  "gridExtra"
)

new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

lapply(packages, library, character.only = TRUE)

#directory
setwd("/path/to/files")
getwd()

#blast data medaka B vs human titin FnIII
fnIIBHu <- read.table("medakaBfnIII_vs_humanfnIII_blast_results.tsv", header=FALSE, sep="\t", 
                    col.names = c("query_id", "subject_id", "per_identity", "alignment_length", "mismatches", 
                                  "gap_openings", "query_start", "query_end", "subject_start", "subject_end", 
                                  "evalue", "bitscore"))

#check
head(fnIIBHu)

#clean
fnIIBHu <- na.omit(fnIIBHu)

fnIIBHu$alignment_length <- as.integer(fnIIBHu$alignment_length)
fnIIBHu$evalue <- as.numeric(fnIIBHu$evalue)
fnIIBHu$bitscore <- as.numeric(fnIIBHu$bitscore)

str(fnIIBHu)

fnIIBHu %>% count(query_id)
fnIIBHu %>% count(subject_id)

#summary
summary(fnIIBHu)
describe(fnIIBHu)

#compare highest
highestfnIII <- fnIIBHu%>%
  group_by(query_id) %>%
  top_n(1, per_identity)  

highestfnIII<- highestfnIII %>%
  ungroup()

str(highestfnIII)
describe(highestfnIII)

fnIII_stats <- fnIIBHu%>%
  summarise(mean_identity = mean(per_identity),
            sd_identity = sd(per_identity))
print(fnIII_stats)

hifnIII_stats <- highestfnIII%>%
  summarise(mean_identity = mean(per_identity),
            sd_identity = sd(per_identity))
print(hifnIII_stats)

t.test(df$Größe~df$Geschlecht, var.equal = TRUE, alternative = "two.sided")

setwd("/Users/stella/Desktop/Master_Projekt/ttn-medaka-comparison")

#blast data medaka A vs B A band only
medAB <- read.table("medakaA_vs_B_blast_results.tsv", header=FALSE, sep="\t", 
                      col.names = c("query_id", "subject_id", "per_identity", "alignment_length", "mismatches", 
                                    "gap_openings", "query_start", "query_end", "subject_start", "subject_end", 
                                    "evalue", "bitscore"))

#check
head(medAB)

medAB <- na.omit(medAB)

medAB$alignment_length <- as.integer(medAB$alignment_length)
medAB$evalue <- as.numeric(medAB$evalue)
medAB$bitscore <- as.numeric(medAB$bitscore)

str(medAB)

medAB %>% count(query_id)
medAB %>% count(subject_id)
#summary
summary(medAB)
describe(medAB)

#compare highest
highestAB <- medAB%>%
  group_by(query_id) %>%
  top_n(1, per_identity)  

highestAB<- highestAB %>%
  ungroup()

str(highestAB)
describe(highestAB)

AB_stats <- medAB%>%
  summarise(mean_identity = mean(per_identity),
            sd_identity = sd(per_identity))
print(AB_stats)

hiAB_stats <- highestAB%>%
  summarise(mean_identity = mean(per_identity),
            sd_identity = sd(per_identity))
print(hiAB_stats)

t.test(df$Größe~df$Geschlecht, var.equal = TRUE, alternative = "two.sided")
setwd("/Users/stella/Desktop/Master_Projekt/ttn-medaka-comparison")

#blast data medaka A vs human A band only
medAHu <- read.table("medakaA_vs_human_blast_results.tsv", header=FALSE, sep="\t", 
                     col.names = c("query_id", "subject_id", "per_identity", "alignment_length", "mismatches", 
                                   "gap_openings", "query_start", "query_end", "subject_start", "subject_end", 
                                   "evalue", "bitscore"))

#check
head(medAHu)

medAHu <- na.omit(medAHu )

medAHu$alignment_length <- as.integer(medAHu$alignment_length)
medAHu$evalue <- as.numeric(medAHu$evalue)
medAHu$bitscore <- as.numeric(medAHu$bitscore)

str(medAHu)

medAHu %>% count(query_id)
medAHu %>% count(subject_id)
#summary
summary(medAHu)
describe(medAHu)

#highest 
highestAHu <- medAHu%>%
  group_by(query_id) %>%
  top_n(1, per_identity)  

highestAHu<- highestAHu %>%
  ungroup()

str(highestAHu)
describe(highestAHu)

AHu_stats <- medAHu%>%
  summarise(mean_identity = mean(per_identity),
            sd_identity = sd(per_identity))
print(AHu_stats)

hiAHu_stats <- highestAHu%>%
  summarise(mean_identity = mean(per_identity),
            sd_identity = sd(per_identity))
print(hiAHu_stats)


#blast data medaka B vs human A band only
medBHu <- read.table("medakaB_vs_human_blast_results.tsv", header=FALSE, sep="\t", 
                     col.names = c("query_id", "subject_id", "per_identity", "alignment_length", "mismatches", 
                                   "gap_openings", "query_start", "query_end", "subject_start", "subject_end", 
                                   "evalue", "bitscore"))

#check
head(medBHu)

medBHu <- na.omit(medBHu )

medBHu$alignment_length <- as.integer(medBHu$alignment_length)
medBHu$evalue <- as.numeric(medBHu$evalue)
medBHu$bitscore <- as.numeric(medBHu$bitscore)

str(medBHu)

medBHu %>% count(query_id)
medBHu %>% count(subject_id)

#summary
summary(medBHu)
describe(medBHu)

#highest 
highestBHu <- medBHu%>%
  group_by(query_id) %>%
  top_n(1, per_identity)  

highestBHu<- highestBHu %>%
  ungroup()

str(highestBHu)
describe(highestBHu)

BHu_stats <- medBHu%>%
  summarise(mean_identity = mean(per_identity),
            sd_identity = sd(per_identity))
print(BHu_stats)

hiBHu_stats <- highestBHu%>%
  summarise(mean_identity = mean(per_identity),
            sd_identity = sd(per_identity))
print(hiBHu_stats)

#combine stats

AB_stats$comparison <- c("AB")
AHu_stats$comparison <- c("AHu")
BHu_stats$comparison <- c("BHu")
fnIII_stats$comparison <- c("FnIII")

AB_stats$type <- c("all")
AHu_stats$type <- c("all")
BHu_stats$type <- c("all")
fnIII_stats$type <- c("all")

hiAB_stats$comparison <- c("AB")
hiAHu_stats$comparison <- c("AHu")
hiBHu_stats$comparison <- c("BHu")
hifnIII_stats$comparison <- c("FnIII")

hiAB_stats$type <- c("highest")
hiAHu_stats$type <- c("highest")
hiBHu_stats$type <- c("highest")
hifnIII_stats$type <- c("highest")

combined_stats <- rbind(AB_stats, AHu_stats, BHu_stats,fnIII_stats, 
                        hiAB_stats, hiAHu_stats, hiBHu_stats, hifnIII_stats)

str(combined_stats)
print(combined_stats)

combined_stats <- combined_stats %>% mutate(fill_group = interaction(type, comparison))

#compare all four

medAB$dataset <- "medAB"
medAHu$dataset <- "medAHu"
medBHu$dataset <- "medBHu"
fnIIBHu$dataset <- "fnIIBHu"

combined_data <- rbind(medAB, medAHu, medBHu, fnIIBHu)

str(combined_data)

plot(combined_data)

#normality check
shapiro.test(medAB$per_identity)
shapiro.test(medAHu$per_identity)
shapiro.test(medBHu$per_identity)
shapiro.test(fnIIBHu$per_identity)

qqnorm(medAB$per_identity)
qqnorm(medAHu$per_identity)
qqnorm(medBHu$per_identity)
qqnorm(fnIIBHu$per_identity)

qqline(medAB$per_identity)
qqline(medAHu$per_identity)
qqline(medBHu$per_identity)
qqline(fnIIBHu$per_identity)

kruskal.test(per_identity ~ dataset, data = combined_data)

#levens check for euql variance
library(car)
leveneTest(per_identity ~ dataset, data = combined_data)

#if normal then anova
anova_result <- aov(per_identity ~ dataset, data = combined_data)
summary(anova_result)

#if significant, t-test 
pairwise_ttest <- pairwise.t.test(combined_data$per_identity, combined_data$dataset, p.adjust.method = "bonferroni")
print(pairwise_ttest)

#if variances are unequal welchs test
t_test_result <- t.test(per_identity ~ dataset, data = combined_data, var.equal = FALSE)
print(t_test_result)

#compare highest

highestAB$dataset <- "medAB"
highestAHu$dataset <- "medAHu"
highestBHu$dataset <- "medBHu"
highestfnIII$dataset <- "fnIIBHu"

hicombined_data <- rbind(highestAB, highestAHu, highestBHu, highestfnIII)

str(hicombined_data)

#normality check
shapiro.test(highestAB$per_identity)
shapiro.test(highestAHu$per_identity)
shapiro.test(highestBHu$per_identity)

qqnorm(highestAB$per_identity)
qqline(highestAB$per_identity)
qqnorm(highestAHu$per_identity)
qqline(highestAHu$per_identity)
qqnorm(highestBHu$per_identity)
qqline(highestBHu$per_identity)


#levens check for euql variance*
leveneTest(per_identity ~ dataset, data = hicombined_data)

#anova
hianova_result <- aov(per_identity ~ dataset, data = hicombined_data)
summary(hianova_result)

#t-test 
hipairwise_ttest <- pairwise.t.test(hicombined_data$per_identity, hicombined_data$dataset, p.adjust.method = "bonferroni")
print(hipairwise_ttest)



#### Different Plot Versions ####

ggplot(combined_stats, aes(x = type, y = mean_identity, fill = type)) +
  geom_bar(stat = "identity",  color = "black", position = "dodge", width = 0.5) + 
  ylim(0, 85) +
  geom_errorbar(aes(ymin = mean_identity - sd_identity , ymax = mean_identity + sd_identity), 
                width = 0.2, position = position_dodge(0.5)) +
  labs(title = "Comparison of % identity: Same vs Different Positions in A-band repeats",
       x = "Position Group",
       y = "Mean identity (%)") +
  scale_fill_manual(values = c("Same Position" = "#b053bb", "Different Position" = "#edb4f3")) +  
  theme_minimal()  

#barplot alternative
print(combined_stats)

ggplot(combined_stats, aes(x = comparison, y = mean_identity, fill = fill_group)) +
  geom_col(position = "dodge") +
  labs(title = "Comparison shared identity Medaka titin A and B and Human titin - A-band only",
       x = "Group",
       y = "% Identity") +
  geom_errorbar(aes(ymin = mean_identity - sd_identity , ymax = mean_identity + sd_identity), 
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(
    name = "Type",
    values = c(
      "all.AB" = "#ffc465",
      "all.AHu" = "#a2ab61",
      "all.BHu" = "#578359",
      "all.FnIII" = "#194f42",
      "highest.AB" = "#ffd898",
      "highest.AHu" = "#bcc483",
      "highest.BHu" = "#7bac77",
      "highest.FnIII" = "#368d73"),
    labels = c(
      "all.AB" = "All A-band Igs A-medaka vs B-medaka",
      "all.AHu" = "All A-band Igs A-medaka vs Human",
      "all.BHu" = "All A-band Igs B-medaka vs Human",
      "all.FnIII" = "All FnIII B-medaka vs Human",
      "highest.AB" = "Homologous A-band Igs A-medaka vs B-medaka",
      "highest.AHu" = "Homologous A-band Igs A-medaka vs Human",
      "highest.BHu" = "Homologous A-band Igs B-medaka vs Human",
      "highest.FnIII" = "Homologous FnIII B-medaka vs Human")
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_x_discrete(labels=c("medaka A vs B", "medaka A vs Human", "medaka B vs Human", "FnIII medaka B vs Human")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#with full sequence 
Sq_FnIII_stats <- cbind(Sq_BHu_stats)

#create columns for all types 
AB_stats$comparison <- c("AB")
AHu_stats$comparison <- c("AHu")
BHu_stats$comparison <- c("BHu")
fnIII_stats$comparison <- c("FnIII")

Sq_AB_stats$comparison <- c("AB")
Sq_AHu_stats$comparison <- c("AHu")
Sq_BHu_stats$comparison <- c("BHu")
Sq_FnIII_stats$comparison <- c("FnIII")

AB_stats$type <- c("all")
AHu_stats$type <- c("all")
BHu_stats$type <- c("all")
fnIII_stats$type <- c("all")

hiAB_stats$comparison <- c("AB")
hiAHu_stats$comparison <- c("AHu")
hiBHu_stats$comparison <- c("BHu")
hifnIII_stats$comparison <- c("FnIII")

hiAB_stats$type <- c("highest")
hiAHu_stats$type <- c("highest")
hiBHu_stats$type <- c("highest")
hifnIII_stats$type <- c("highest")


Sq_AB_stats$type <- c("full-Sq")
Sq_AHu_stats$type <- c("full-Sq")
Sq_BHu_stats$type <- c("full-Sq")
Sq_FnIII_stats$type <- c("full-Sq")



combined_stats <- rbind(Sq_AB_stats, Sq_AHu_stats, Sq_BHu_stats, Sq_FnIII_stats, 
                        AB_stats, AHu_stats, BHu_stats, fnIII_stats,
                        hiAB_stats, hiAHu_stats, hiBHu_stats, hifnIII_stats)

str(combined_stats)
print(combined_stats)

#create fill_group for each bar
combined_stats <- combined_stats %>%
  mutate(fill_group = interaction(type, comparison))

#make sure that the bars are ordered like in the dataframe

combined_stats$fill_group <- factor(combined_stats$fill_group, 
                                    levels = unique(combined_stats$fill_group))

#barplot alternative
print(combined_stats)

p <- ggplot(combined_stats, aes(x = comparison, y = mean_identity, fill = fill_group)) +
  geom_col(position = "dodge") +
  labs(title = "Comparison shared identity Medaka titin A and B and Human titin",
       x = "Group",
       y = "% Identity") +
  geom_errorbar(aes(ymin = mean_identity - sd_identity , ymax = mean_identity + sd_identity), 
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(
    name = "Type",
    values = c(
      "full-Sq.AB" = "#e2a544",
      "full-Sq.AHu" = "#8b9441",
      "full-Sq.BHu" = "#41724a",
      "full-Sq.FnIII" = "#14423b",
      "all.AB" = "#ffc465",
      "all.AHu" = "#a2ab61",
      "all.BHu" = "#578359",
      "all.FnIII" = "#194f42",
      "highest.AB" = "#ffd898",
      "highest.AHu" = "#bcc483",
      "highest.BHu" = "#7bac77",
      "highest.FnIII" = "#368d73"),
    labels = c(
      "full-Sq.AB" = "Full protein Seqeunce A-medaka vs B-medaka",
      "full-Sq.AHu" = "Full protein Seqeunce A-medaka vs Human",
      "full-Sq.BHu" = "Full protein Seqeunce B-medaka vs Human",
      "full-Sq.FnIII" = "Full protein Seqeunce B-medaka vs Human",
      "all.AB" = "All A-band Igs A-medaka vs B-medaka",
      "all.AHu" = "All A-band Igs A-medaka vs Human",
      "all.BHu" = "All A-band Igs B-medaka vs Human",
      "all.FnIII" = "All FnIII B-medaka vs Human",
      "highest.AB" = "Homologous A-band Igs A-medaka vs B-medaka",
      "highest.AHu" = "Homologous A-band Igs A-medaka vs Human",
      "highest.BHu" = "Homologous A-band Igs B-medaka vs Human",
      "highest.FnIII" = "Homologous FnIII B-medaka vs Human")
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_x_discrete(labels=c("medaka A vs B", "medaka A vs Human", "medaka B vs Human", "FnIII medaka B vs Human")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, size = 15, hjust = 0.5,  vjust = 1), 
        legend.text = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 22, hjust = 0.4),
        plot.title = element_text(size = 26, hjust = - 0.1 ),
        axis.text.y = element_text(size = 16))
#theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 10))

print(p)

ggsave("medaka vs human titn quadruple plot.png", plot = p, width = 16, height = 13)

#legend text organized differently 

library(ggplot2)
library(dplyr)
library(scales)

#fill_group for each bar, type first for bar order
combined_stats <- combined_stats %>%
  mutate(fill_group = interaction(type, comparison))

#order bars  as in dataframe (type first, comparison second)
combined_stats$fill_group <- factor(combined_stats$fill_group, 
                                    levels = unique(combined_stats$fill_group))

# Plot
p <- ggplot(combined_stats, aes(x = comparison, y = mean_identity, fill = fill_group)) +
  geom_col(position = "dodge") +
  labs(title = "Comparison of Protein Seqeunce Identity: Medaka Titin A, Medaka Titin B and  Human Titin",
       x = "Group",
       y = "% Identity") +
  geom_errorbar(aes(ymin = mean_identity - sd_identity , ymax = mean_identity + sd_identity), 
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(
    name = "Type",
    values = c(
      "full-Sq.AB" = "#e2a544",
      "full-Sq.AHu" = "#8b9441",
      "full-Sq.BHu" = "#41724a",
      "full-Sq.FnIII" = "#14423b",
      "all.AB" = "#ffc465",
      "all.AHu" = "#a2ab61",
      "all.BHu" = "#578359",
      "all.FnIII" = "#194f42",
      "highest.AB" = "#ffd898",
      "highest.AHu" = "#bcc483",
      "highest.BHu" = "#7bac77",
      "highest.FnIII" = "#368d73"
    ),
    labels = c(
      "full-Sq.AB" = "Full protein Seqeunce A-medaka vs B-medaka",
      "full-Sq.AHu" = "Full protein Seqeunce A-medaka vs Human",
      "full-Sq.BHu" = "Full protein Seqeunce B-medaka vs Human",
      "full-Sq.FnIII" = "Full protein Seqeunce B-medaka vs Human",
      "all.AB" = "All A-band Igs A-medaka vs B-medaka",
      "all.AHu" = "All A-band Igs A-medaka vs Human",
      "all.BHu" = "All A-band Igs B-medaka vs Human",
      "all.FnIII" = "All FnIII B-medaka vs Human",
      "highest.AB" = "Homologous A-band Igs A-medaka vs B-medaka",
      "highest.AHu" = "Homologous A-band Igs A-medaka vs Human",
      "highest.BHu" = "Homologous A-band Igs B-medaka vs Human",
      "highest.FnIII" = "Homologous FnIII B-medaka vs Human"
    ),
    # This controls legend order: comparison first, then type
    breaks = c(
      "full-Sq.AB",
      "all.AB",
      "highest.AB",
      "full-Sq.AHu",
      "all.AHu",
      "highest.AHu",
      "full-Sq.BHu",
      "all.BHu",
      "highest.BHu",
      "full-Sq.FnIII",
      "all.FnIII",
      "highest.FnIII"
    )
  ) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_x_discrete(labels = c("medaka A vs B", "medaka A vs Human", "medaka B vs Human", "FnIII: medaka B vs Human")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, size = 15, hjust = 0.5, vjust = 1), 
    legend.text = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, hjust = 0.4),
    plot.title = element_text(size = 22),
    axis.text.y = element_text(size = 16)
  )

print(p)

ggsave("medaka vs human titn quadruple plot.png", plot = p, width = 16, height = 13)

combined_stats <- rbind(Sq_AB_stats, Sq_AHu_stats, Sq_BHu_stats, Sq_FnIII_stats,
                        hiAB_stats, hiAHu_stats, hiBHu_stats, hifnIII_stats)

combined_stats <- combined_stats %>%
  mutate(fill_group = interaction(type, comparison))
combined_stats$fill_group <- factor(combined_stats$fill_group, 
                                    levels = unique(combined_stats$fill_group))



combined_stats


library(scales)
library(stringr)
library(dplyr)

# Plot
p <- ggplot(combined_stats, aes(x = comparison, y = mean_identity, fill = fill_group)) +
  geom_col(position = "dodge") +
  labs(title = "Comparison of Protein Seqeunce Identity: Medaka Titin A, Medaka Titin B and  Human Titin",
       x = "Group",
       y = "% Identity") +
  geom_errorbar(aes(ymin = mean_identity - sd_identity , ymax = mean_identity + sd_identity), 
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(
    name = "Type",
    values = c(
      "all.AB" = "#e2a544",
      "all.AHu" = "#8b9441",
      "all.BHu" = "#41724a",
      "all.FnIII" = "#14423b",
      "highest.AB" = "#ffd898",
      "highest.AHu" = "#bcc483",
      "highest.BHu" = "#7bac77",
      "highest.FnIII" = "#368d73"
    ),
    labels = c(
      "all.AB" = "A-band Igs A-medaka vs B-medaka",
      "all.AHu" = "A-band Igs A-medaka vs Human",
      "all.BHu" = "A-band Igs Seqeunce B-medaka vs Human",
      "all.FnIII" = "FnIII A-band  B-medaka vs Human",
      "highest.AB" = "Homologous A-band Igs A-medaka vs B-medaka",
      "highest.AHu" = "Homologous A-band Igs A-medaka vs Human",
      "highest.BHu" = "Homologous A-band Igs B-medaka vs Human",
      "highest.FnIII" = "Homologous FnIII B-medaka vs Human"
      ),
    # This controls legend order: comparison first, then type
    breaks = c(
      "all.AB",
      "highest.AB",
      "all.AHu",
      "highest.AHu",
      "all.BHu",
      "highest.BHu",
      "all.FnIII",
      "highest.FnIII"
    )
  ) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_x_discrete(labels = str_wrap(c("Ig: medaka A vs B", "Ig: medaka A vs Human", 
                                       "Ig: medaka B vs Human", "FnIII: medaka B vs Human"), 
                                     width = 8)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, size = 12, hjust = 0.5, vjust = 1), 
    legend.text = element_text(size = 10, margin = margin(l = 5, b = 5, t = 5, unit = "pt")),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, hjust = 0.4),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 8),
      legend.spacing.x = unit(1.0, 'cm')
  ) +
  guides(fill = guide_legend(label.theme = element_text(lineheight = 0.9)))


p

ggsave("medaka vs human titn quadruple plot.png", plot = p, width = 16, height = 13)

#### end ####

#Barplot preperations

#use gridExtra to create a grid
grid.arrange(p, p2, ncol = 2)   

#make a shared legend 
#source for this part: https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p)


#to align x axis text
common_margin <- margin(t = 10, r = 10, b = 20, l = 10)
p  <- p  + theme(plot.margin = common_margin, legend.position = "none")
p2 <- p2 + theme(plot.margin = common_margin, legend.position = "none")

#combine plots
grid.arrange(arrangeGrob(p, p2, nrow=1), 
             widths=c(4,1),
             mylegend, 
             ncol=2, 
             heights=c(10, 1))


summary(highestAB)
summary(highestAHu)
summary(highestBHu)



library(stringr)
library(scales)
library(ggplot2)



### This was the finla plot used in my thesis, 

#changing legend breaks/text wrap
p <- ggplot(combined_stats, aes(x = comparison, y = mean_identity, fill = fill_group)) +
  geom_col(position = "dodge") +
  labs(title = "",
       x = "Group",
       y = "% Identity") +
  geom_errorbar(aes(ymin = mean_identity - sd_identity , ymax = mean_identity + sd_identity), 
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(
    name = "Comparison Group",
    values = c(
      "all.AB" = "#e2a544",
      "all.AHu" = "#8b9441",
      "all.BHu" = "#41724a",
      "all.FnIII" = "#14423b",
      "highest.AB" = "#ffd898",
      "highest.AHu" = "#bcc483",
      "highest.BHu" = "#7bac77",
      "highest.FnIII" = "#368d73"
    ),
    labels = c(
      "all.AB" = "Full protein Seqeunce\nA-medaka vs B-medaka",
      "all.AHu" = "Full protein Seqeunce\nA-medaka vs Human",
      "all.BHu" = "Full protein Seqeunce\nB-medaka vs Human",
      "all.FnIII" = "Full protein Seqeunce\nB-medaka vs Human",
      "highest.AB" = "Homologous A-band Igs\nA-medaka vs B-medaka",
      "highest.AHu" = "Homologous A-band Igs\nA-medaka vs Human",
      "highest.BHu" = "Homologous A-band Igs\nB-medaka vs Human",
      "highest.FnIII" = "Homologous FnIII\nB-medaka vs Human"
    ),
    breaks = c(
      "all.AB",
      "highest.AB",
      "all.AHu",
      "highest.AHu",
      "all.BHu",
      "highest.BHu",
      "all.FnIII",
      "highest.FnIII"
    )
  ) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_x_discrete(labels = str_wrap(c("Ig: medaka A vs B", "Ig: medaka A vs Human", 
                                       "Ig: medaka B vs Human", "FnIII: medaka B vs Human"), 
                                     width = 8)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, hjust = 0.5),
    legend.text = element_text(size = 10, lineheight = 0.9, 
                               margin = margin(l = 5, b = 5, t = 5, unit = "pt")),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.5, "cm"), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, hjust = 0.4),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 8),
    legend.spacing.x = unit(1, 'cm')
  ) +
  guides(fill = guide_legend(label.theme = element_text(lineheight = 1.0)))

p


