# Titin Ig/Fn-III Domain Numbering Script
# =======================================
#This script was made to convert Ig and FnIII annotations from csv files into a modified system specific to titin for structural and comparative analyses 
#original annotations IgStrand System (Tawfeeq et al., 2024) https://github.com/ncbi/icn3d.git on the 2025-05-02 using both experimental pdb and alphafold models
#Numbering Explentation:
# - Each residue is annotated as [a,b,c]:
#     a = β-strand number (single-digit; titin Ig domains have 9 strands, Fn-III domains 7)
#     b = region type: 0 = β-strand, 1 = loop/coiled region
#     c = position relative to anchor (0 = anchor, negative = N-terminal, positive = C-terminal)
# - Loop numbering follows the same logic as IgStrand System: first half relative to preceding strand,
#   second half relative to the following strand.
#In this script new annotaions are calculated with function "calculate_new_annotation":
# x -> strand number (a)
# y -> region type (b)
# z -> position relative to anchor (c)
#2025 Stella Hillman R version 4.4.2
#===================================

#packages
library(dplyr)

#set directory
setwd("/path/to/directory")
getwd()

#load data
Igs <- read.csv("numbered-Ig-ttn.csv")
head(Igs)
str(Igs)

#final version: function to compute x, y, z

calculate_new_annotation <- function(annotations) {
  #remove whitespace
  annotations <- trimws(annotations)
  
  #extract first uppercase letter (β-strand A,A',B etc) 
  prefix <- sub("^([A-Z]'?).*$", "\\1", annotations)
  #extract _loop sufixes for residues in loop regions
  sufix  <- sub("^.*?(_loop)", "\\1", annotations)
  #extract digists 3-5 for residue position
  number <- as.integer(sub("^.*?(\\d{3,5}).*$", "\\1", annotations))
  
  #make sure extraction worked
  print("cleaned prefix and numbers")
  print(data.frame(annotations, prefix, number))
  
  # x = strand number according to modified IgStrand
  x <- case_when(
    prefix == "A"  ~ "1",
    prefix == "A'" ~ "2",
    prefix == "B"  ~ "3",
    prefix == "C"  ~ "4",
    prefix == "C'" ~ "5",
    prefix == "D"  ~ "6",
    prefix == "E"  ~ "7",
    prefix == "F"  ~ "8",
    prefix == "G"  ~ "9",
    TRUE           ~ "NA" 
  )
  
  # anchor-based numbering for z
  bucket <- floor(number / 100)
  anchor <- bucket * 100 + 50
  z <- number - anchor
  
  # y = loop/coil or strand differentitiation
  y <- case_when(
    sufix == "_loop" ~ "1",
    TRUE             ~ "0"
  )
  
  #check if it worked
  print("computed x, y, z")
  print(data.frame(annotations, x, y, z))
  
  #return data.frame with separate columns
  return(data.frame(x = x, y = y, z = z))
}


#apply final function to all *_annot columns

annot_cols <- grep("_annot$", names(Igs), value = TRUE)
for (col in annot_cols) {
  new_values <- calculate_new_annotation(Igs[[col]])
  
  #separate x, y, z columns
  Igs[[paste0(col, "_x")]] <- new_values$x
  Igs[[paste0(col, "_y")]] <- new_values$y
  Igs[[paste0(col, "_z")]] <- new_values$z
  
  #combined annotation column
  Igs[[paste0(col, "_new")]] <- paste0("(", new_values$x, ", ", new_values$y, ", ", new_values$z, ")")
}

#save annotated data
write.csv(Igs,
          file='/Users/stella/Desktop/Master_Projekt/IgStrand/FnIII-numbering_new-annotations.csv',
          na='')

##### previous version and simpler apply loops #####
#apply to all *_annot columns (basic version)
annot_cols <- grep("_annot$", names(Igs), value = TRUE)
for (col in annot_cols) {
  new_col <- paste0(col, "_new")
  Igs[[new_col]] <- calculate_new_annotation(Igs[[col]])
}


##### things used to check and debugg ####
str(Igs)
head(select(Igs, contains("")))
sapply(colnames(Igs), print)
colnames_df <- colnames(Igs)
str(colnames_df)
domains <- sub("_.*", "", colnames_df)
str(domains)
print_column_names <- function(Igs) { for (col in colnames(Igs)) print(col) }
print_column_names(Igs)



##### create newIg dataframe grouped by first row (not used in final) ####
newIg <- data.frame()
firstrow <- Igs[,1]
names <- unique(firstrow)
for (n in names) {
newIg <- bind_cols(newIg, frame[3,which(firstrow == n)])
 }
colnames() <- names;

##### Optional: reorder columns by domain based names ####
all_cols <- colnames(Igs)
base_names <- sub("(_annot(_[a-z]+)?|_annot_new)?$", "", all_cols)
unique_bases_in_order <- unique(base_names)
reordered_cols <- unlist(lapply(unique_bases_in_order, function(b) {
  all_cols[base_names == b]
}))
Igs_reordered <- Igs[, reordered_cols]

#save reordered dataframe
write.csv(Igs_reordered,
          file='/Users/stella/Desktop/Master_Projekt/IgStrand/Ig-numbering_new-annotations.csv',
          na='')



