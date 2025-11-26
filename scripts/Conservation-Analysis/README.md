Comparative Analysis of Human and Medaka Titin Domain
=====================================================

sequence identities between a-band domains of human and medaka titin A and B were calculated using blast+ v2.16.0 (camacho et al., 2009). 

Databases made with blast makeblastdb -in domain-sequences.fasta -dbtype prot -out output_db: BLAST-dbs
% Identity calculated with: blastp -query domain-sequences.fasta -db output_db -out results.tsv -outfmt 6: BLAST-Results

Statistical analysis and Barplots made with: BLAST-results-comparison-and-plotting.R
Boxplots with Boxplots.R 


