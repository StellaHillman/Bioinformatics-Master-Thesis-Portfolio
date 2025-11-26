Variant Databases
=========================

Scripts used to convert csv to vcf: "convert_csv_to_vcf.py"
As these formats weren't good yet they had to be fixed using: "fix-vcf.py" and "clean-vcf-file.py", indivdual fixes were made for each file 


ClinVar:
- Source: https://github.com/MorellThomas/JalVar 
- Original source: https://www.ncbi.nlm.nih.gov/clinvar/ 
- TTN VCF creation date: 2024-04-16
- OBSCN VCF creation date: 2024-08-05

gnomAD (v4.1.0):
- TTN (ENSG00000155657.29): exported 2024-12-04 as CSV, converted to VCF
- OBSCN (ENSG00000154358.23): exported 2025-03-31 as CSV, converted to VCF

ABraOM:
- TTN: downloaded 2024-12-03 as CSV, converted to VCF
- OBSCN: downloaded 2025-04-24 as CSV, converted to VCF

Non-human species:
- Dummy VCF files created manually 

Variant Distribution Comparison
------------------------------
Venn diagram 
- counted overlaps using: "vcf_intersections.py"
- plotted using: "Venn_Diagram.R" (R v4.4.2) 


SNV Distribution Analysis:
- counted SNVs within defined genomic limits via command line; saved as "Dataset-Distribution_comparison.csv"
- statistical analysis Chi-squre-datasets.R (R v4.4.2) 
- plotted using a modified version of `plot_vcf.R`. Script source: https://doi.org/10.5281/zenodo.11453080 (accessed: 2025-04-23, Kamaraj & Sinha, 2024)
