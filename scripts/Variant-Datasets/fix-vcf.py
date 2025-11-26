#this was made to fix the files that was converted to vcf so that it would work for JalVar
#this script was made with help of chatgtp

import re

#files
vcf_input_path = "OBSCN-ABraOM_mSNVs.vcf"
vcf_output_path = "OBSCN-ABraOM_mSNVs-fix.vcf"

#define the header lines that are missing
new_header_lines = [
    "##fileformat=VCFv4.2",  
    "##fileDate=2025-04-24",  
    "##source=ABraOM",  
    "##reference=grch38:homo_sapiens",  
    "##contig=<ID=1,length=248956422>",  
    "##contig=<ID=2,length=243199373>",  
    "##INFO=<ID=FREQ,Number=A,Type=Float,Description=\"Frequency of variant\">",
    "##INFO=<ID=COHORT,Number=1,Type=String,Description=\"Cohort description\">",
    "##INFO=<ID=RSID,Number=1,Type=String,Description=\"RSID from dbSNP\">",
    "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene associated with the variant\">",
    "##INFO=<ID=ANNOTATION,Number=1,Type=String,Description=\"Annotation of the variant\">"
]

#function to fix the VCF file
def fix_vcf(vcf_input_path, vcf_output_path):
    with open(vcf_input_path, 'r') as infile, open(vcf_output_path, 'w') as outfile:
        lines = infile.readlines()
        new_header = []
        header_found = False
        fileformat_found = False
        info_tags_added = False

        for line in lines:
            if line.startswith("##"):
                #skip duplicate fileformat line
                if line.startswith("##fileformat=") and fileformat_found:
                    continue
                if line.startswith("##fileformat="):
                    fileformat_found = True
                
                #skip duplicate contig lines
                if "##contig=" in line:
                    if any(contig in line for contig in ["ID=1", "ID=2"]):
                        new_header.append(line)
                        continue
                    #skip if contig line for CHROM 1 or CHROM 2 is already in the header
                    if "ID=1" not in line and "ID=2" not in line:
                        new_header.extend(new_header_lines)
                else:
                    new_header.append(line)

            elif line.startswith("#CHROM"):
                #add the new INFO fields if missing before the #CHROM line
                if not info_tags_added:
                    new_header.extend(new_header_lines)
                    info_tags_added = True
                
                #add the #CHROM line
                new_header.append(line)
                header_found = True

            else:
                #foer the data lines (rows)
                fields = line.strip().split("\t")

                #mking sure there are exactly 8 columns, add empty values if not
                if len(fields) < 8:
                    print(f"Warning: Line has fewer than 8 columns: {line.strip()}")
                    while len(fields) < 8:
                        fields.append(".")
                    line = "\t".join(fields) + "\n"

                #foe cases where REF is '-' (deletion), replace with 'N'
                ref_allele = fields[3]
                if ref_allele == "-":
                    print(f"Warning: Found REF allele '-' at {fields[0]}:{fields[1]}.")
                    fields[3] = "N"
                    line = "\t".join(fields) + "\n"

                #append the corrected data line to the header
                new_header.append(line)

        #add a blank line before the #CHROM header if its missing
        if len(new_header) > 0 and not new_header[-1].endswith("\n"):
            new_header.append("\n")

        #write the modified content to the output file
        outfile.writelines(new_header)

#run the function to fix the VCF
fix_vcf(vcf_input_path, vcf_output_path)
print(f"File has been fixed and saved as {vcf_output_path}.")

