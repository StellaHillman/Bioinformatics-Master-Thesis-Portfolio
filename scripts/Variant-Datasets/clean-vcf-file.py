#this script is made to remove whitespace from the vcf files as this was a problem
#modiefied for each vcf file
#this script was made with help of chatgtp

import re

#files
vcf_input_path = "OBSCN-ABraOM_mSNVs-fix.vcf"
vcf_output_path = "OBSCN-ABraOM_mSNVs-nowhite.vcf"

#function to fix the vcf file
def fix_vcf(vcf_input_path, vcf_output_path):
    with open(vcf_input_path, 'r') as infile, open(vcf_output_path, 'w') as outfile:
        lines = infile.readlines()
        new_header = []
        header_found = False
        info_tags_added = False

        for line in lines:
            if line.startswith("##"):
                new_header.append(line)

            elif line.startswith("#CHROM"):
                new_header.append(line)  #add the column names line (e.g., #CHROM)
                header_found = True

            else:
                #for the data lines (rows)
                fields = line.strip().split("\t")

                #make sure there are exactly 8 columns. if not add emtpy values
                if len(fields) < 8:
                    print(f"Warning: Line has fewer than 8 columns: {line.strip()}")
                    while len(fields) < 8:
                        fields.append(".")
                    line = "\t".join(fields) + "\n"

                #clean up the ANNOTATION field to remove extra spaces
                info_field = fields[7]
                if "ANNOTATION=" in info_field:
                    # Remove extra space between nonsynonymous and SNV
                    info_field = re.sub(r"(ANNOTATION=nonsynonymous)\s+SNV", r"\1SNV", info_field)
                    fields[7] = info_field

                #join the fields back into a line
                line = "\t".join(fields) + "\n"
                new_header.append(line)

        #write to output file
        outfile.writelines(new_header)

#run function to clean file
fix_vcf(vcf_input_path, vcf_output_path)
print(f"file saved {vcf_output_path}.")

