#this script was made to convert csv files to vcf files
#it is made for csv files saved in excel that have this hidden character in the header: '\ufeff'

import csv
from datetime import datetime

#files
csv_file = 'OBSCN-SNVs.csv'  #csv file
vcf_file = 'OBSCN-ABraOM.vcf'  #vcf file

#create vcf header
with open(vcf_file, 'w') as vcf:
    vcf.write('##fileformat=VCFv4.1\n')
    vcf.write(f'##fileDate={datetime.now().strftime("%Y-%m-%d")}\n')
    vcf.write('##source=ABraOM\n')
    vcf.write('##reference=grch38:homo_sapiens\n')
    vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

    #open csv file 
    with open(csv_file, 'r', encoding='utf-8-sig') as infile:
        reader = csv.DictReader(infile)

        #clean filed names
        reader.fieldnames = [field.strip() for field in reader.fieldnames]
        
        # print headers (debugging step)
        print("Cleaned Field Names:", reader.fieldnames)

        for row in reader:
            try:
                #extrct info for INFO field
                chrom = row['Chr']  
                pos = row['Position']  
                ref = row['Ref']  
                alt = row['Alt']  
                rsid = row.get('dbSNP150', '.')  #'.' if missing
                freq = row.get('Frequency', '.')  
                cohort = row.get('Cohort', '.')  

                #add gene and annotation to INFO field
                gene = row.get('Gene', '.')  
                annotation = row.get('Annotation', '.')  

            
                info = f'RSID={rsid};FREQ={freq};COHORT={cohort};GENE={gene};ANNOTATION={annotation}'

                #write vcf
                vcf.write(f'{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t{info}\n')

            except KeyError as e:
                print(f"missing field in row: {e}") #deubigging step
                continue

print(f"VCF file '{vcf_file}' has been created.")

