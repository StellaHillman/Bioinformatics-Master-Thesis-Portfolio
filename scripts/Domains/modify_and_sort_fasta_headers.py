import re
#to sort domains and name the Ig domains seqeuntially 

def modify_fasta_headers(input_fasta, output_fasta):
   
    records = []
    current_header = None
    current_sequence = []

    #read the input FASTA file
    with open(input_fasta, 'r') as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('>'):  #header line
                if current_header:
                    #save the previous record
                    records.append((current_header, ''.join(current_sequence)))
                #extract and clean the header
                match = re.search(r'olaTTN-Ig\d+', line)
                if match:
                    current_header = match.group(0)  #keep only 'olaTTN-IgXX'
                else:
                    current_header = None
                current_sequence = []  #reset sequence
            else:
                #sequence line
                current_sequence.append(line)

        #save the last record
        if current_header:
            records.append((current_header, ''.join(current_sequence)))

    #sort records numerically by 'IgXX' part
    sorted_records = sorted(
        records, key=lambda x: int(re.search(r'\d+', x[0]).group(0))
    )

    #write to the output file
    with open(output_fasta, 'w') as outfile:
        for header, sequence in sorted_records:
            outfile.write(f">{header}\n")
            #split the sequence 
            for i in range(0, len(sequence), 60):
                outfile.write(sequence[i:i+60] + '\n')

#files
input_fasta = 'olaTTN_all_manual.fasta'
output_fasta = 'olaTTN_all_ig_manual_sorted.fasta'
modify_fasta_headers(input_fasta, output_fasta)

