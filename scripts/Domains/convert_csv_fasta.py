import csv

def csv_to_fasta(csv_file, fasta_file):
    """Converts a CSV file with domain sequences to a FASTA file."""
    with open(csv_file, 'r') as infile, open(fasta_file, 'w') as outfile:
        reader = csv.DictReader(infile)
        
        for row in reader:
            #get the domain name and sequence
             header = row['sequence_id']
             sequence = row['sequence']
            
            #write to FASTA format
            fasta_out.write(f'>{header}\n')
        
            #split sequence for better readability
            for i in range(0, len(sequence), 60):
                outfile.write(sequence[i:i+60] + "\n")

    print(f"FASTA file {fasta_file} created successfully!")

if __name__ == "__main__":
    # Replace 'titin_domain.csv' and 'output.fasta' with your actual file paths
    csv_file = 'MSA-new.csv' # Input CSV file
    fasta_file = 'MSA_mafft_all_igs_trimmed_species-type_manually-adjusted_NEW-4-6-25.fasta'  # Output FASTA file

    csv_to_fasta(csv_file, fasta_file)
