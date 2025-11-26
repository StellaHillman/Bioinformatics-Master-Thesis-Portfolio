import re

#this is for the sorting of domains

def parse_fasta(filename):
    with open(filename, 'r') as file:
        fasta_entries = []
        entry = {"header": "", "sequence": ""}
        
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if entry["header"]:  #save the previous entry
                    fasta_entries.append(entry)
                entry = {"header": line, "sequence": ""}
            else:
                entry["sequence"] += line
        #add the last entry after loop
        if entry["header"]:
            fasta_entries.append(entry)
    
    return fasta_entries

def get_number_from_header(header):
    #use regex to extract the number after e.g. 'FN3_'
    match = re.search(r"FN3_(\d+)", header)
    if match:
        return int(match.group(1))  #return the number as an integer for sorting
    else:
        return float('inf')  #in case the number is missing, place it at the end

def sort_fasta_by_number(fasta_entries):
    #sort entries by the extracted number from the header
    sorted_entries = sorted(fasta_entries, key=lambda entry: get_number_from_header(entry["header"]))
    return sorted_entries

def write_fasta(filename, sorted_entries):
    with open(filename, 'w') as file:
        for entry in sorted_entries:
            file.write(f"{entry['header']}\n")
            #wrap sequence at 80 characters per line
            sequence = entry['sequence']
            for i in range(0, len(sequence), 80):
                file.write(f"{sequence[i:i+80]}\n")

#file paths
input_file = "twc_fn3_manually-adjusted_domains.fasta"
output_file = "Twc_fn3-domains_sorted.fasta"

#process the file
fasta_entries = parse_fasta(input_file)
sorted_entries = sort_fasta_by_number(fasta_entries)
write_fasta(output_file, sorted_entries)

print(f"Sorted sequences have been saved to {output_file}.")

