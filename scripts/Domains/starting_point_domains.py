from Bio import SeqIO

#this script is to make extended domain seqeunces to then manually define the boundries
#for this you need align the initial domains to the full protein sequence
#this can be done using mafft --addfragments 

def extend_domain_boundaries(alignment_file, output_file, n_extend=10, c_extend=14): #define how many amino acids the sequences should be extended by
    #read the seqeunces
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(alignment_file, "fasta")}

    #identify the full-length protein sequence by looking for the base name e.g. 'Twc-seq'
    full_length_id = next(
        (key for key in sequences if key.startswith("Twc-seq")),  
        None
    )
    if not full_length_id:
        raise ValueError("Full-length titin sequence not found in the alignment file.")
    
    full_length_seq = sequences[full_length_id].replace("-", "")  #remove alignment gaps    
    #store extended domains
    extended_domains = []
    
    for record_id, seq in sequences.items():
        if record_id == full_length_id:
            continue  #skip the full-length sequence
        
        #find the start and end positions of the domain in the aligned sequence
        domain_start = seq.find("-")  #adjust for leading gaps
        domain_end = seq.rfind("-")
        if domain_start == -1 or domain_end == -1:
            continue  #skip records without domain content
        
        #trim gaps to find the actual positions in the full-length sequence
        domain_seq = seq.replace("-", "")
        full_start = full_length_seq.find(domain_seq)
        full_end = full_start + len(domain_seq)
        
        if full_start == -1:
            print(f"Domain sequence for {record_id} not found in full-length sequence.")
            continue
        
        #extend boundaries
        extended_start = max(0, full_start - n_extend)
        extended_end = min(len(full_length_seq), full_end + c_extend)
        extended_domain_seq = full_length_seq[extended_start:extended_end]
        
        #save the extended domain
        extended_domains.append((record_id, extended_domain_seq))
    
    #write extended domains to output file
    with open(output_file, "w") as out_f:
        for record_id, seq in extended_domains:
            out_f.write(f">{record_id}_extended\n{seq}\n")
    print(f"Extended domains written to {output_file}")

#example usage
alignment_file = "Twc_fn3_add-fragmetns-mafft.fasta"  #alignment file
output_file = "twc_fn3_starting-point_domains.fasta"  #output file 
extend_domain_boundaries(alignment_file, output_file)

