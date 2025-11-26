Domain Boundary and Sequence Sources
====================================

Here are all the fasta files containing the domains used to make the MSA for the
JalVar alignment. The directory "pre-adjusted_domains" contains the domains before defining the boundaries manually. Aligned to the full protein sequence (referred to as main). These files were used to extend the boundaries so that they could be manually adjusted (see python script "startg_point_domains.py). The "MedakaA-A-band_Ig-Domains.fasta" were only used for the comparative analysis of the titin A-band and not in the JalVar analysis.  

Details
-------------------------------------

Human Titin:
- Initial domain sequences from TitinDB (accessed: 2024-11-20)
- Fibronectin III domains: directly from TitinDB
- Immunoglobulin (Ig) domains: manually using from NCBI RefSeq NP_001254479.2
  (accessed: 2025-07-17)

Obscurin:
- Reference: UniProt Q5VST9
- Domain boundaries: from InterPro (accessed: 2024-11-20)
- Manually refined

Medaka (Oryzias latipes):
- Two titin genes: A and B
- MedakaB titin: XP_023806503.1, domains from InterPro (2024-11-20), manually refined
- MedakaA titin: XP_023806459.1, domains aligned to MedakaB full-length sequence

Twitchin:
- Reference: UniProt Q23551
- Domain boundaries: InterPro + Fong et al. (1996), manually refined

Domain Boundary Manual Adjustment:
- Ig domain boundaries aligned to full-length protein
- Extended by 10–15 residues (N- and C-termini) using Python v3.12.2
- Alignment to known structures using MUSCLE (in Jalview)
  - Parameters: BLOSUM62, max iterations: 16, gap open: -12, gap extension: -1
- Domain boundaries trimmed based on structural alignment
- Realigned to full protein to ensure no overlap between adjacent domains

Naming:
- Human TTN: domain names follow Su et al. (2022)
- Other proteins: domains numbered sequentially based on order in sequence
