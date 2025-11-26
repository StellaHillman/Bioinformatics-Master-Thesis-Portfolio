#this script was made to get number of shared variants for Venn diagram
#there are also different tools that do the same thing and I used these to double check e.g. vcf-compare(v0.1.14-12-gcdb80b8)

from cyvcf2 import VCF

def extract_variants(file):
    return set((v.CHROM, v.POS, v.REF, tuple(v.ALT)) for v in VCF(file))

a = extract_variants("TTN-ABraOM-mSNVs.vcf.gz")
b = extract_variants("TTN-ClinVar-SNVs.vcf.gz")
c = extract_variants("TTN-gnomAD_mSNVs.sorted.vcf.gz")

only_a = a - b - c
only_b = b - a - c
only_c = c - a - b
a_b = (a & b) - c
a_c = (a & c) - b
b_c = (b & c) - a
abc  = a & b & c

print(f"A only: {len(only_a)}")
print(f"B only: {len(only_b)}")
print(f"C only: {len(only_c)}")
print(f"A & B: {len(a_b)}")
print(f"A & C: {len(a_c)}")
print(f"B & C: {len(b_c)}")
print(f"All three: {len(abc)}")

