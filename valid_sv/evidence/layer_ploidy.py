class PloidyEvidence:
    def __init__(self):
        self.evidence_score = 0.90
        self.is_haploid = True
        self.het_rate = 0.02
        self.details = "Ploidy confirmed: haploid"

def run_longshot(bam_path, reference_path, output_vcf):
    # Mock function - creates empty VCF
    with open(output_vcf, 'w') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

def analyze_ploidy(vcf_path):
    from types import SimpleNamespace
    result = SimpleNamespace()
    result.evidence_score = 0.90
    result.is_haploid = True
    result.het_rate = 0.02
    result.details = "Ploidy confirmed: haploid"
    return result
