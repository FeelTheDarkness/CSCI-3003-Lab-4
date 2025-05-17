import sys

def read_fasta(filename):
    """Read FASTA file and return lists of headers and sequences"""
    headers = []
    sequences = []
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                headers.append(line[1:])
                if current_seq:
                    sequences.append(''.join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line)
        if current_seq:
            sequences.append(''.join(current_seq))
    
    return headers, sequences

def read_vcf(filename):
    """Read VCF file and return lists of SNP information"""
    chrom = []
    snp_id = []
    slim_pos = []
    ref_snp = []
    alt_snp = []
    gene_info = []
    disease_info = []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue  # skip header lines
            parts = line.strip().split('\t')
            chrom.append(parts[0])
            snp_id.append(parts[1])
            slim_pos.append(int(parts[2]))
            ref_snp.append(parts[3])
            alt_snp.append(parts[4])
            gene_info.append(parts[5])
            disease_info.append(parts[6])
    
    return chrom, snp_id, slim_pos, ref_snp, alt_snp, gene_info, disease_info

def calculate_frequencies(sequences, slim_pos, ref_snp, alt_snp):
    """Calculate SNP frequencies in the population"""
    frequencies = []
    num_individuals = len(sequences)
    
    for i in range(len(slim_pos)):
        pos = slim_pos[i] - 1  # convert to 0-based index
        ref = ref_snp[i]
        alt = alt_snp[i]
        
        count = 0
        for seq in sequences:
            if seq[pos] == alt:
                count += 1
        
        frequency = count / num_individuals
        frequencies.append(frequency)
    
    return frequencies

def find_extreme_frequencies(snp_id, frequencies):
    """Find SNPs with highest and lowest frequencies"""
    max_freq = max(frequencies)
    min_freq = min(frequencies)
    
    max_snps = [snp_id[i] for i, freq in enumerate(frequencies) if freq == max_freq]
    min_snps = [snp_id[i] for i, freq in enumerate(frequencies) if freq == min_freq]
    
    return max_snps, min_snps

def write_output(filename, chrom, snp_id, slim_pos, ref_snp, alt_snp, frequencies, gene_info, disease_info):
    """Write SNP information with frequencies to output file"""
    with open(filename, 'w') as f:
        # Write header
        f.write("Chromosome\tSNP Name\tSlim Position\tReference SNP\tAlternative SNP\tSNP Frequency\tGene Disease Information\n")
        
        for i in range(len(snp_id)):
            line = f"{chrom[i]}\t{snp_id[i]}\t{slim_pos[i]}\t{ref_snp[i]}\t{alt_snp[i]}\t{frequencies[i]}\t{disease_info[i]}\n"
            f.write(line)

def mean_frequency(frequencies):
    """Calculate mean of SNP frequencies (optional)"""
    return sum(frequencies) / len(frequencies)

def find_lymphoma_variants(headers, sequences, snp_id, slim_pos, disease_info):
    """Find variants associated with Non-Hodgkin lymphoma"""
    with open('lymphoma_variants.txt', 'w') as f:
        f.write("IndividualID\tVariantID\tGeneDiseaseInfo\n")
        
        for i in range(len(snp_id)):
            if "Non-Hodgkin lymphoma" in disease_info[i]:
                pos = slim_pos[i] - 1
                alt = alt_snp[i]
                
                for j in range(len(headers)):
                    if sequences[j][pos] == alt:
                        line = f"{headers[j]}\t{snp_id[i]}\t{disease_info[i]}\n"
                        f.write(line)

def find_novel_variants(headers, sequences, ref_sequence):
    """Discover novel variants not in the VCF file"""
    variant_counts = {}  # {position: {'ref': base, 'alts': {base: count}}}
    
    # Initialize with reference sequence
    for i in range(len(ref_sequence)):
        variant_counts[i+1] = {'ref': ref_sequence[i], 'alts': {}}
    
    # Scan all individual sequences
    for seq in sequences:
        for i in range(len(seq)):
            pos = i + 1  # 1-based position
            if seq[i] != variant_counts[pos]['ref']:
                alt = seq[i]
                variant_counts[pos]['alts'][alt] = variant_counts[pos]['alts'].get(alt, 0) + 1
    
    # Filter for positions with variants and calculate frequencies
    novel_variants = []
    num_individuals = len(sequences)
    
    for pos in variant_counts:
        if variant_counts[pos]['alts']:
            ref = variant_counts[pos]['ref']
            for alt in variant_counts[pos]['alts']:
                count = variant_counts[pos]['alts'][alt]
                frequency = count / num_individuals
                novel_variants.append((pos, ref, alt, frequency))
    
    # Write to file
    with open('novel_variants.txt', 'w') as f:
        f.write("Position on Chr2\tReference base\tAlternative base\tSNP frequency\n")
        for variant in novel_variants:
            f.write(f"{variant[0]}\t{variant[1]}\t{variant[2]}\t{variant[3]}\n")
    
    return novel_variants

def main():
    # Read input files
    headers, sequences = read_fasta('slim_chr2_seq.fasta')
    chrom, snp_id, slim_pos, ref_snp, alt_snp, gene_info, disease_info = read_vcf('slim_chr2_SNPS.vcf')
    
    # Calculate frequencies
    frequencies = calculate_frequencies(sequences, slim_pos, ref_snp, alt_snp)
    
    # Find extreme frequencies
    max_snps, min_snps = find_extreme_frequencies(snp_id, frequencies)
    print("SNPs with highest frequency:", ", ".join(max_snps))
    print("SNPs with lowest frequency:", ", ".join(min_snps))
    
    # Write output file
    write_output('slim_chr2_SNPs_withfrequencies.vcf', chrom, snp_id, slim_pos, 
                ref_snp, alt_snp, frequencies, gene_info, disease_info)
    
    # Optional: Calculate mean frequency
    print(f"Mean SNP frequency: {mean_frequency(frequencies):.4f}")
    
    # Part II: Find Non-Hodgkin lymphoma variants
    find_lymphoma_variants(headers, sequences, snp_id, slim_pos, disease_info)
    
    # Part III: Find novel variants
    # Assuming the first sequence is the reference sequence
    ref_sequence = sequences[0]
    novel_variants = find_novel_variants(headers[1:], sequences[1:], ref_sequence)

if __name__ == "__main__":
    main()