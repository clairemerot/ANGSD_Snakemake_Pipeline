import sys
import random

def main():
    # --- 1. Argument Parsing ---
    # When called via shell, we retrieve arguments from the command line
    # sys.argv[1] = input .fai file
    # sys.argv[2] = output file path
    # sys.argv[3] = total number of sites to sample
    # sys.argv[4] = comma-separated string of chromosomes
    
    try:
        input_fai    = sys.argv[1]
        output_sites = sys.argv[2]
        total_sites  = int(sys.argv[3])
        chroms_list  = sys.argv[4].split(",")
    except IndexError:
        print("Usage: python3 select_sites.py <input.fai> <output.sites> <total_sites> <chrom1,chrom2,...>")
        sys.exit(1)

    # --- 2. Read Reference Index (.fai) ---
    # The .fai format is: Chromosome_Name \t Length \t ...
    with open(input_fai, "r") as f:
        lines = [line.strip().split("\t") for line in f if line.strip()]

    # --- 3. Filter and Calculate Total Genome Length ---
    # We only consider chromosomes specified in the 'chroms_list'
    total_length = 0
    valid_chroms = []
    
    for sl in lines:
        chrom = sl[0]
        length = int(sl[1])
        if chrom in chroms_list:
            total_length += length
            valid_chroms.append((chrom, length))

    # --- 4. Proportional Sampling ---
    # We set a seed for reproducibility (crucial for bioinformatics pipelines)
    random.seed(42)
    
    with open(output_sites, "w") as fout:
        for chrom, length in valid_chroms:
            # Calculate how many sites to pick for this specific chromosome
            # based on its proportion of the total selected genome length
            theoretical_nb = round(total_sites * (length / total_length))
            
            # SAFETY CHECK: 
            # We cannot sample more sites than the chromosome length.
            # min(theoretical_nb, length) ensures we stay within bounds.
            nb_to_pick = min(theoretical_nb, length)

            if nb_to_pick > 0:
                # random.sample picks unique integers without replacement.
                # range(1, length + 1) covers 1-based genomic coordinates.
                # Sorting ensures the output file is ordered by position.
                sampled_positions = sorted(random.sample(range(1, length + 1), nb_to_pick))
                
                # Write results to the output file in Tab-Separated format
                for pos in sampled_positions:
                    fout.write(f"{chrom}\t{pos}\n")

if __name__ == "__main__":
    main()
