# ==============================
# CONFIG
# ==============================
INPUT_FASTA = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz2_run2/nr_run/ftsz2.nr.fasta"
OUTPUT_FASTA = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz2_run2/nr_run/ftsz2_cterm.nr.fasta"
N = 50  # number of residues to extract from the end


# ==============================
# FUNCTION
# ==============================
def extract_last_n_residues(input_fasta, output_fasta, n):
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        
        header = None
        seq_lines = []

        for line in infile:
            line = line.strip()

            if line.startswith(">"):
                # Process previous sequence
                if header is not None:
                    sequence = "".join(seq_lines)
                    last_n = sequence[-n:] if len(sequence) >= n else sequence
                    
                    outfile.write(header + "\n")
                    outfile.write(last_n + "\n")

                # Start new sequence
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        # Process last sequence
        if header is not None:
            sequence = "".join(seq_lines)
            last_n = sequence[-n:] if len(sequence) >= n else sequence
            
            outfile.write(header + "\n")
            outfile.write(last_n + "\n")


# ==============================
# RUN
# ==============================
if __name__ == "__main__":
    extract_last_n_residues(INPUT_FASTA, OUTPUT_FASTA, N)
    print(f"Done! Output written to {OUTPUT_FASTA}")