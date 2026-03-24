import sys
import os
from Bio import SeqIO

review_file = sys.argv[1]
gbff_file   = sys.argv[2]
outdir      = sys.argv[3]

FLANK = 1500

# -------------------------------
# Load reviewed IDs + optional labels
# -------------------------------
id_map = {}

with open(review_file) as f:
    for line in f:
        if not line.strip():
            continue
        parts = line.strip().split()
        pid   = parts[0]
        label = parts[1] if len(parts) > 1 else None
        id_map[pid] = label

print(f"[INFO] Reviewed IDs: {len(id_map)}")

os.makedirs(outdir, exist_ok=True)

# -------------------------------
# Parse GBFF
# -------------------------------
count = 0

for record in SeqIO.parse(gbff_file, "genbank"):
    genome_seq = record.seq

    for feature in record.features:
        if feature.type != "CDS":
            continue

        qualifiers = feature.qualifiers
        protein_id = (qualifiers.get("protein_id", [None])[0] or
                      qualifiers.get("locus_tag",  [None])[0])

        if protein_id not in id_map:
            continue

        seq_label   = id_map[protein_id]
        base_header = f"{protein_id}|{seq_label}" if seq_label else protein_id
        safe_header = base_header.replace("|", "_")
        count += 1

        g_start = int(feature.location.start)
        g_end   = int(feature.location.end)
        strand  = feature.location.strand

        # -------------------------------
        # CDS: spliced, feature.extract handles strand + joining
        # -------------------------------
        cds_seq = feature.extract(genome_seq)

        # -------------------------------
        # Gene: full unspliced genomic span
        # -------------------------------
        gene_seq = genome_seq[g_start:g_end]

        # -------------------------------
        # Flank: genomic span + FLANK on each side
        # -------------------------------
        f_start   = max(0, g_start - FLANK)
        f_end     = min(len(genome_seq), g_end + FLANK)
        flank_seq = genome_seq[f_start:f_end]

        # -------------------------------
        # Build exon intervals relative to gene start (forward strand coords)
        # feature.location.parts gives individual exon ranges in genomic coords
        # -------------------------------
        gene_len  = g_end - g_start
        flank_len = f_end - f_start

        # Each part is one exon; collect (start, end) relative to the gene's
        # forward-strand origin at g_start.
        exon_coords_gene  = []   # relative to gene sequence start
        exon_coords_flank = []   # relative to flank sequence start

        for part in feature.location.parts:
            ps = int(part.start)
            pe = int(part.end)

            # Relative to gene start
            rel_s = ps - g_start
            rel_e = pe - g_start

            # Relative to flank start
            fl_s = ps - f_start
            fl_e = pe - f_start

            exon_coords_gene.append((rel_s, rel_e))
            exon_coords_flank.append((fl_s, fl_e))

        # -------------------------------
        # RC: flip gene and flank sequences.
        # Also flip the exon intervals to match the RC'd coordinate space.
        # After RC, position p in the original becomes (length - 1 - p),
        # so interval (s, e) becomes (length - e, length - s).
        # -------------------------------
        if strand == -1:
            gene_seq  = gene_seq.reverse_complement()
            flank_seq = flank_seq.reverse_complement()

            exon_coords_gene = [
                (gene_len - e, gene_len - s)
                for s, e in exon_coords_gene
            ]
            exon_coords_flank = [
                (flank_len - e, flank_len - s)
                for s, e in exon_coords_flank
            ]

        # Sort intervals so the coords file reads 5'→3'
        exon_coords_gene.sort()
        exon_coords_flank.sort()

        # -------------------------------
        # Write FASTA files
        # -------------------------------
        with open(f"{outdir}/{safe_header}_cds.fasta", "w") as out:
            out.write(f">{base_header}_cds\n{cds_seq}\n")

        with open(f"{outdir}/{safe_header}_gene.fasta", "w") as out:
            out.write(f">{base_header}_gene\n{gene_seq}\n")

        with open(f"{outdir}/{safe_header}_flank.fasta", "w") as out:
            out.write(f">{base_header}_flank\n{flank_seq}\n")

        # -------------------------------
        # Write coords files
        # Format: one exon per line, tab-separated start and end
        # Coords are 0-based, half-open [start, end)
        # -------------------------------
        with open(f"{outdir}/{safe_header}_gene.coords", "w") as out:
            for s, e in exon_coords_gene:
                out.write(f"{s}\t{e}\n")

        with open(f"{outdir}/{safe_header}_flank.coords", "w") as out:
            for s, e in exon_coords_flank:
                out.write(f"{s}\t{e}\n")

print(f"[INFO] Extracted {count} loci")