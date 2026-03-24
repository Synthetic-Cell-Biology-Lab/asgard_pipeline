import sys
import os
import glob
from Bio import SeqIO
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch

# -------------------------------
# Inputs
# -------------------------------
seq_dir    = sys.argv[1]
output_pdf = sys.argv[2]

LINE_WIDTH = 80

# -------------------------------
# Load FASTA records from a glob pattern
# Strips a known suffix from rec.id to produce the lookup key
# -------------------------------
def load_fasta_records(pattern, strip_suffix):
    records = {}
    for f in glob.glob(pattern):
        rec = next(SeqIO.parse(f, "fasta"))
        key = rec.id
        if key.endswith(strip_suffix):
            key = key[: -len(strip_suffix)]
        records[key] = rec
    return records

gene_records  = load_fasta_records(f"{seq_dir}/*_gene.fasta",  "_gene")
cds_records   = load_fasta_records(f"{seq_dir}/*_cds.fasta",   "_cds")
flank_records = load_fasta_records(f"{seq_dir}/*_flank.fasta", "_flank")

print(f"[INFO] Loaded:")
print(f"  Gene:  {len(gene_records)}")
print(f"  CDS:   {len(cds_records)}")
print(f"  Flank: {len(flank_records)}")

# -------------------------------
# Load exon coordinate intervals from a .coords file
# Returns a sorted list of (start, end) tuples, or None if missing
# Coords are 0-based half-open [start, end)
# -------------------------------
def load_coords(coords_path):
    if not os.path.exists(coords_path):
        return None
    intervals = []
    with open(coords_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            s, e = line.split("\t")
            intervals.append((int(s), int(e)))
    intervals.sort()
    return intervals

# -------------------------------
# Build a boolean lookup array for fast per-position CDS testing
# Much faster than calling any() over intervals for every base
# -------------------------------
def build_cds_mask(seq_len, intervals):
    mask = bytearray(seq_len)   # 0 = not CDS, 1 = CDS
    if intervals:
        for s, e in intervals:
            s = max(0, s)
            e = min(seq_len, e)
            for i in range(s, e):
                mask[i] = 1
    return mask

# -------------------------------
# Collapse consecutive same-color bases into a single run.
# Emitting one <font> tag per base produces enormous XML that
# slows ReportLab significantly for long sequences.
# -------------------------------
def build_colored_lines(seq_str, cds_mask, non_cds_color, line_width=LINE_WIDTH):
    """
    seq_str       : plain DNA string
    cds_mask      : bytearray, 1 where position is CDS
    non_cds_color : "black" for gene (introns), "red" for flanks
    Returns a list of strings, one per line, with ReportLab XML markup.
    """
    lines = []
    n = len(seq_str)

    for i in range(0, n, line_width):
        chunk     = seq_str[i : i + line_width]
        chunk_len = len(chunk)

        parts = []
        run_start = 0

        while run_start < chunk_len:
            is_cds = bool(cds_mask[i + run_start])
            run_end = run_start + 1

            # Extend the run while color stays the same
            while run_end < chunk_len and bool(cds_mask[i + run_end]) == is_cds:
                run_end += 1

            bases = chunk[run_start:run_end]

            if is_cds:
                parts.append(f'<font color="blue">{bases}</font>')
            elif non_cds_color == "red":
                parts.append(f'<font color="red">{bases}</font>')
            else:
                parts.append(bases)   # plain black — no tag needed

            run_start = run_end

        lines.append("".join(parts))

    return lines

# -------------------------------
# PDF setup
# -------------------------------
doc    = SimpleDocTemplate(output_pdf, pagesize=letter)
styles = getSampleStyleSheet()

mono_style = ParagraphStyle(
    "SeqMono",
    parent   = styles["Normal"],
    fontName = "Courier",
    fontSize = 8,
    leading  = 10,
    spaceAfter = 0,
)

header_style = ParagraphStyle(
    "SeqHeader",
    parent     = styles["Normal"],
    fontName   = "Helvetica-Bold",
    fontSize   = 9,
    leading    = 12,
    spaceAfter = 2,
)

elements = []

# -------------------------------
# Helper: add a labelled sequence block to the PDF
# -------------------------------
def add_sequence_block(header, lines):
    elements.append(Paragraph(f"&gt;{header}", header_style))
    elements.append(Spacer(1, 0.05 * inch))
    for line in lines:
        elements.append(Paragraph(line, mono_style))
    elements.append(Spacer(1, 0.25 * inch))

# -------------------------------
# Generate report
# -------------------------------
all_ids = sorted(set(gene_records) | set(flank_records))

for pid in all_ids:

    # Derive the safe filename prefix (| replaced by _ in filenames)
    safe_pid = pid.replace("|", "_")

    # --------------------------------------------------
    # Gene sequence
    # CDS positions (exons) → blue
    # Non-CDS positions (introns/UTRs) → black
    # --------------------------------------------------
    if pid in gene_records:
        gene     = gene_records[pid]
        gene_str = str(gene.seq)
        seq_len  = len(gene_str)

        coords_path = os.path.join(seq_dir, f"{safe_pid}_gene.coords")
        intervals   = load_coords(coords_path)

        if intervals is None:
            print(f"[WARN] No coords file for {pid} gene — rendering all black")

        mask  = build_cds_mask(seq_len, intervals or [])
        lines = build_colored_lines(gene_str, mask, non_cds_color="black")
        add_sequence_block(f"{pid}_gene", lines)

    # --------------------------------------------------
    # Flank sequence
    # CDS positions (exons present in flank) → blue
    # Non-CDS positions (flanking genomic DNA) → red
    # --------------------------------------------------
    if pid in flank_records:
        flank     = flank_records[pid]
        flank_str = str(flank.seq)
        seq_len   = len(flank_str)

        coords_path = os.path.join(seq_dir, f"{safe_pid}_flank.coords")
        intervals   = load_coords(coords_path)

        if intervals is None:
            print(f"[WARN] No coords file for {pid} flank — rendering all red")

        mask  = build_cds_mask(seq_len, intervals or [])
        lines = build_colored_lines(flank_str, mask, non_cds_color="red")
        add_sequence_block(f"{pid}_flank", lines)

# -------------------------------
# Build PDF
# -------------------------------
doc.build(elements)
print(f"[INFO] PDF report generated: {output_pdf}")