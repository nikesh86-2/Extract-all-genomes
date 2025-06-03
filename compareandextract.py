import argparse
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
import subprocess
import os
import tempfile

def parse_tbl(tbl_file):
    hits = []
    with open(tbl_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split()
            try:
                seq_id = fields[0]
                start = int(fields[7])
                end = int(fields[8])
                strand = fields[9]
                hits.append((seq_id, start, end, strand))
            except (IndexError, ValueError):
                print(f"[!] Skipping line: {line.strip()}")
    return hits

def extract_hits(hits, genome_fasta, extracted_fa):
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    records = []

    for i, (seq_id, start, end, strand) in enumerate(hits):
        if seq_id not in genome:
            print(f"[!] Sequence ID {seq_id} not found.")
            continue
        full_seq = genome[seq_id].seq
        subseq = full_seq[start-1:end]
        if strand == "-":
            subseq = subseq.reverse_complement()
        records.append(SeqRecord(subseq, id=f"{seq_id}_{start}_{end}_{strand}", description=""))

    SeqIO.write(records, extracted_fa, "fasta")
    print(f"[✓] Wrote {len(records)} sequences to {extracted_fa}")
    return extracted_fa

def run_clustal(infile, outfile):
    print(f"[⏳] Running Clustal Omega...")
    subprocess.run([
        "clustalo", "-i", infile, "-o", outfile,
        "--outfmt=fasta", "--force", "--seqtype=DNA", "--threads=8"
    ], check=True)
    print(f"[✓] Alignment complete: {outfile}")

def plot_logo(aligned_fasta, output_png="sequence_logo.svg"):
    alignment = AlignIO.read(aligned_fasta, "fasta")
    aln_len = alignment.get_alignment_length()

    pfm = []
    for i in range(aln_len):
        column = [record.seq[i].upper() for record in alignment]
        counts = Counter(column)
        row = {nt: counts.get(nt, 0) for nt in "ACGT"}
        pfm.append(row)

    df = pd.DataFrame(pfm)
    df.index.name = "position"

    # Define color scheme by base
    color_scheme = {
        'A': '#64b964',  # green
        'C': '#4c9be8',  # blue
        'G': '#f5b942',  # orange
        'T': '#e04f5f',  # red
    }

    logo = logomaker.Logo(df, color_scheme=color_scheme, font_name="Arial")
    logo.style_spines(visible=False)
    logo.ax.set_title("Sequence Logo of Extracted Hits")
    logo.ax.set_xlabel("Position")
    logo.ax.set_ylabel("Count")

    plt.tight_layout()
    plt.savefig(output_png)
    print(f"[✓] Logo saved as {output_png}")
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Extract, align, and visualize RNA hits.")
    parser.add_argument("tbl", help="Path to results.tbl from cmsearch")
    parser.add_argument("genome", help="FASTA genome file")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        extracted = os.path.join(tmpdir, "hits.fa")
        aligned = os.path.join(tmpdir, "aligned.fa")

        hits = parse_tbl(args.tbl)
        extract_hits(hits, args.genome, extracted)
        run_clustal(extracted, aligned)
        plot_logo(aligned)

if __name__ == "__main__":
    main()
