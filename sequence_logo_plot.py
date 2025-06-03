import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, AlignIO
import logomaker
import matplotlib.pyplot as plt
import os
import tkinter as tk
from tkinter import simpledialog, messagebox

# --- GUI Setup ---
root = tk.Tk()
root.withdraw()  # Hide main window

# Prompt for query sequence
query_sequence = simpledialog.askstring(
    "Input", "Enter the query nucleotide sequence (A/C/G/T only):"
)

if not query_sequence or not query_sequence.strip():
    messagebox.showerror("Error", "No query sequence entered. Exiting.")
    exit()

query_sequence = query_sequence.strip().upper()

# Prompt for output PNG filename
output_png = simpledialog.askstring(
    "Input", "Enter output PNG filename (e.g., sequence_logo.png):", initialvalue="sequence_logo.png"
)

if not output_png or not output_png.strip():
    messagebox.showerror("Error", "No output filename entered. Exiting.")
    exit()

output_png = output_png.strip()

# Check for existing file
if os.path.exists(output_png):
    if not messagebox.askyesno("Overwrite?", f"'{output_png}' already exists. Overwrite?"):
        messagebox.showinfo("Aborted", "Operation cancelled by user.")
        exit()

# save query to fasta - important for downstream processing (clustal)
with open("query.fasta", "w") as f:
    f.write(">query\n" + query_sequence + "\n")

#blast and extract sequences
# Output: sseqid and the full matching subject sequence
result = subprocess.run([
    "blastn",
    "-query", "query.fasta",
    "-db", "hbv_db",
    "-outfmt", "6 sseqid sseq",
    "-perc_identity", "80",
    "-evalue", "0.001",
    "-num_alignments", "1000000",  # High max
    "-max_hsps", "1",               # One HSP per hit
    "-out", "blast_ids.txt"
], capture_output=True, text=True)

# convert the BLAST output to sequences ===
input_txt = "blast_ids.txt"
query_id = "query"

# convert the TSV file to FASTA records
records = []
with open(input_txt) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) != 2:
            continue
        seq_id, seq = parts
        records.append(SeqRecord(Seq(seq), id=seq_id, description=""))

records.append(SeqRecord(Seq(query_sequence), id=query_id, description=""))

#Save to FASTA
SeqIO.write(records, "matches.fasta", "fasta")

if len(records) < 2:
    messagebox.showerror("Error", "Too few sequences to align. Exiting.")
    exit()

#Run Clustal Omega
subprocess.run([
    "clustalo",
    "-i", "matches.fasta",
    "-o", "aligned.fasta",
    "--force",
    "--threads=8"
])
 #load alignment
alignment = AlignIO.read("aligned.fasta", "fasta")

# find alignment within sequence
query_rec = next((rec for rec in alignment if query_id in rec.id.lower()), None)
if query_rec is None:
    messagebox.showerror("Error", "Query sequence not found in alignment. Exiting.")
    exit()

query_str = str(query_rec.seq)

# Find start of query within alignment
query_clean = query_sequence.replace("-", "")
start_idx = query_str.find(query_clean[:10])
if start_idx == -1:
    messagebox.showerror("Error", "Unable to locate query in alignment. Exiting.")
    exit()

end_idx = start_idx + len(query_clean)

#Extract alignment region from all records
window = []
for rec in alignment:
    subseq = str(rec.seq[start_idx:end_idx])
    if "-" not in subseq:
        window.append(subseq)

if len(window) < 2:
    messagebox.showerror("Error", "Not enough aligned sequences without gaps. Exiting.")
    exit()

#Plot sequence logo
counts_df = logomaker.alignment_to_matrix(window)

plt.figure(figsize=(12, 4))
logo = logomaker.Logo(
    counts_df,
    color_scheme="classic",  # Nucleotide coloring (A: green, C: blue, G: orange, T: red)
    shade_below=0.5,
    fade_below=0.5,
    font_name="Arial Rounded MT Bold"
)
plt.title("Sequence Logo (Query Region)")
plt.xlabel("Position in Query")
plt.ylabel("Information Content")
plt.tight_layout()
plt.savefig(output_png, dpi=300)
plt.show()

messagebox.showinfo("Done", f"Logo saved to: {output_png}")
