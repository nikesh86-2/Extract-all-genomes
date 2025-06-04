from Bio import Entrez
import time
import tkinter as tk
from tkinter import simpledialog, messagebox

# --- GUI Setup ---
root = tk.Tk()
root.withdraw()  # Hide main window

# Prompt for email
user_email = simpledialog.askstring(
    "Input", "Enter user email:"
)

if not user_email or not user_email.strip():
    messagebox.showerror("Error", "No query sequence entered. Exiting.")
    exit()

user_email = user_email.strip().upper()


# Always include your email or ncbi will limit your search
Entrez.email = user_email

# Prompt for organism
user_organism = simpledialog.askstring(
    "Input", "Enter target organism:"
)
user_organism = user_organism.strip().title()

if not user_organism or not user_organism.strip():
    messagebox.showerror("Error", "No query sequence entered. Exiting.")
    exit()

# Always include your email or ncbi will limit your search
Entrez.email = user_email

# Search term for full-length genomes (adjust as needed)
search_term = f'"{user_organism}"[Organism] AND "complete genome"[Title]' 
batch_size = 1000

def fetch_full_genomes(term, batch_size=1000): #adjust for any organism?
    # Step 1: Search NCBI for matching sequences
    search_handle = Entrez.esearch(
        db="nucleotide",
        term=term,
        usehistory="y"  # Enables use of WebEnv for batch downloads
    )
    search_results = Entrez.read(search_handle)
    count = int(search_results["Count"])
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]

    print(f"Total sequences found: {count}")

    if count == 0:
        messagebox.showinfo("No Results", "No genome sequences found for this organism.")
        return []

    all_data = []

    # Step 2: Fetch in chunks using retstart/retmax
    for start in range(0, count, batch_size):
        print(f"Fetching records {start} to {min(start + batch_size, count)}...")
        fetch_handle = Entrez.efetch(
            db="nucleotide",
            rettype="fasta",  # or "gb" if you want gb format
            retmode="text",
            retstart=start,
            retmax=batch_size,
            webenv=webenv,
            query_key=query_key
        )
        data = fetch_handle.read()
        all_data.append(data)
        fetch_handle.close()
        time.sleep(0.5)  # careful not to upset the NCBI servers
# have to go in batches, otherwise you get limited to 500
    return all_data

# Run the function and save to a file
results = fetch_full_genomes(search_term)

with open(str(user_organism) + ".fasta", "w") as out_handle:
    out_handle.write("".join(results))

print("Download complete! Sequences saved to " + str(user_organism) + ".fasta")
