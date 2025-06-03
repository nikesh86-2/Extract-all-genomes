from Bio import Entrez
import time

# Always include your email or ncbi will limit your search
Entrez.email = "fbsnpat@leeds.ac.uk"

# Search term for full-length genomes (adjust as needed)
search_term = '"Hepatitis B virus"[Organism] AND "complete genome"[Title]' # rewrite so that prompt for organism?
batch_size = 1000

def fetch_full_hbv_genomes(term, batch_size=1000): #adjust for any organism?
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
results = fetch_full_hbv_genomes(search_term)

with open("hbv_full_genomes.fasta", "w") as out_handle:
    out_handle.write("".join(results))

print("Download complete! Sequences saved to hbv_full_genomes.fasta")

#if rewrite for user prompt, these lines will need adjusting (47-52)