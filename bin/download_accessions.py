#!/usr/bin/env python

from Bio import Entrez, SeqIO
import sys
import os
import time

def download_fasta(accession, email):
    Entrez.email = email  # Always tell NCBI who you are
    try:

        # handle = Entrez.efetch(db='nuccore', id='NC_019843.3', format='fasta', rettype='fasta')
            
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")
        fasta_record = handle.read()
        handle.close()

        # Define the file name
        file_name = f"{accession}.fasta"

        # Save the FASTA record to a file
        with open(file_name, "w") as file:
            file.write(fasta_record)
        
        # Check if the file is created
        while not os.path.exists(file_name):
            print(f"Waiting for {file_name} to be created...")
            time.sleep(1)

        print(f"FASTA file for accession {accession} has been saved as {file_name}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python download_accessions.py <accession> <email>")
    else:
        download_fasta(sys.argv[1], sys.argv[2])