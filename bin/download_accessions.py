#!/usr/bin/env python

from Bio import Entrez, SeqIO
import sys

def download_fasta(accession, email):
    Entrez.email = email  # Always tell NCBI who you are
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    fasta_record = handle.read()
    handle.close()

    # Save the FASTA record to a file
    with open(f"{accession}.fasta", "w") as file:
        file.write(fasta_record)

    print(f"FASTA file for accession {accession} has been saved as {accession}.fasta")

# Example usage
print(sys.argv[1])
print(sys.argv[2])
# download_fasta(sys.argv[1], sys.argv[2])


# from Bio import Entrez
# import sys

# Entrez.email = 'nfellaby@gmail.com'

# # fasta_handle = Entrez.efetch(
# #         db="nucleotide", id=str('"')+str(sys.argv[1])+str('"'), rettype="fasta", retmode="text"
# #     )


# handle = Entrez.efetch(db="nucleotide", id=sys.argv[1], rettype="fasta", retmode="text")

# with open(f"{sys.argv[1]}_genomic.fna", "wt") as fh:
#           for line in handle:
#             fh.write(line)