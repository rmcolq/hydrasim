#!/usr/bin/env python


from Bio import Entrez
import sys

fasta_handle = Entrez.efetch(
        db="nucleotide", id=gi, rettype="fasta", retmode="text"
    )

with open(f"{sys.argv[1]}_genomic.fna", "wt") as fh:
          for line in fasta_handle:
            fh.write(line)