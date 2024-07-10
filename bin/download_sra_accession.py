#!/usr/bin/env python
import subprocess
import os
import argparse


def fetch_data(sra_id):
    # Execute the command (prefetching) and capture its output
    prefetch_cmd = ["prefetch", sra_id]
    print(f"\nRunning: {' '.join(prefetch_cmd)}")
    prefetch = subprocess.Popen(prefetch_cmd, stdout=subprocess.PIPE)
    output, _ = prefetch.communicate()

    # Convert the output bytes to a string
    output = output.decode('utf-8')
    print(output)

    # Split the output into words
    words = output.split()

    # Find the word that starts with "SRR"
    srr_id = None
    for word in words:
        if word.startswith("'SRR"):
            srr_id = word.split("'")[1]
            break

    # fetch for data
    if srr_id:
        # Now download the fastq and capture its output
        fetch_cmd = ["fastq-dump", "--outdir", "raw", "--gzip", "--skip-technical",  "--readids", "--read-filter", "pass", "--dumpbase", "--split-3", "--clip", f"{srr_id}/{srr_id}.sra"]
        print(f"\nRunning: {' '.join(fetch_cmd)}")
        fetch = subprocess.Popen(fetch_cmd, stdout=subprocess.PIPE)
        output, _ = fetch.communicate()

        # Convert the output bytes to a string
        output = output.decode('utf-8')
        print(output)
    else:
        print("No word starting with 'SRR' found.")


def subsample(sra_id, paired, coverage):
    if paired:
        subsample_cmd = ["rasusa", "--bases", coverage, "-i", f"raw/{sra_id}*.fastq.gz", "-o", f"subsampled/{sra_id}_1.fastq.gz", f"subsampled/{sra_id}_2.fastq.gz"]
    else:
        subsample_cmd = ["rasusa", "--bases", coverage, "-i", f"raw/{sra_id}*.fastq.gz", "-o", f"subsampled/{sra_id}.fastq.gz"]

    print(f"\nRunning: {' '.join(subsample_cmd)}")
    subsample = subprocess.Popen(subsample_cmd, stdout=subprocess.PIPE)
    output, _ = subsample.communicate()

    # Convert the output bytes to a string
    output = output.decode('utf-8')
    print(output)

    os.remove(f"raw/{sra_id}*")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download SRA accession and subsample.')
    parser.add_argument("--accession", dest='accession', help='SRA accession')
    parser.add_argument('--paired', action="store_true", help='Should rasusa use paired mode?')
    parser.add_argument("--coverage", dest='coverage', help='Read coverage to return')

    args = parser.parse_args()

    print(f"Fetching {args.accession}...")
    fetch_data(args.accession)
    os.remove(f"{args.accession}/{args.accession}.sra")
    os.rmdir(args.accession)
    subsample(args.accession, args.paired, args.coverage)
