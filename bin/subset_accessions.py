#!/usr/bin/env python

import csv
import random
import argparse
from collections import defaultdict
import urllib2

def select_random_rows(input_csv, output_csv, category_column, sample_size):
    # check  if input_csv is URL (from testing profile)
    if input_csv.startswith('https://'):
        input_csv = urllib2.urlopen(input_csv)
    # Read the CSV file
    grouped_data = defaultdict(list)
    with open(input_csv, mode='r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            grouped_data[row[category_column]].append(row)
    
    # Randomly select 10 rows from each group
    selected_rows = []
    for category, rows in grouped_data.items():
        category_rows = []
        print (len(rows), sample_size, type(rows))
        num_without_replacement = min(len(rows), sample_size)
        chosen = random.sample(range(len(rows)), num_without_replacement)
        category_rows.extend([rows[i] for i in chosen])

        num_with_replacement = sample_size - num_without_replacement
        if num_with_replacement > 0:
            chosen = random.choices(range(len(rows)), k=num_with_replacement)
            category_rows.extend([rows[i].copy() for i in chosen])

        # add an index
        i = 0
        for row in category_rows:
            row["index"] = i
            i+=1

        selected_rows.extend(category_rows)

    # Write the selected rows to a new CSV file
    with open(output_csv, mode='w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames + ["index"])
        writer.writeheader()
        writer.writerows(selected_rows)

    for row in selected_rows:
        print(row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Randomly select n rows from each group in the category column of a CSV file.')
    parser.add_argument('input_csv', help='Input CSV file path')
    parser.add_argument('output_csv', help='Output CSV file path')
    parser.add_argument('category_column', help='Category column to group by')
    parser.add_argument('sample_size', type=int, help="Size of random sample")

    args = parser.parse_args()
    select_random_rows(args.input_csv, args.output_csv, args.category_column, args.sample_size)

