#!/usr/bin/env python

import sys
import pandas as pd
import csv

# Get the input file and gem_id from command-line arguments
input_file = sys.argv[1]
gem_id = sys.argv[2]

# Specify the output file path and name
output_file_path = "/home/groups/singlecell/mabdalfttah/projects/DOLSORI_05/jobs/{}/config.csv".format(gem_id)

# Read the CSV file into a DataFrame
df = pd.read_csv(input_file)

# Filter the DataFrame based on the gem_id value
filtered_df = df[df['gem_id'] == gem_id]

# Update the 'fastq_id' column with the gem_id value
filtered_df['fastq_id'] = gem_id

# Data for the [gene-expression] section
gene_expression_data = [
    ['[gene-expression]'],
    ['reference', '/scratch/groups/singlecell/data/reference/refdata-gex-GRCh38-2020-A'],
    ['cmo-set', '/home/groups/singlecell/mabdalfttah/projects/data/CMO_reference.csv'],
    ['expect-cells', '30000'],
    ['chemistry', 'SC3Pv3'],
    ['no-secondary', 'true'],
    ['no-bam', 'false']
]

# Data for the [libraries] section
libraries_data = [
    ['[libraries]'],
    ['fastq_id', 'fastqs', 'feature_types'],
    [gem_id, '/home/groups/singlecell/mabdalfttah/projects/DOLSORI_06/jobs/{}/fastq'.format(gem_id), 'Multiplexing Capture'],
    [gem_id, '/home/groups/singlecell/mabdalfttah/projects/DOLSORI_05/jobs/{}/fastq'.format(gem_id), 'Gene Expression']
]

# Data for the [samples] section
samples_data = [
    ['[samples]'],
    ['sample_id', 'cmo_ids'],
] + filtered_df[['sample_id', 'CMO_id']].values.tolist()

# Open the file in write mode
with open(output_file_path, mode='w', newline='') as file:
    # Create a CSV writer object
    writer = csv.writer(file)

    # Write the [gene-expression] section to the CSV file
    writer.writerows(gene_expression_data)
    
    # Add an empty line between sections
    writer.writerow([])
    
    # Write the [libraries] section to the CSV file
    writer.writerows(libraries_data)
    
    # Add an empty line between sections
    writer.writerow([])

    # Write the [samples] section to the CSV file
    writer.writerows(samples_data)

print('CSV file created successfully at {}'.format(output_file_path))
