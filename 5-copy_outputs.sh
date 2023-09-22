#!/bin/bash

# Check if the number of arguments is correct
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <destination_path> <gem_id>"
    exit 1
fi

# Extract the arguments
destination_dir=$1
gem_id=$2

# Create the "reports" folder in the destination path
reports_dir="${destination_dir}/${gem_id}_reports"
mkdir -p "$reports_dir"

# Define the base source path
base_source_path="/home/groups/singlecell/mabdalfttah/projects/DOLSORI_05/jobs/${gem_id}/${gem_id}/outs/per_sample_outs/"

# Get a list of sample folders
sample_folders=$(ls "$base_source_path")

# Loop over each sample folder
for sample in $sample_folders; do
    # Define the source and destination paths for the current sample
    source_path="${base_source_path}${sample}/web_summary.html"
    new_filename="${gem_id}_${sample}_web_summary.html"
    destination_path="${reports_dir}/${new_filename}"

    # Copy the file
    cp "$source_path" "$destination_path"

    echo "Copied $source_path to $destination_path"
done
