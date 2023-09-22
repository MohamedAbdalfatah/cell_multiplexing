# This script initializes the filesystem of this project:
# It creates a "jobs" folder which contains as many subdirectories as samples it has
# For each sample directory, it creates the following files/folders:
# 1. fastq: dir with the symlinks pointing to the fastq files
# 2. log: dir which contains standard error and output of cellranger
# 3. (sample_id).cmd: job script to compute the features-barcode matrix using cellranger


# Import required packages
import numpy as np
import pandas as pd
import os
import argparse
import subprocess
import re
import sys
import config_vars as cfg
from utils import *


# Define command-line arguments
parser = argparse.ArgumentParser(description = "options to initialize the filesystem and scripts of this project")
parser.add_argument("--subproject",
                    dest = "subproject",
                    action = "store",
                    default = None,
                    help = "Subproject we are working on (i.e. BCLLATLAS_10)")
parser.add_argument("--gem_id",
                    dest = "gem_id",
                    action = "store",
                    default = None,
                    help = "Gel Beads in Emulsion id")
parser.add_argument("--verbose",
                    dest = "verbose",
                    action = "store_true",
                    default = False,
                    help = "Print log in standard error")
parser.add_argument("--metadata",
                    dest = "metadata",
                    action = "store",
                    default = None,
                    help = "Metadata csv file for the tonsil atlas project")
parser.add_argument("--fastq_paths",
                    dest = "fastq_paths",
                    action = "store",
                    default = None,
                    help = "File that contains the paths of the fastqs for the subproject libraries")
                    

def create_fastq_symlink_nh(gem_id, fastq_path_df, symlink_path):
    """Creates a symbolic link pointing to a fastq file using cellranger notation

    Args:
      gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
      fastq_path_df: pandas dataframe with the fastq paths for that gem_id
      symlink_path: string specifying where to create the symlinks

    Returns:
      None
    """
    pair_ids = np.unique(fastq_path_df["pair_id"])
    for i in range(len(pair_ids)):
        filt = (fastq_path_df["pair_id"] == pair_ids[i])
        pair_df = fastq_path_df.loc[filt, :]
        for j in pair_df.index:
            fastq_path = pair_df.loc[j, "fastq_path"]
            lane = str(i + 1)
            read = pair_df.loc[j, "read"]
            read = read.replace("R", "")
            subprocess.run(["ln", "-s", fastq_path, "{}/{}_S1_L00{}_R{}_001.fastq.gz".format(symlink_path, gem_id, lane, read)])


def make_cellranger_nh(gem_id, jobscript_path, fastq_path, expected_cells):
    """Creates a cellranger script for a non-hashed GEM well

    Args:
      gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
      jobscript_path: path to save the jobscript
      fastq_path: path to the fastq files
      expected_cells: expected number of high-quality cells in this experiment

    Returns:
      None
    """
    job_script_file = open("{}/{}.cmd".format(jobscript_path, gem_id), "w")
    job_script = """#!/bin/bash
#SBATCH --job-name={}

#SBATCH --mail-type=all        # send email when job begins, ends, or fails
#SBATCH --mail-user=mohamed.abdalfttah@cnag.crg.eu

#SBATCH --output=%x.slurm.%J.out        # define where our output and error from the job will be stored
#SBATCH --error=%x.slurm.%J.err

#SBATCH --time=11:00:00 # set a maximum time that the job will take HH:MM:SS (process will be terminated after this is reached)

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=normal
#SBATCH --partition=genD
#SBATCH --mem=32G

echo [`date "+%Y-%m-%d %T"`] started job on $HOSTNAME

export TENX_IGNORE_DEPRECATED_OS=1
export HDF5_USE_FILE_LOCKING=FALSE


{} multi --id {} --csv config.csv  --jobmode /scratch/groups/singlecell/software/cellranger/6.1.1/external/martian/jobmanagers/new_cluster/slurm.template;
""".format(gem_id, cfg.cellranger_path, gem_id)
    job_script_file.write(job_script)
    job_script_file.close()


options = parser.parse_args()
subproject = options.subproject
gem_id = options.gem_id
metadata_path = options.metadata
fastq_paths = options.fastq_paths


# Read data
project_dir = "/home/groups/singlecell/mabdalfttah/projects/{}".format(subproject)
fastq_path_df = pd.read_csv(fastq_paths, sep = "\t", header = 0)
metadata_df = pd.read_csv(metadata_path, sep = ",", header = 0)
if options.verbose:
    sys.stderr.write("Files read successfully!\n")


# For each sample, create directories and jobscript
if not os.path.exists("{}/jobs".format(project_dir)):
    os.mkdir("{}/jobs".format(project_dir))
filt = (metadata_df["gem_id"] == gem_id)
metadata_df = metadata_df.loc[filt]


# Create directories
subproject_dir = "{}/jobs/{}".format(project_dir, gem_id)
fastq_dir = "{}/fastq".format(subproject_dir)
log_dir = "{}/log".format(subproject_dir)
for direct in [subproject_dir, fastq_dir, log_dir]:
    if not os.path.exists(direct):
        os.mkdir(direct)


# Define variables and subset dataframes
library_id = metadata_df.loc[filt, "library_id"]
fastq_sub_df = fastq_path_df.loc[fastq_path_df["library_id"].isin(library_id), :]
type = metadata_df["type"]
type = type.values[0]

# Create symmlinks to fastq files
create_fastq_symlink_nh(gem_id, fastq_sub_df, fastq_dir)
# Create cellranger script
fastq_dir_abs = "{}{}".format(cfg.project_path, fastq_dir)
make_cellranger_nh(gem_id, subproject_dir, fastq_dir_abs, 5000)


