"""
The purpose of this python3 script is to implement the main API functions of vmwhere.
"""

import pysam
import argparse
import pandas as pd
import glob
import os
from Bio.Seq import Seq
from Bio import SeqIO
import multiprocessing as mp
from multiprocessing import Pool
from datetime import date
import Levenshtein
import re


parser = argparse.ArgumentParser(description="Process BAM files and extract reads overlapping repeat motif regions.")
parser.add_argument("-s", "--sample_id", type=str, required=True, help="Unique sample identifier")
parser.add_argument("-b", "--bam_file", type=str, required=True, help="Path to the BAM file")
parser.add_argument("-m", "--motif", type=str, required=True, help="Sequence of interest")
parser.add_argument("-f", "--fasta", type=str, default="/proj/dllab/human_refs/ncbi_dataset/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna", help="Reference FASTA file (default is T2T CHM12v2.0")
parser.add_argument("-c", "--chr_map", type=str, default="/proj/dllab/human_refs/ncbi_dataset/data/GCF_009914755.1/chr_mapping_simple.txt", help="Chromosome name mapping file")
parser.add_argument("-cl", "--cluster_dist", type=int, default=4, help="Distance for clustering reads")
parser.add_argument("-mt", "--minor_thresh", type=float, default=0.15, help="Minor allele threshold [default: 0.15]")
parser.add_argument("-ht", "--homozygous_thresh", type=float, default=0.8, help="Homozygous threshold [default: 0.8]")
parser.add_argument("-d", "--bed_file", type=str, required=True, help="Headerless tab separated bed file")
parser.add_argument("-o", "--output_dir", type=str, required=True, help="Output directory to save results")


def run_pipeline():


def run_process_reads():


def run_cluster_reads():


def run_call_alleles():



def run_visualize_microsatellite():


