#!/bin/bash
#SBATCH --job-name=vmwhere_find_test
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00


# Check if vmwhere is installed
if ! command -v vmwhere &> /dev/null; then
    echo "Error: vmwhere is not installed. Please run 'pip install -e .' from the project root."
    exit 1
fi


# kick off genotyping function
vmwhere find \
	--motif GGAA \
	--fasta data/GCF_009914755.1_T2T-CHM13v2.0_chr6_chr10.fasta \
	--perfect_repeats 4 \
	--max_gap 50 \
	--buffer_size 50 \
	--output_dir output/
