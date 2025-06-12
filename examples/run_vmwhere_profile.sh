#!/bin/bash
#SBATCH --job-name=vmwhere_test
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=4:00:00
#SBATCH --output=vmwhere_test_%j.out
#SBATCH --error=vmwhere_test_%j.err



# Example usage of vmwhere with test data
echo "Running vmwhere profile with example data..."


# Check if vmwhere is installed
if ! command -v vmwhere &> /dev/null; then
    echo "Error: vmwhere is not installed. Please run 'pip install -e .' from the project root."
    exit 1
fi


# kick off profiling function
vmwhere profile \
    --sample_id "example_sample" \
    --bam_file data/A673_sampled_reads.sorted.bam \
    --motif "GGAA" \
    --fasta data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
    --cluster_distance 4 \
    --minor_threshold 0.20 \
    --major_threshold 0.80 \
    --bed_file data/T2T_regions.bed \
    --output_dir output/

echo "Example run complete! Check output/ directory for example results."


