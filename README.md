# VMwhere? VMhere!

VariantMotifwhere(VMwhere) performs several key analyses in tandem repeat regions of the genome including (1) repeat genoptying to identify expansion and contraction of motifs (2) motif centric sequence decompsition to determine motif variants and their location witin the tandem repeat (3) visualization of sequence resolved alleles at the repeat region.


## Installation 

VMwhere is available as a python package

```
pip install vmwhere
``` 

## Input Requirements

- **BAM file**: Aligned sequencing reads (`--bam_file`)
- **BED file**: Tab-delimited file of motif regions (no header) (`--bed_file`)
- **Motif**: Nucleotide sequence of the motif of interest (e.g., `GGAA`) (`--motif`)
- **Reference genome**: FASTA file (`--fasta`)
- **Sample ID**: Unique identifier for the sample (`--sample_id`)

---

## Example Usage 


### 1. Profile repeat regions

```bash
vmwhere profile \
  --sample_id example_sample \
  --bam_file data/A673_sampled_reads.sorted.bam \
  --motif GGAA \
  --fasta data/GCF_009914755.1_T2T-CHM13v2.0_chr6_chr10.fasta \
  --cluster_distance 4 \
  --minor_threshold 0.20 \
  --major_threshold 0.80 \
  --bed_file data/T2T_regions.bed \
  --output_dir output/
```

This will:
- Extract and decompose reads overlapping motif regions
- Cluster read structures by edit distance
- Call alleles with user-defined thresholds
- Output `*_clustered_results.csv` and `*_allele_results.csv`


Shell script to run provided example on two regions: [`run_vmwhere_profile.sh`]


### 2. Visualize a specific region

```bash
vmwhere visualize \
  --sample_csv output/example_sample_allele_results.csv \
  --chr chr6 \
  --start 6706603 \
  --output_pdf output/chr6_example_region_visualization.pdf
```

This will:
- Generate a stacked bar plot showing distinct read structures (alleles) and their read count


Shell script to run provided example on one profiled region: [`run_vmwhere_visualize.sh`]

---

## Input Requirements

- **BAM file**: Aligned and indexed sequencing reads (`--bam_file`)
- **BED file**: Tab-delimited file of motif regions (no header) (`--bed_file`)
- **Motif**: Nucleotide sequence of the motif of interest (e.g., `GGAA`) (`--motif`)
- **Reference genome**: FASTA file (`--fasta`)
---

##  Outputs

- `*_clustered_results.csv`: All sequence-resolved reads and clusters
- `*_allele_results.csv`: Filtered alleles per region with support
- `.pdf`: Optional per-region visualizations

