# VMwhere? VMhere!

VariantMotifwhere(VMwhere) performs several key analyses in tandem repeat regions of the genome including (1) repeat genoptying to identify expansion and contraction of motifs (2) motif centric sequence decompsition to determine motif variants and their location witin the tandem repeat (3) visualization of sequence resolved alleles at the repeat region.

---

## Installation 

VMwhere is available as a python package

```bash
pip install vmwhere
``` 

---
## Requirements

- **Python** ≥ 3.11

- **R** (≥ 4.0 recommended)  
  Visualization depends on `visualize_region.R`. Install required R packages using:

  ```bash
  Rscript -e "install.packages(readLines('requirements-r.txt'))"
  ```


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
- Output two files (1) `*_clustered_results.csv`: all sequence resolved reads and their clusters and (2) `*_allele_results.csv`: alleles called based on read support and major/minor thresholds


Shell script to run provided example on two regions: [`run_vmwhere_profile.sh`]


#### Input Parameter Details

**Required Parameters:**
- `--sample_id`: Unique identifier for your sample (used in output filenames).
- `--bam_file`: Path to the **sorted** BAM file containing aligned reads (with indexed bam file in same directory).
- `--motif`: The repeat motif to genotype (e.g., `GGAA`, `GCAT`).
- `--fasta`: Path to the reference genome FASTA file (used for sequence context).
- `--bed_file`: BED file with regions to query. Must have at least 4 columns: `chr`, `start`, `end`, and `region_id`.
- `--output_dir`: Directory to write output files.

**Optional Parameters:**
- `--cluster_distance`: Maximum Levenshtein distance allowed to group reads into a cluster (default = motif length).
- `--minor_threshold`: Minimum read support fraction required to call a **minor** allele (default = 0.20).
- `--major_threshold`: Minimum read support fraction required to call a **major** allele (default = 0.80).

---


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
- Output `.pdf`: visualization per region


Shell script to run provided example on one profiled region: [`run_vmwhere_visualize.sh`]


