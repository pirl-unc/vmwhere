# VMwhere? VMhere!

VariantMotifwhere(VMwhere) performs several key analyses in tandem repeat regions of the genome from long-read sequencing data including (1) repeat genoptying to identify allele length and repeat length and (2) motif centric sequence decompsition to determine motif variants and their location witin the tandem repeat and (3) provides the ability to visualize sequence resolved alleles at the repeat region.


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
  Rscript -e "install.packages(readLines('requirements-r.txt'), repos='https://cloud.r-project.org')"
  ```


---


## Example Usage 

### 1. Find microsatellites

```bash
vmwhere find \
  --
  --

```


### 2. Genotype microsatellites

```bash
vmwhere genotype \
  --sample_id example_sample \
  --bed_file data/T2T_regions.bed \
  --bam_file data/A673_sampled_reads.sorted.bam \
  --fasta data/GCF_009914755.1_T2T-CHM13v2.0_chr6_chr10.fasta \
  --cluster_distance 4 \
  --minor_threshold 0.20 \
  --major_threshold 0.80 \
  --output_dir output/
  --num_processes 24 
```

This will:
- Extract and decompose reads overlapping the microsatellites in the bed file
- Cluster read structures by edit distance
- Call alleles based on read support with user-defined thresholds
- Output a TSV file with the genotype of each provided microsatellite


Shell script to run provided example on two regions: [`run_vmwhere_profile.sh`]


#### Input Parameter Details

**Required Parameters:**
- `--sample_id`: Unique identifier for your sample (used in output filename).
- `--bam_file`: Path to the **sorted** BAM file containing aligned reads (with indexed bam file in same directory).
- `--fasta`: Path to the reference genome FASTA.
- `--bed_file`: Headerless BED file with microsatellites to genotype. Must have at least 5 columns: `chr`, `start`, `end`, `region_id`, `motif`.
- `--output_dir`: Directory to write output files (file will be named sample_id_vmwhere_results.tsv).

**Optional Parameters:**
- `--cluster_distance`: Maximum Levenshtein distance allowed to group reads into a cluster (default =  canonical motif length).
- `--minor_threshold`: Minimum read support fraction required to call a **minor** allele (default = 0.20).
- `--major_threshold`: Minimum read support fraction required to call a **major** allele (default = 0.80).

---

#### Output TSV Details

The output has the standard columns of a VCF file, but is returned as a TSV for simpler parsing. 

`CHROM` : chromosome  
`POS` : starting index of the microsatellite  
`ID` : unique microsatellite ID provided in the input bed file  
-- `REF` : microsatellite reference sequence
-- `ALT` : variant microsatellite sequence(s) in read (. if none present)
-- `END` : ending index of the microsatellites
-- `MOTIF`: canonical motif
-- `GT` : genotype (e.g., 0/0 only reference alleles found, 0/1 one reference one alt allele found etc.)
-- `AL` : the length of the microsatellite allele in base pairs
-- `CN` : the copy number/repeat number of the canonical motif (consecutive and non-consecutive)
-- `CNM`: the maximum consecutive copy number/repeat number of the canonical motif (e.g., 6GGAA_1GGAT_2GGAA = 6)
-- `MD`: the motif density/purity of the allele (fraction of base pairs contributing to canonical motif)
-- `DS_READ` : decomposed sequence of the microsatellite alleles 
-- `DS_REF`: decomposed sequence of the reference sequence
-- `RS`: read support for the allele

### 2. Visualize alleles at microsatellites of interest

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


