import argparse
from pipeline.py import run_pipeline

def vmhere():
    parser = argparse.ArgumentParser(
        description="Genotyping and sequence-resolved allele calling for tandem repeat regions"
    )
    parser.add_argument("-s", "--sample_id", required=True, type=str, help="Sample ID")
    parser.add_argument("-b", "--bam_file", required=True, type=str, help="Path to sorted BAM file")
    parser.add_argument("-m", "--motif", required=True, type=str, help="Motif sequence to analyze")
    parser.add_argument("-f", "--fasta", type=str, help="Reference FASTA file")
    parser.add_argument("-c", "--chr_map", type=str, default="chr_mapping_simple.txt", help="Chromosome mapping file to convert between chromosome nomenclature")
    parser.add_argument("-cl", "--cluster_dist", type=int, default=4, help="Edit distance threshold for clustering reads")
    parser.add_argument("-mt", "--minor_thresh", type=float, default=0.15, help="Minor allele frequency threshold")
    parser.add_argument("-ht", "--homozygous_thresh", type=float, default=0.8, help="Homozygous threshold")
    parser.add_argument("-d", "--bed_file", required=True, type=str, help="BED file with regions of interest")
    parser.add_argument("-o", "--output_dir", required=True, type=str, help="Directory for output files")

    args = parser.parse_args()

    run_pipeline(
        sample_id=args.sample_id,
        bam_file=args.bam_file,
        motif=args.motif,
        fasta=args.fasta,
        chr_map=args.chr_map,
        cluster_dist=args.cluster_dist,
        minor_thresh=args.minor_thresh,
        homozygous_thresh=args.homozygous_thresh,
        bed_file=args.bed_file,
        output_dir=args.output_dir
    )


