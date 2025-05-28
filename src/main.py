import argparse
from pipeline import run_pipeline

def profile_repeats():
    parser = argparse.ArgumentParser(
        description="Genotyping and sequence-resolved allele calling for tandem repeat regions"
    )
    parser.add_argument("-s", "--sample_id", required=True, type=str, help="Sample ID")
    parser.add_argument("-b", "--bam_file", required=True, type=str, help="Path to sorted BAM file")
    parser.add_argument("-m", "--motif", required=True, type=str, help="Motif sequence to analyze")
    parser.add_argument("-f", "--fasta", type=str, help="Reference FASTA file")
    parser.add_argument("-cl", "--cluster_dist", type=int, default=None, help="Edit distance threshold for clustering reads (default: length of motif)" )
    parser.add_argument("-mt", "--minor_thresh", type=float, default=0.15, help="Minor allele frequency threshold")
    parser.add_argument("-ht", "--homozygous_thresh", type=float, default=0.8, help="Homozygous threshold")
    parser.add_argument("-d", "--bed_file", required=True, type=str, help="BED file with regions of interest")
    parser.add_argument("-o", "--output_dir", required=True, type=str, help="Directory for output files")

    args = parser.parse_args()

    # Set dynamic default for cluster_dist based on motif length
    if args.cluster_dist is None:
        args.cluster_dist = len(args.motif)

    run_pipeline(
        sample_id=args.sample_id,
        bam_file=args.bam_file,
        motif=args.motif,
        fasta=args.fasta,
        cluster_dist=args.cluster_dist,
        minor_thresh=args.minor_thresh,
        homozygous_thresh=args.homozygous_thresh,
        bed_file=args.bed_file,
        output_dir=args.output_dir
    )


