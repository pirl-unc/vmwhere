import argparse
from .pipeline import run_pipeline

def profile_repeats(args):
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


def visualize_region(args):
    # Placeholder function — replace with real logic
    print(f"Visualizing region {args.region_id} in sample {args.sample_id}")
    # TODO: call your actual visualization function

def main():
    parser = argparse.ArgumentParser(
            description="vmwhere: tool for microsatellite genotyping and sequence resolved allele calling"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- Subcommand: profile ---
    profile_parser = subparsers.add_parser("profile", help="Profile motif regions from BAM + BED")
    profile_parser.add_argument("-s", "--sample_id", required=True)
    profile_parser.add_argument("-b", "--bam_file", required=True)
    profile_parser.add_argument("-m", "--motif", required=True)
    profile_parser.add_argument("-f", "--fasta", required=True)
    profile_parser.add_argument("-cl", "--cluster_distance", type=int, default=None)
    profile_parser.add_argument("-mit", "--minor_threshold", type=float, default=0.20)
    profile_parser.add_argument("-mat", "--major_threshold", type=float, default=0.8)
    profile_parser.add_argument("-d", "--bed_file", required=True)
    profile_parser.add_argument("-o", "--output_dir", required=True)
    profile_parser.set_defaults(func=profile_repeats)

    # --- Subcommand: visualize ---
    vis_parser = subparsers.add_parser("visualize", help="Visualize results for a specific region")
    vis_parser.add_argument("-s", "--sample_id", required=True, help="Name of Sample to visualize region from")
    vis_parser.add_argument("-r", "--region_id", required=True, help="Region ID to visualize")
    vis_parser.set_defaults(func=visualize_region)

    # Parse args and dispatch
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()

