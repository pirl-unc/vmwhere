import argparse
import subprocess
from importlib.resources import files
from src.pipeline import run_pipeline

def profile_repeats(args):
    if args.cluster_distance is None:
        args.cluster_distance = len(args.motif)

    run_pipeline(
        sample_id=args.sample_id,
        bam_file=args.bam_file,
        motif=args.motif,
        fasta=args.fasta,
        cluster_dist=args.cluster_distance,
        minor_threshold=args.minor_threshold,
        major_threshold=args.major_threshold,
        bed_file=args.bed_file,
        output_dir=args.output_dir
    )


def visualize_region(args):
    r_script_path = files("src").joinpath("region_visualization.R")
    subprocess.run([
        "Rscript",
        str(r_script_path),
        args.sample_csv,
        args.chr,
        str(args.start),
        args.output_pdf], check=True)
    
    

def main():
    parser = argparse.ArgumentParser(
            description="vmwhere: tool for microsatellite genotyping and sequence resolved allele calling"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- Subcommand: profile ---
    profile_parser = subparsers.add_parser("profile", help="Profile reads at repeat regions given a sample BAM file")
    profile_parser.add_argument("-s", "--sample_id", required=True)
    profile_parser.add_argument("-b", "--bam_file", required=True)
    profile_parser.add_argument("-m", "--motif", required=True)
    profile_parser.add_argument("-f", "--fasta", required=True)
    profile_parser.add_argument("-cl", "--cluster_distance", type=int, default=None)
    profile_parser.add_argument("-mit", "--minor_threshold", type=float, default=0.20)
    profile_parser.add_argument("-mat", "--major_threshold", type=float, default=0.80)
    profile_parser.add_argument("-d", "--bed_file", required=True)
    profile_parser.add_argument("-o", "--output_dir", required=True)
    profile_parser.set_defaults(func=profile_repeats)

    # --- Subcommand: visualize ---
    vis_parser = subparsers.add_parser("visualize", help="Visualize results for a specific region as stacked bar plot")
    vis_parser.add_argument("--sample_csv", required=True, help="Path to vmwhere output for the sample")
    vis_parser.add_argument("--chr", required=True, help="Chromosome the region of interest is on ie chr1")
    vis_parser.add_argument("--start", required=True, type=int, help = "The start index of the region")
    vis_parser.add_argument("--output_pdf", required=True, help = "Output file name and path")
    vis_parser.set_defaults(func=visualize_region)

    # Parse args and dispatch
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()

