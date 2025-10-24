import argparse
import subprocess
from importlib.resources import files
from .find import run_find
from .genotyper import run_genotyper


def find_microsatellites(args):

    run_find(
            motif=args.motif,
            fasta_file=args.fasta,
            repeats=args.perfect_repeats,
            gap=args.max_gap,
            output_dir=args.output_dir,
            buffer=args.buffer_size
            )


def genotype_microsatellites(args):

    run_genotyper(
        sample_id=args.sample_id,
        bam_file=args.bam_file,
        fasta=args.fasta,
        cluster_distance=args.cluster_distance,
        minor_threshold=args.minor_threshold,
        major_threshold=args.major_threshold,
        bed_file=args.bed_file,
        output_dir=args.output_dir
    )


def visualize_microsatellite(args):
    r_script_path = files("vmwhere").joinpath("visualize_region.R")

    subprocess.run([
        "Rscript",
        str(r_script_path),
        args.sample_csv,
        args.chr,
        str(args.start),
        args.output_pdf], check=True)
    
    

def main():
    parser = argparse.ArgumentParser(
            description="vmwhere: microsatellite reference identification, sample genotyping, sequence decomposition, and visualization from long-read data"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    
    # --- Subcommand: find ---
    find_parser = subparsers.add_parser("find", help="Identify genomic coordinates of repeat microsatellite sequences based on a reference")
    
    find_parser.add_argument("-m", "--motif", required=True)
    find_parser.add_argument("-r", "--perfect_repeats", type=int, default=2)
    find_parser.add_argument("-g", "--max_gap", type=int, default=50)
    find_parser.add_argument("-b", "--buffer_size", type=int, default=50)
    find_parser.add_argument("-o", "--output_dir", required=True)
    find_parser.add_argument("-f", "--fasta", required=True)
    find_parser.set_defaults(func=find_microsatellites)


    # --- Subcommand: genotype ---
    profile_parser = subparsers.add_parser("genotype", help="Genotype microsatellites given a sample BAM file")
    
    profile_parser.add_argument("-s", "--sample_id", required=True)
    profile_parser.add_argument("-b", "--bam_file", required=True)
    profile_parser.add_argument("-f", "--fasta", required=True)
    profile_parser.add_argument("-cl", "--cluster_distance", type=int, default=0)
    profile_parser.add_argument("-mit", "--minor_threshold", type=float, default=0.20)
    profile_parser.add_argument("-mat", "--major_threshold", type=float, default=0.80)
    profile_parser.add_argument("-d", "--bed_file", required=True)
    profile_parser.add_argument("-o", "--output_dir", required=True)
    profile_parser.set_defaults(func=genotype_microsatellites)

    # --- Subcommand: visualize ---
    vis_parser = subparsers.add_parser("visualize", help="Visualize sequence resolved alleles for a specific region")
    
    vis_parser.add_argument("--sample_tsv", required=True, help="Name for the sample")
    vis_parser.add_argument("--chr", required=True, help="Chromosome the region of interest is on ie chr1")
    vis_parser.add_argument("--start", required=True, type=int, help = "The start index of the region")
    vis_parser.add_argument("--output_pdf", required=True, help = "Output file name and path")
    vis_parser.set_defaults(func=visualize_microsatellite)

    # Parse args and dispatch
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()

