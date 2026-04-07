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
        output_dir=args.output_dir,
        num_processes=args.num_processes
    )


def visualize_microsatellite(args):
    r_script_path = files("vmwhere").joinpath("visualize_region.R")

    subprocess.run([
        "Rscript",
        str(r_script_path),
        args.genotype_tsv,
        args.microsatellite_id,
        args.min_allele_count,
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
    
    profile_parser.add_argument("--sample_id", required=True, help="Output file will be sample_id_vmwhere_results.tsv")
    profile_parser.add_argument("--bed_file", required=True, help="Header free bed file with columns chr start end region_id motif")
    profile_parser.add_argument("--bam_file", required=True, help="Sorted, indexed, sample bam file")
    profile_parser.add_argument("--fasta", required=True, help="Path to reference fasta file")
    profile_parser.add_argument("--cluster_distance", type=int, default=0, help="Edit distance to use when clustering reads prior to allele calling")
    profile_parser.add_argument("--minor_threshold", type=float, default=0.20, help="Minimium locus read support (fraction) to be called an allele")
    profile_parser.add_argument("--major_threshold", type=float, default=0.80, help="Read support (fraction) for calling homozygous microsatellites")
    profile_parser.add_argument("--output_dir", required=True, help="Parent directory for genotyping results")
    profile_parser.add_argument("--num_processes", type=int, default=24)
    profile_parser.set_defaults(func=genotype_microsatellites)

    # --- Subcommand: visualize ---
    vis_parser = subparsers.add_parser("visualize", help="Visualize sequence resolved alleles for a specific region")
    
    vis_parser.add_argument("-g", "--genotype_tsv", required=True, help="genotype results in vmwhere tsv format")
    vis_parser.add_argument("-m", "--microsatellite_id", required=True, help="The unique region_id in the genotype output file")
    vis_parser.add_argument("-o", "--output_pdf", required=True, help = "Output file path and name")
    vis_parser.add_argument("-c", "--min_allele_count", required=True, default=0, help="Filter out low frequency alleles")
    vis_parser.set_defaults(func=visualize_microsatellite)

    # Parse args and dispatch
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()

