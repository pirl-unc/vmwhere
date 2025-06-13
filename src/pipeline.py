import pysam
import argparse
import pandas as pd
from importlib.resources import files
import glob 
import os
from Bio.Seq import Seq
from Bio import SeqIO
import multiprocessing as mp
from multiprocessing import Pool
from datetime import date
import Levenshtein
import re


## define functions
def fetch_reference_sequence(fasta_file):
    """
    Fetch the reference chromosome sequence

    Parameters:
      fasta_file : FASTA file
      chrom_names_mapping_file : maps chromosome names in ref seq style to ucsc style

    Returns:
      Dict[chromosome name, chromosome sequence]

    """

    reference_sequence = {} # dictionary to store the chromosome name and the reference sequence

    for seq_record in SeqIO.parse(fasta_file, 'fasta'): # each seq_record is a chromosome, loop over each seq_record
        chr_name_key = seq_record.id # get chromosome name
        reference_sequence[chr_name_key] = seq_record.seq # generates the dict for that chr using the ucsc chr name where seq_record.seq pulls the sequence
    return reference_sequence

def count_max_consecutive_motifs(sequence, motif):
    """
    Counts the number of consecutive occurrences of a motif in a given sequence.

    Parameters:
        sequence (str): The sequence to search within.
        motif (str): The motif to look for.

    Returns:
        int: The number of consecutive motif occurrences, or 0 if the motif is not found.
    """

    max_sequence_repeats = 0
    start_index = 0

    while start_index < len(sequence):
        # find starting index of motif
        start_index = sequence.find(motif, start_index)
        if start_index == -1:
            break

        # count consecutive occurrences from this position
        motif_count_consecutive = 1
        while sequence[start_index + len(motif)*motif_count_consecutive : start_index + len(motif) * (motif_count_consecutive + 1)] == motif:
            motif_count_consecutive += 1

        # update maximum repeats
        max_sequence_repeats = max(max_sequence_repeats, motif_count_consecutive) # compare current motif repeat to previously found motif repeat, retain larger value
        start_index += len(motif)*motif_count_consecutive  # move to end of motif to keep searching sequenceuntil you reach the end of the sequence

    return max_sequence_repeats

def count_total_motifs(sub_sequence, motif):
    """
    Counts the total number of occurrences of a motif in a given sequence (consecutive and non-consecutive).

    Parameters:
        sequence (str): The sequence to search within.
        motif (str): The motif to look fo in the string.

    Returns:
        int: The total number of motif occurrences in the sequence.
    """
    total_repeat_num = sub_sequence.upper().count(motif)
    return total_repeat_num

def hamming_distance(s1, s2):
    """Calculate the Hamming distance between two sequences of equal length."""
    if len(s1) != len(s2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def is_single_nucleotide_variant(seq_chunk, motif_seq):
    """Check if a sequence is a single nucleotide variant of the motif (ie an snv). Return True if it is."""
    return hamming_distance(seq_chunk, motif_seq) == 1

def is_homopolymer(seq):
    """check if a sequence is a homopolymer (all bases are the same). Return True if it is."""
    return len(seq) > 0 and all(bp == seq[0] for bp in seq)


def recompose_string_from_structure(decomposed_read_structure_to_check):

    recomposed_read_seq = []
    for decomposed_string in decomposed_read_structure_to_check.split('_'):
        regex_decomp = re.findall(r"(\d+)([a-zA-Z]+)", decomposed_string)
        recomposed_sub_string_seq = int(regex_decomp[0][0])*regex_decomp[0][1]
        recomposed_read_seq.append(recomposed_sub_string_seq)

    return ''.join(recomposed_read_seq)

def decompose_seq_with_motif_anchors(trim_seq, motif):
    """
    Determines the structure of a read sequence by first finding the exact matches to the motif, and then decomposing the sequences in between exact motifs by minimizing the edit distance to a perfect repeat structure.

    Parameters:
        read_query_sequence: str, sequence of the microsatellite(starting at first motif and ending at the lastoccurance) that overlaps the region
        motif: str, motif of interest (e.g., "GGAA")

    Returns:
        str: Structured summary in format like "3GGAA_1TC_2GAAA_3GGAA"
    """
    # check that there are any motifs present in the read
    if trim_seq.find(motif) == -1:
        return None, 0, '' # no motif so we return None as the structure, 0 for the total length, and an empty string as the sequence

    motif_len = len(motif)
    n = len(trim_seq)

    # first, locate all exact motif positions
    motif_start_positions = []
    i = 0
    while i <= n - motif_len:
        if trim_seq[i:i + motif_len] == motif:
            motif_start_positions.append(i)
            i += motif_len
        else:
            i += 1

    read_structure = []
    idx = 0
    while idx < len(motif_start_positions):
        start = motif_start_positions[idx]
        count = 1
        # second, pull together consecutive motifs based on their start positions
        while (idx + count < len(motif_start_positions) and motif_start_positions[idx + count] == start + count * motif_len):
            count += 1
        read_structure.append(f"{count}{motif}")

        # third, process the gap sequence between motif repeats
        # Intervening sequence between current motif block and next
        next_idx = idx + count
        if next_idx < len(motif_start_positions):
            gap_start = start + count * motif_len
            gap_end = motif_start_positions[next_idx]
            inter_seq = trim_seq[gap_start:gap_end]
        # decompose the intervening sequence using edit distance minimzation
            if inter_seq:
                if len(inter_seq) <= 4:
                    read_structure.append(f"1{inter_seq}")
                elif is_homopolymer(inter_seq):
                    read_structure.append(f"1{inter_seq}")
                else:
                    parsed = further_parse_intervening(inter_seq, motif)
                    read_structure.append(parsed)

        idx += count

    return "_".join(read_structure), n, trim_seq

def further_parse_intervening(intervening_seq, motif):
    """
    Further parse intervening sequences that are greater than 4 bps by finding segments that minimize
    edit distance to the motif.

    Parameters:
        intervening_seq: str, the intervening sequence to parse (must be greater than 4 bps)
        motif: str, the motif to look for

    Returns:
        str: Structured representation of the intervening sequence
    """

    n = len(intervening_seq)
    motif_len = len(motif)

    # Initialize dynamic programming table for minimizing edit distance
    # dp[i] = minimum edit distance up to position i
    dp = [float('inf')] * (n + 1)
    dp[0] = 0

    # Store backtracking information
    prev_pos = [0] * (n + 1)
    segment_len = [0] * (n + 1)  # Length of the segment ending at position i

    # Fill the dp table
    for i in range(1, n + 1):
        # Try segments of different lengths
        for seg_len in range(1, min(i, motif_len) + 1):
            # Get current segment
            segment = intervening_seq[i-seg_len:i]

            # Calculate edit distance to the motif if segment length equals motif length
            # Otherwise, just count each character as a 1-cost gap
            if seg_len == motif_len:
                edit_dist = hamming_distance(segment, motif)
            else:
                edit_dist = seg_len  # Each character costs 1

            # Update if this gives a better solution
            if dp[i-seg_len] + edit_dist < dp[i]:
                dp[i] = dp[i-seg_len] + edit_dist
                prev_pos[i] = i-seg_len
                segment_len[i] = seg_len

# Reconstruct the optimal decomposition by backtracking
    structure = []
    current_pos = n
    last_segment = ""
    count = 0

    while current_pos > 0:
        seg_len = segment_len[current_pos]
        segment = intervening_seq[current_pos - seg_len:current_pos]

        if seg_len == motif_len:
            if segment == last_segment:
                count += 1
            else:
                if last_segment:
                    structure.insert(0, f"{count}{last_segment}")
                last_segment = segment
                count = 1
        else:
            if last_segment:
                structure.insert(0, f"{count}{last_segment}")
                last_segment = ""
                count = 0
            structure.insert(0, f"1{segment}")

        current_pos = prev_pos[current_pos]

    if last_segment:
        structure.insert(0, f"{count}{last_segment}")

    return "_".join(structure)

def classify_variants_from_structure(structure, motif):
    """
    Given a read structure string (e.g. "3GGAA_1TC_2GAAA_3GGAA"),
    classify and count:
      - SNV motif variants: segments that are the same length as the motif but differ by 1 bp.
      - multi nucleotide variants: segments that are more than 1 bp different from the motif but the same length of the motif.
      - Single nucleotide in-between events: exactly 1 bp between two motifs/motif variants.
      - Multi nucleotide in-between events: more than 1 bp < len(motif) bp between two motifs/motif variants.

    Returns:
      tuple: (snp_variant_count, multi_variant_count, single_bp_inbetween_count, multi_bp_inbetween_count)
    """
    perfect_motif_count_decomp = 0
    consecutive_motif_repeats_decomp = 0
    count_max_consecutive_motif_repeats_decomp = 0
    snp_variant_count = 0
    multi_bp_variant_count = 0
    single_bp_inbetween_count = 0
    multi_bp_inbetween_count = 0

    if structure is None:
        return 0, 0, 0, 0, 0, 0

    # Split the structure string on underscores.
    seq_split = structure.split('_')
    for sub_string in seq_split:
        regex = re.findall(r"(\d+)([a-zA-Z]+)", sub_string)
        count_prefix = int(regex[0][0])
        seq_string = regex[0][1]

        # characerize blocks that are the same length as the motif
        if len(seq_string) == len(motif):
            if hamming_distance(seq_string, motif) == 0:
            # if difference is 0, count as a perfect match
                perfect_motif_count_decomp += count_prefix
                consecutive_motif_repeats_decomp = count_prefix
                count_max_consecutive_motif_repeats_decomp = max(count_max_consecutive_motif_repeats_decomp, consecutive_motif_repeats_decomp)
            # If the difference is 1, count it as a SNP variant.
            elif hamming_distance(seq_string, motif) == 1:
                snp_variant_count += count_prefix
            else:
                # If the difference is more than 1, count it as a structural variant.
                multi_bp_variant_count += count_prefix
        else:
            # For sub strings that are not the same length as the motif:
            if len(seq_string) == 1:
                single_bp_inbetween_count += count_prefix
            else:
                multi_bp_inbetween_count += count_prefix

    return perfect_motif_count_decomp, count_max_consecutive_motif_repeats_decomp, snp_variant_count, multi_bp_variant_count,single_bp_inbetween_count, multi_bp_inbetween_count

def filter_reads_with_low_support(read_structure_occurance_list, region_read_results_df):
    """
    Keeps only entries where the parsed sequence string has at least 2 supporting reads per that region.

    Parameters:
    motif_occurance_list (list): List of sequence structure occurrences for that region.
    region_read_results (df): dataframe of resuts for that region, where each read is a row

    Returns:
    filtered_region_read_results (df): Filtered dataframe for that region, where each read is a row and low support parsed stuctures have been removed.
    """
    # Find sequences with at least 2 supporting reads
    supported_read_structures = {occurance for occurance in read_structure_occurance_list if read_structure_occurance_list.count(occurance) >= 2}

    # Filter the region read results to keep only entries with valid motif read numbers
    filtered_region_read_results = region_read_results_df[region_read_results_df['read_structure'].isin(supported_read_structures)]

    return filtered_region_read_results

def consolidate_reads_per_region(filtered_region_results):
    """
    Consolidates read results for each region by grouping based on unique read structures.

    Parameters:
        filtered_region_results (pd.DataFrame): DataFrame with filtered read results

    Returns:
        pd.DataFrame: Consolidated DataFrame where each row represents a unique read structure in the region and the read support is counted.
    """
    # Count occurrences of each unique read structure within each region
    read_counts = filtered_region_results.groupby(['read_structure', 'region_id']).size().reset_index(name='read_support')

    # Get unique combinations of read_structure and region_id with their attributes
    consolidated_results = filtered_region_results.drop_duplicates(subset=['read_structure', 'region_id'])

    # Merge the count information with the unique read structures
    consolidated_results = consolidated_results.merge(read_counts, on=['read_structure', 'region_id'], how='left')

    return consolidated_results


def cluster_consolidated_reads_by_edit_distance(consolidated_results, cluster_distance):
    """
    Clusters read_sequences based on Levenshtein distance and read support.

    Parameters:
        df (pd.DataFrame): Must contain 'read_sequence' and 'read_support'.
        cluster_distance (int): Max allowed distance to join a cluster.

    Returns:
        pd.DataFrame: DataFrame with cluster assignment and updated supporting reads for the cluster
    """
    results_sorted = consolidated_results.sort_values(by='read_support', ascending=False).reset_index(drop=True)
    clusters = []  # Each cluster is {'centroid': str, 'members': [int], 'ids': [int]}
    cluster_assignments = [-1] * len(results_sorted)

    for i, row in results_sorted.iterrows():
        sequence = row['read_sequence']
        best_dist = float('inf')
        best_cluster_idx = -1

        # Compare to each existing cluster centroid
        for idx, cluster in enumerate(clusters):
            dist = Levenshtein.distance(sequence, cluster['centroid'])
            if dist < best_dist and dist <= cluster_distance:
                best_dist = dist
                best_cluster_idx = idx

        if best_cluster_idx == -1:
            # Create new cluster with this sequence as centroid
            clusters.append({'centroid': sequence, 'members': [i]})
            cluster_assignments[i] = len(clusters) - 1
        else:
            # Assign to closest cluster
            clusters[best_cluster_idx]['members'].append(i)
            cluster_assignments[i] = best_cluster_idx

    clustered = results_sorted.copy()
    clustered['cluster_id'] = cluster_assignments
    # identify clustered total read support for each cluster
    cluster_read_support = clustered.groupby('cluster_id')['read_support'].sum()
    # identify the index of the centroid for each cluster
    centroids_idx = clustered.groupby('cluster_id')['read_support'].idxmax()
    # create dataframe with only centroids and drop read support column since it's wrong
    centroids_df = clustered.loc[centroids_idx]
    centroids_df = centroids_df.drop(columns=['read_support'])
    # update total read support for each cluster
    centroids_df_final = centroids_df.merge(cluster_read_support, on='cluster_id', how='left')

    return centroids_df_final

def process_reads_overlapping_regions(motif, chr_map, regions, bamfile, reference_sequence, cluster_distance):
    """
    Processes reads for regions specified in a BED-like dataframe.

    Parameters:
      regions : DataFrame with tab delineated columns [chromosome, start, end].
      bamfile : Sorted BAM file for fetching reads.
      reference_sequence (dict): Dictionary with chromosome names as keys and their sequences as values.

    Returns:
      pd.DataFrame: DataFrame containing sequence resolved structure of reads overlapping each provided region.
    """
    # determine which chromosome bedfile we are processing
    chrom = regions.iloc[0,0] # column zero contains the chromosome un ucsc naming

    ## determine ucsc or refseq naming of chromosomes in bam file
    if bamfile.references[0] in chr_map.keys():
        chrom_bam = chrom #chrom needs to stay in ucsc naming
    else:
        chrom_bam = chr_map.get(chrom) # chrom  needs to get swtiched to ref_seq naming

    ## Determine ucsc or freseq naming of chromosomes in fasta reference dictionary
    if chrom in reference_sequence.keys():
        chrom_reference = chrom
    else:
        chrom_reference = chr_map[chrom] # convert to proper nomenclature

    # Initialize an empty dictionary to store results for each chromosome
    chromosome_query_read_results = {
        'region_id' : [],
        'chr': [],
        'start': [],
        'end': [],
        'read_structure': [],
        'reference_structure':[],
        'read_support': [],
        'cluster_id': [],
        'max_consecutive_repeats': [],
        'total_motif_repeats': [],
        'snv_motifs': [],
        'mnv_motifs': [],
        'single_bp_between': [],
        'multi_bp_inbetween' : [],
        'total_length':[],
        'motif_density' : [],
        'read_sequence': []
        }

    discarded_reads=0
    reverse_motif = str(Seq(motif).reverse_complement())

    # analyse one region at a time
    for i in range(len(regions)):
        if chrom == "chrM": #skip chrM
            continue

        region_id = regions.iloc[i,3] # column 3 contains the unique region identifier
        start_coord = int(regions.iloc[i, 1])  # Column 1 contains the start position
        end_coord = int(regions.iloc[i, 2])  # Column 2 contains the end position

        # initialize an empty dictionary to store results for each region

        # verify that there are enough reads overlapping the region
        read_count = 0
        for read in bamfile.fetch(chrom_bam, start_coord, end_coord):
            if read.is_secondary:
                continue
            read_count += 1
            if read_count >= 10: # don't care how many reads, just care more than 10
                break

        # skip the region if not enough reads
        if read_count < 10:
            print(f"skipping region {chrom_bam}:{start_coord}-{end_coord} ({region_id}) due to low read count of {read_count}")
            continue

        # initialize dictionary to store results for the region
        region_query_read_results = {
                'chr': [],
                'start': [],
                'end': [],
                'read_structure': [],
                'max_consecutive_repeats': [],
                'total_motif_repeats': [],
                'total_length':[],
                'motif_density' : [],
                'snv_motifs': [],
                'single_bp_between': [],
                'mnv_motifs': [],
                'multi_bp_inbetween' : [],
                'reference_structure':[],
                'read_sequence': [],
                'region_id' : []
                }

        # Fetch reads from BAM file overlapping current region
        for read in bamfile.fetch(chrom_bam, start_coord, end_coord):
            # Only count reads where this is the primary alignment (skip secondary alignments)
            if read.is_secondary:
                continue

            qstart = None
            qend = None

            aligned_pairs = read.get_aligned_pairs()

            # Process the aligned pairs to determine overlap
            for query_pos, ref_pos in aligned_pairs:
                if ref_pos == start_coord:
                    ref_start = ref_pos
                    qstart = query_pos  # First location of overlap of the query (read) with the reference
                if ref_pos == end_coord:
                    qend = query_pos  # Last location of overlap of the query (read) with the reference
                    ref_end = ref_pos
                    break  # Exit the loop once end is found

            if qstart is not None and qend is not None:
                # trim reference to be only the part that overlaps with the region
                ref_seq = str(reference_sequence[chrom_reference][ref_start:ref_end]).upper()
                # trim read to be only part that overlaps with the region
                query_seq = read.query_sequence[qstart:qend].upper()


                # check if we are looking for the motif or its reverse complement
                ref_forward_total = count_total_motifs(ref_seq, motif)
                ref_reverse_total = count_total_motifs(ref_seq, reverse_motif)


                # if reverse motif is dominant, convert ref string and read string to reverse complement so we can find the motif
                if ref_forward_total < ref_reverse_total:
                    ref_seq = str(Seq(ref_seq).reverse_complement())
                    query_seq = str(Seq(query_seq).reverse_complement())

                ## identify start and end indices of first and last motif to determine microsatellite coordinates/sequence
                ref_seq_motif_start = ref_seq.find(motif)
                ref_seq_motif_end = ref_seq.rfind(motif)
                ref_seq_microsat = ref_seq[ref_seq_motif_start:ref_seq_motif_end + len(motif)]
                ## pull reference sequence
                ref_structure, total_length_ref, ref_microsat_seq = decompose_seq_with_motif_anchors(ref_seq_microsat, motif)

                query_seq_motif_start = query_seq.find(motif)
                query_seq_motif_end = query_seq.rfind(motif)
                query_seq_microsat = query_seq[query_seq_motif_start:query_seq_motif_end + len(motif)]

                 # perform decomposition of the read sequence and identify variants present
                read_structure, total_length, read_microsat_seq = decompose_seq_with_motif_anchors(query_seq_microsat, motif)
                perfect_motif_count_decomp, count_max_consecutive_motif_repeats_decomp, snp_variant_count, multi_bp_variant_count,single_bp_inbetween_count, multi_bp_inbetween_count = classify_variants_from_structure(read_structure, motif)

                # determine total number of prefect motifs present and max consecutive by counting them in the original string
                max_consecutive_repeats = count_max_consecutive_motifs(query_seq_microsat, motif)
                total_repeats = count_total_motifs(query_seq_microsat, motif)

                # compare motifs found in decomposed sequence to those found in the original sequence to ensure the decompositions is correct
                assert perfect_motif_count_decomp == total_repeats, f"Decomposed motif count {perfect_motif_count_decomp} does not match original sequence motif count of {total_repeats} in {read_structure} in {region_id} for read sequence {query_seq_microsat}"
                assert count_max_consecutive_motif_repeats_decomp == max_consecutive_repeats, f"Decomposed motif count {count_max_consecutive_motif_repeats_decomp} does not match original sequence motif count of {max_consecutive_repeats} in {read_structure} in {region_id} for read sequence {query_seq_microsat}"

                # recomposed the decomposed sequence and check to ensure it matches the original sequence that was decomposed
                if read_structure is not None:
                    reconstructed_seq_check = recompose_string_from_structure(read_structure)
                    assert reconstructed_seq_check == query_seq_microsat, f"Recomposed sequence {reconstructed_seq_check} from decomposed sequence {read_structure} does not match original sequence {query_seq_microsat} for {region_id}"
                    motif_density = (total_repeats*len(motif))/total_length
                else:
                    motif_density = 0

                # record read characteristics
                region_query_read_results['chr'].append(chrom)
                region_query_read_results['start'].append(start_coord)
                region_query_read_results['end'].append(end_coord)
                region_query_read_results['read_structure'].append(read_structure if read_structure is not None else "<NONE>")
                region_query_read_results['max_consecutive_repeats'].append(max_consecutive_repeats)
                region_query_read_results['total_motif_repeats'].append(total_repeats)
                region_query_read_results['total_length'].append(total_length)
                region_query_read_results['motif_density'].append(motif_density)
                region_query_read_results['snv_motifs'].append(snp_variant_count)
                region_query_read_results['single_bp_between'].append(single_bp_inbetween_count)
                region_query_read_results['mnv_motifs'].append(multi_bp_variant_count)
                region_query_read_results['multi_bp_inbetween'].append(multi_bp_inbetween_count)
                region_query_read_results['reference_structure'].append(ref_structure)
                region_query_read_results['read_sequence'].append(read_microsat_seq)
                region_query_read_results['region_id'].append(region_id)
            else:
                discarded_reads += 1

        # collect all read_structures
        motif_occurance_list = region_query_read_results['read_structure']
       # turn into dataframe for downstream filtering
        region_query_read_results_df = pd.DataFrame(region_query_read_results)

        # filter out singleton reads
        filtered_region_read_results = filter_reads_with_low_support(motif_occurance_list, region_query_read_results_df)
        filtered_count = len(filtered_region_read_results) # check that there are reads after verifying each read structure has at least 2 supporting reads

        if filtered_count == 0:
            print(f"skipping region {chrom_bam}:{start_coord}-{end_coord} ({region_id}) due to no concordance in reads")
            continue

        # consolidate read results for each region where each row is a unique read structure and the read support is counted
        consolidated_region_read_results = consolidate_reads_per_region(filtered_region_read_results)

        # cluster reads that are within the specific edit distance to each other
        clustered_region_read_results = cluster_consolidated_reads_by_edit_distance(consolidated_region_read_results, cluster_distance)

        # add filtered region results to the complete chromosome dictionary results
        for column in clustered_region_read_results.columns:
            chromosome_query_read_results[column].extend(clustered_region_read_results[column].tolist())

    # provide a summary output for each chromosome in the log file
    print(f"Processed {len(regions)} regions on {chrom}")

    return pd.DataFrame(chromosome_query_read_results)

def process_single_chromosome(args):
    """
    Process a single chromosome file at a time, with outputs coming from the function that processes each region. This function will be called by the multiprocessing pool so each chromosome is run in parallel.

    Parameters:
        args (tuple): Contains (chromosome that is processed, regions on that chromosome, reference_sequence)

    Returns:
        str: Name of the output file that was created for that chromosome
    """
    chromosome_curr, regions_df, MOTIF, BAM_FILE, chr_map, reference_sequence, cluster_distance = args

    # Open BAM file
    try:
        bam = pysam.AlignmentFile(BAM_FILE, "rb")
    except FileNotFoundError:
        print(f"Error: BAM File {BAM_FILE} not found")

    # Process the regions on the chromosome
    query_read_results_df = process_reads_overlapping_regions(
        motif=MOTIF,
        chr_map=chr_map,
        regions=regions_df,
        bamfile=bam,
        reference_sequence=reference_sequence,
        cluster_distance = cluster_distance
    )

    # Close the BAM file
    bam.close()

    # return the chromosome df
    return chromosome_curr, query_read_results_df


def identify_alleles(read_summary_df, minor_threshold, major_threshold):
    """Determine alleles based on read support."""

    read_summary_df['total_region_reads'] = read_summary_df.groupby('region_id')['read_support'].transform('sum')
    read_summary_df['percent_support'] = read_summary_df['read_support'] / read_summary_df['total_region_reads']

    allele_dfs = []

    for region_id, group in read_summary_df.groupby('region_id'):
        max_support = group['percent_support'].max()
        if max_support >= major_threshold:
            filtered = group[group['percent_support'] >= major_threshold]
        else:
            filtered = group[group['percent_support'] >= minor_threshold]
        allele_dfs.append(filtered)

    return pd.concat(allele_dfs, ignore_index=True)


def run_pipeline(
    sample_id,
    bam_file,
    motif,
    fasta,
    cluster_distance,
    minor_threshold,
    major_threshold,
    bed_file,
    output_dir
):
    # Define constants
    SAMPLE_ID = sample_id
    BAM_FILE = bam_file
    MOTIF = motif
    FASTA_FILE = fasta
    BED_FILE = bed_file
    OUTPUT_DIR = output_dir

    # Load chromosome name mappings
    with files("src").joinpath("chr_mapping_simple.txt").open("r") as f:
        chr_names = pd.read_csv(f, header=None, sep=' ')
    chr_map = dict(zip(chr_names[1], chr_names[0]))

    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Fetch the reference sequence
    reference_sequence = fetch_reference_sequence(FASTA_FILE)

    # Load BED file with all regions
    regions_df = pd.read_csv(BED_FILE, sep='\t', header=None)
    chromosomes = regions_df[0].unique()

    # Prepare args for multiprocessing
    process_args = []
    for chrom in chromosomes:
        chrom_regions = regions_df[regions_df[0] == chrom].copy()
        process_args.append((chrom, chrom_regions, MOTIF, BAM_FILE, chr_map, reference_sequence, cluster_distance))

    # Run each chromosome in parallel
    n_processes = len(chromosomes)
    pool = mp.Pool(processes=n_processes)
    results = pool.map(process_single_chromosome, process_args)
    pool.close()
    pool.join()

    # Combine results
    results_dict = {chromosome: df for chromosome, df in results if not df.empty}
    combined_results = pd.concat(results_dict.values(), ignore_index=True)

    # Save all read-level results
    output_file_all_reads = os.path.join(OUTPUT_DIR, f'{SAMPLE_ID}_clustered_results.csv')
    combined_results.to_csv(output_file_all_reads, sep=',', index=False)

    # Call alleles from clustered reads
    sample_allele_df = identify_alleles(combined_results, minor_threshold, major_threshold)
    output_file_alleles = os.path.join(OUTPUT_DIR, f'{SAMPLE_ID}_allele_results.csv')
    sample_allele_df.to_csv(output_file_alleles, sep=',', index=False)



