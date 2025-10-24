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

def find_microsatellite_boundaries(sequence, motif):
    """Return start and end of region flanked by ≥2 consecutive repeats of motif at both ends."""
    # Pattern to match ≥2 perfect repeats
    pattern = f'(?:{motif}){{{2},}}'

    # All matches of ≥2 tandem motifs
    matches = list(re.finditer(pattern, sequence))

    if len(matches) == 0:
        # No matches found
        return -1, -1
    elif len(matches) == 1:
        # Only one match, keep it
        return matches[0].start(), matches[-1].end()
    else:  # len(matches) > 1
        # Calculate repeat numbers for each match
        match_repeat_nums = []
        for match in matches:
            repeat_num = (match.end() - match.start()) // len(motif)
            match_repeat_nums.append(repeat_num)
        
        max_repeats = max(match_repeat_nums)
        
        if max_repeats >= 4:
            # Check if any locations with fewer repeats are within 20 bp of another match
            filtered_matches = []
            
            for i, match in enumerate(matches):
                current_repeats = match_repeat_nums[i]
                
                if current_repeats >= 4:
                    # Keep matches with maximum repeats
                    filtered_matches.append(match)
                else:
                    # Check if this match is within 20 bp of any other match
                    within_20bp = False
                    for j, other_match in enumerate(matches):
                        if i != j:  # Don't compare with itself
                            # Calculate distance between matches
                            distance = min(
                                abs(match.start() - other_match.end()),
                                abs(match.end() - other_match.start())
                            )
                            if distance <= 20:
                                within_20bp = True
                                break
                    
                    if within_20bp:
                        filtered_matches.append(match)
            
            matches = filtered_matches
        else:
            # Max repeats < 4, keep all matches since we don't know which are microsat vs whihc are buffer
            pass  # matches remains unchanged
    
    # Return start of first match and end of last match
    if len(matches) == 0:
        return -1, -1
    else:
        return matches[0].start(), matches[-1].end()

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

def decompose_seq_with_motif_anchors(sequence, motif):
    """
    Determines the structure of a complete sequence by finding exact matches to the motif
    and decomposing all regions by minimizing edit distance.

    Parameters:
        sequence: str, complete sequence to decompose
        motif: str, motif of interest (e.g., "GGAA")

    Returns:
        tuple: (structure_string, allele_length, original_sequence)
    """
    if not sequence:
        return None, 0, ''
    
    # Check if motif exists in sequence
    if sequence.find(motif) == -1:
        # No motifs - decompose entire sequence
        if len(sequence) <= 4 or is_homopolymer(sequence):
            return f"1{sequence}", len(sequence), sequence
        else:
            parsed = further_parse_intervening(sequence, motif)
            return parsed, len(sequence), sequence
    
    # Find motif positions as anchors
    motif_len = len(motif)
    motif_start_positions = []
    start = 0
    
    while start <= len(sequence) - motif_len:
        pos = sequence.find(motif, start)
        if pos == -1:
            break
        motif_start_positions.append(pos)
        start = pos + motif_len
    
    # Build structure with all segments decomposed around motifs
    read_structure = []
    current_pos = 0
    idx = 0
    
    while idx < len(motif_start_positions):
        motif_pos = motif_start_positions[idx]
        
        # Process sequence before this motif (prefix or intervening)
        if current_pos < motif_pos:
            prefix_seq = sequence[current_pos:motif_pos]
            prefix_structure = decompose_non_motif_region(prefix_seq, motif)
            if prefix_structure:
                read_structure.append(prefix_structure)
        
        # Count consecutive motifs
        count = 1
        while (idx + count < len(motif_start_positions) and 
               motif_start_positions[idx + count] == motif_pos + count * motif_len):
            count += 1
        
        read_structure.append(f"{count}{motif}")
        current_pos = motif_pos + count * motif_len
        idx += count
    
    # Process remaining sequence (suffix)
    if current_pos < len(sequence):
        suffix_seq = sequence[current_pos:]
        suffix_structure = decompose_non_motif_region(suffix_seq, motif)
        if suffix_structure:
            read_structure.append(suffix_structure)
    
    return "_".join(read_structure), len(sequence), sequence


def decompose_non_motif_region(region_seq, motif):
    """
    Decompose a region that doesn't start with a motif.
    
    Parameters:
        region_seq: str, sequence region to decompose
        motif: str, reference motif
        
    Returns:
        str: Structured representation of the region
    """
    if not region_seq:
        return ""
    
    if len(region_seq) <= 4:
        return f"1{region_seq}"
    elif is_homopolymer(region_seq):
        return f"1{region_seq}"
    else:
        # Use the same parsing logic for any non-motif region
        return further_parse_intervening(region_seq, motif)


def further_parse_intervening(intervening_seq, motif):
    """
    Parse any sequence by finding segments that minimize edit distance to the motif.
    Used for all non-motif regions regardless of position (prefix, intervening, suffix).
    
    Parameters:
        intervening_seq: str, the sequence to parse
        motif: str, the motif to look for
        
    Returns:
        str: Structured representation of the sequence
    """
    n = len(intervening_seq)
    motif_len = len(motif)
    
    gap_penalty = 2

    # Initialize dynamic programming table for minimizing edit distance
    dp = [float('inf')] * (n + 1)
    dp[0] = 0
    
    # Store backtracking information
    prev_pos = [0] * (n + 1)
    segment_len = [0] * (n + 1)
    
    # Fill the dp table
    for i in range(1, n + 1):
        # Try segments of different lengths up to motif length
        for seg_len in range(1, min(i, motif_len) + 1):
            segment = intervening_seq[i-seg_len:i]
            
            # Calculate edit distance
            if seg_len == motif_len:
                edit_dist = hamming_distance(segment, motif)
            else:
                # Short segment - treat as single unit with penalty equal to length + gap penalty
                edit_dist = seg_len + gap_penalty  
            
            # Update if this gives a better solution
            if dp[i-seg_len] + edit_dist < dp[i]:
                dp[i] = dp[i-seg_len] + edit_dist
                prev_pos[i] = i-seg_len
                segment_len[i] = seg_len
    
    # Reconstruct the optimal decomposition
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
            elif is_single_nucleotide_variant(seq_string, motif):
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
    filtered_region_read_results = region_read_results_df[region_read_results_df['DS_READ'].isin(supported_read_structures)]

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
    read_counts = filtered_region_results.groupby(['DS_READ', 'ID']).size().reset_index(name='RS')
    
    # Get unique combinations of read_structure and region_id with their attributes
    consolidated_results = filtered_region_results.drop_duplicates(subset=['DS_READ', 'ID'])
    
    # Merge the count information with the unique read structures
    consolidated_results = consolidated_results.merge(read_counts, on=['DS_READ', 'ID'], how='left')
    
    return consolidated_results


def cluster_consolidated_reads_by_edit_distance(consolidated_results, cluster_dist):
    """
    Clusters read_sequences based on Levenshtein distance and read support.

    Parameters:
        df (pd.DataFrame): Must contain 'read_sequence' and 'read_support'.
        cluster_distance (int): Max allowed distance to join a cluster.

    Returns:
        pd.DataFrame: DataFrame with cluster assignment and updated supporting reads for the cluster
    """
    results_sorted = consolidated_results.sort_values(by='RS', ascending=False).reset_index(drop=True)
    clusters = []  # Each cluster is {'centroid': str, 'members': [int], 'ids': [int]}
    cluster_assignments = [-1] * len(results_sorted)

    for i, row in results_sorted.iterrows():
        sequence = row['ALT']
        best_dist = float('inf')
        best_cluster_idx = -1

        # Compare to each existing cluster centroid
        for idx, cluster in enumerate(clusters):
            dist = Levenshtein.distance(sequence, cluster['centroid'])
            if dist < best_dist and dist <= cluster_dist:
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
    cluster_read_support = clustered.groupby('cluster_id')['RS'].sum()
    # identify the index of the centroid for each cluster
    centroids_idx = clustered.groupby('cluster_id')['RS'].idxmax()
    # create dataframe with only centroids and drop read support column 
    centroids_df = clustered.loc[centroids_idx]
    centroids_df = centroids_df.drop(columns=['RS'])
    # update total read support for each cluster 
    centroids_df_final = centroids_df.merge(cluster_read_support, on='cluster_id', how='left')
    # we no longer care about cluster id so drop it
    centroids_df_final = centroids_df_final.drop(columns=['cluster_id'])
    
    return centroids_df_final

def process_reads_overlapping_regions(chr_map, regions, bamfile, reference_sequence, cluster_distance):
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
    chrom = regions.iloc[0,0] # column zero contains the chromosome in ucsc naming
    # determine which motif we are processing
    motif = regions.iloc[0,4] # column 4 contains the motif

    ## determine ucsc or refseq naming of chromosomes in bam file
    if bamfile.references[0] in chr_map.keys():
        chrom_bam = chrom #chrom needs to stay in ucsc naming 
    else:
        chrom_bam = chr_map.get(chrom) # chrom  needs to get swtiched to ref_seq naming

    ## Determine ucsc or refseq naming of chromosomes in fasta reference dictionary
    if chrom in reference_sequence.keys():
        chrom_reference = chrom
    else:
        chrom_reference = chr_map[chrom] # convert to proper nomenclature
    
    # Initialize an empty dictionary to store results for each chromosome
    chromosome_query_read_results = {
        'CHROM': [],
        'POS': [],
        'ID': [],
        'REF' : [],
        'ALT': [],
        'END': [],
        'MOTIF' : [],
        'AL':[],
        'CN': [],
        'CNM': [],
        'MD' : [],
        'DS_READ': [],
        'DS_REF':[],
        'RS': []
        }

    suss_dels=0
    reverse_motif = str(Seq(motif).reverse_complement())
        
    # analyse one region at a time
    for i in range(len(regions)):
        if chrom == "chrM": #skip chrM
            continue
        
        region_id = regions.iloc[i,3] # column 3 contains the unique region identifier
        start_coord = int(regions.iloc[i, 1])  # Column 1 contains the start position
        end_coord = int(regions.iloc[i, 2])  # Column 2 contains the end position

        primary_reads = 0
        for read in bamfile.fetch(chrom_bam, start_coord, end_coord):
            if not read.is_secondary or not read.is_unmapped:
                primary_reads += 1
                if primary_reads >= 3:
                    break
        
        if primary_reads < 3:
            print(f"Skipping region {region_id} due to insufficient primary reads ({primary_reads})")
            continue

        # initialize an empty dictionary to store results for each region
        region_query_read_results = {
        'CHROM': [],
        'POS': [],
        'ID': [],
        'REF' : [],
        'ALT': [],
        'END': [],
        'MOTIF': [],
        'AL':[],
        'CN' : [],
        'CNM':[],
        'MD' : [],
        'DS_READ': [],
        'DS_REF' : []
        }

        # Fetch reads from BAM file overlapping current region
        for read in bamfile.fetch(chrom_bam, start_coord, end_coord):
            # Only count reads where this is the primary alignment (skip secondary alignments)
            if read.is_secondary or read.is_unmapped:
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
                use_reverse_complement = ref_forward_total < ref_reverse_total

                # if reverse motif is dominant, find where in the reverse complement the motif is located and convert for forward coordinates
                if use_reverse_complement:
                    working_ref_seq = str(Seq(ref_seq).reverse_complement()) 
                    working_query_seq = str(Seq(query_seq).reverse_complement())
                else:
                    working_ref_seq = ref_seq
                    working_query_seq = query_seq

                ref_motif_start, ref_motif_end = find_microsatellite_boundaries(working_ref_seq, motif)
                query_motif_start, query_motif_end = find_microsatellite_boundaries(working_query_seq, motif)

                if use_reverse_complement:
                    # convert revese complement coordinate for forward coordinates
                    if ref_motif_start != -1:
                        ref_seq_motif_start = len(ref_seq) - ref_motif_end
                        ref_seq_motif_end = len(ref_seq) - ref_motif_start
                    else:
                        ref_seq_motif_start, ref_seq_motif_end = -1, -1
                    
                    if query_motif_start != -1:
                        query_seq_motif_start = len(query_seq) - query_motif_end
                        query_seq_motif_end = len(query_seq) - query_motif_start
                    else:
                        query_seq_motif_start, query_seq_motif_end = -1, -1
                else:
                    ref_seq_motif_start, ref_seq_motif_end = ref_motif_start, ref_motif_end
                    query_seq_motif_start, query_seq_motif_end = query_motif_start, query_motif_end

                
                ref_global_start = ref_start + ref_seq_motif_start
                ref_global_end = ref_start + ref_seq_motif_end

                # handle situation where we have a motif in the reference and the read
                if query_seq_motif_start != -1 and query_seq_motif_end != -1:
                    
                    # convert to genomic global coordinates (we were working in local coords from get_aligned_pairs)
                    query_global_start = qstart + query_seq_motif_start
                    query_global_end = qstart + query_seq_motif_end

                    # Try direct mapping first
                    query_mapped_ref_start = None
                    query_mapped_ref_end = None

                    for q_pos, r_pos in aligned_pairs:
                        if q_pos == query_global_start and r_pos is not None:
                            query_mapped_ref_start = r_pos
                        if q_pos == query_global_end - 1 and r_pos is not None:  # -1 because end is exclusive
                            query_mapped_ref_end = r_pos + 1

                    # extract sequences to compare lengths
                    ref_seq_microsat_initial = str(reference_sequence[chrom_reference][ref_global_start:ref_global_end]).upper()
                    query_seq_microsat_initial = read.query_sequence[query_global_start:query_global_end].upper()

                    # apply reverse complement if we need to 
                    if use_reverse_complement:
                        ref_seq_microsat_initial = str(Seq(ref_seq_microsat_initial).reverse_complement())
                        query_seq_microsat_initial = str(Seq(query_seq_microsat_initial).reverse_complement())
                    
                    # calculate preliminary lengths and compare them
                    ref_length = len(ref_seq_microsat_initial)
                    query_length = len(query_seq_microsat_initial)
                    length_diff = abs(ref_length - query_length)

                    # Determind the boundaries based on the difference in lengths to infer what type of expansion or contraction event we have 
                    if length_diff > 0: # diff lengths - either expansion or contraction
                        if query_mapped_ref_start is not None and query_mapped_ref_end is not None:
                            # compare starting and ending positions to help determine mechanism
                            start_diff = abs(ref_global_start - query_mapped_ref_start)
                            end_diff = abs(ref_global_end - query_mapped_ref_end)

                            if start_diff <= 5 or end_diff <= 5: # small difference likely insertion/deletion
                                # use original boundaires since nothing needs to be expanded or contracted
                                final_ref_start = ref_global_start
                                final_ref_end = ref_global_end
                                final_query_start = query_global_start
                                final_query_end = query_global_end
                            else: # large difference likely substitution
                                # use union boundaries to expand/contract to be properly comparing read seq and ref seqs
                                final_ref_start = min(ref_global_start, query_mapped_ref_start)
                                final_ref_end = max(ref_global_end, query_mapped_ref_end)
                                final_query_start = query_global_start
                                final_query_end = query_global_end
                        else: # there is an insetion/deletion at a msat boundary so can't be mapped, just use the boundaries we already found
                            final_ref_start = ref_global_start
                            final_ref_end = ref_global_end
                            final_query_start = query_global_start 
                            final_query_end = query_global_end
                    else: # same length - no expansion or contraction
                        # Use the boundaries we already identified
                        final_ref_start = ref_global_start
                        final_ref_end = ref_global_end
                        final_query_start = query_global_start
                        final_query_end = query_global_end

                else: # query does not have a microsatellite (due to deletion or substitution
                    final_ref_start = ref_global_start
                    final_ref_end = ref_global_end                            

                    # map reference coordinate to query coordiantes
                    final_query_start = None
                    final_query_end = None

                    for q_pos, r_pos in aligned_pairs:
                        if r_pos == final_ref_start and q_pos is not None:
                            final_query_start = q_pos
                        if r_pos == final_ref_end - 1 and q_pos is not None:  # -1 because end is exclusive
                            final_query_end = q_pos + 1
                    
                    # if mapping fails we can just set coordinates that will result in an empty query sequence and we can manually inspect this region in igv later
                    if final_query_start is None or final_query_end is None:
                        final_query_start = 0
                        final_query_end = 0

            
                # Extract final sequences
                query_seq_microsat = read.query_sequence[final_query_start:final_query_end]
                ref_seq_microsat = str(reference_sequence[chrom_reference][final_ref_start:final_ref_end]).upper()

                # Apply reverse complement if needed
                if use_reverse_complement:
                    query_seq_microsat = str(Seq(query_seq_microsat).reverse_complement())
                    ref_seq_microsat = str(Seq(ref_seq_microsat).reverse_complement())
    
            else:
                continue # skip read since it doesn't cover the whole region

            # get decomposition of the reference
            ref_structure, allele_length_ref, ref_microsat_seq = decompose_seq_with_motif_anchors(ref_seq_microsat, motif)
                
            # perform decomposition of the read sequence and identify variants present 
            read_structure, allele_length, read_microsat_seq = decompose_seq_with_motif_anchors(query_seq_microsat, motif)
            perfect_motif_count_decomp, count_max_consecutive_motif_repeats_decomp, snp_variant_count, multi_bp_variant_count, single_bp_inbetween_count, multi_bp_inbetween_count = classify_variants_from_structure(read_structure, motif)

            # determine total number of prefect motifs present and max consecutive by counting them in the original string
            max_consecutive_repeats = count_max_consecutive_motifs(query_seq_microsat, motif)
            total_repeats = count_total_motifs(query_seq_microsat, motif)

            # compare motifs found in decomposed sequence to those found in the original sequence to ensure the decompositions is correct
            assert perfect_motif_count_decomp == total_repeats, f"Decomposed motif count {perfect_motif_count_decomp} does not match original sequence motif count of {total_repeats} in {read_structure} in {region_id} for read sequence {query_seq_microsat}"
            assert count_max_consecutive_motif_repeats_decomp == max_consecutive_repeats, f"Decomposed motif count {count_max_consecutive_motif_repeats_decomp} does not match original sequence motif count of {max_consecutive_repeats} in {read_structure} in {region_id} for read sequence {query_seq_microsat}"
                
            # recomposed the decomposed sequence and check to ensure it matches the original sequence that was decomposed
            if read_structure is not None:
                reconstructed_seq_check = recompose_string_from_structure(read_structure)
                assert reconstructed_seq_check == query_seq_microsat, f"Recomposed sequence {reconstructed_seq_check} from decomposed sequence {read_structure} does not match original sequence {query_seq_microsat} in read {read.query_name} for region {region_id}"
                motif_density = (total_repeats*len(motif))/allele_length
            else:
                motif_density = 0

            # record read characteristics
            region_query_read_results['CHROM'].append(chrom)
            region_query_read_results['POS'].append(start_coord)
            region_query_read_results['ID'].append(region_id)
            region_query_read_results['REF'].append(ref_seq_microsat)
            region_query_read_results['ALT'].append(read_microsat_seq)
            region_query_read_results['MOTIF'].append(motif)
            region_query_read_results['END'].append(end_coord) 
            region_query_read_results['DS_READ'].append(read_structure if read_structure is not None else "<NONE>")
            region_query_read_results['CNM'].append(max_consecutive_repeats)
            region_query_read_results['CN'].append(total_repeats)
            region_query_read_results['AL'].append(allele_length)
            region_query_read_results['MD'].append(motif_density)
            region_query_read_results['DS_REF'].append(ref_structure)

        # collect all read_structures for the region
        motif_occurance_list = region_query_read_results['DS_READ']
        # turn region result into dataframe that we will use to filter low support reads
        region_query_read_results_df = pd.DataFrame(region_query_read_results)
        
        # filter out singleton reads by checking how many times the decomposed sequence occurs in the list of reads for that region
        filtered_region_read_results = filter_reads_with_low_support(motif_occurance_list, region_query_read_results_df)
        
        if len(filtered_region_read_results) == 0:
            print(f"Microsatellite {region_id} skipped due to no concordance in reads")
            continue
        else:
            # consolidate read results for each region where each row is a unique read structure and the read support is counted 
            consolidated_region_read_results = consolidate_reads_per_region(filtered_region_read_results)

            # cluster reads that are within the specific edit distance to each other
            clustered_region_read_results = cluster_consolidated_reads_by_edit_distance(consolidated_region_read_results, cluster_distance)

            # add filtered region results to the complete chromosome dictionary results
            for column in clustered_region_read_results.columns:
                chromosome_query_read_results[column].extend(clustered_region_read_results[column].tolist())

    # provide a summary output for each chromosome in the log file
    print(f"Checked {len(regions)} regions on {chrom}")
    print(f"Successfully genotyped {len(chromosome_query_read_results['CHROM'])} microsatellites on {chrom}")

    return pd.DataFrame(chromosome_query_read_results)

def process_single_chromosome(args):
    """
    Process a single chromosome at a time, with outputs coming from the function that processes each region. This function will be called by the multiprocessing pool so each chromosome is run in parallel.
    
    Parameters:
        args (tuple): Contains (chromosome that is processed, regions on that chromosome, reference_sequence)
        
    Returns:
        str: Name of the output file that was created for that chromosome
    """
    motif, chromosome_curr, regions_df, BAM_FILE, chr_map, reference_sequence, cluster_distance = args

    # Open BAM file
    try:
        bam = pysam.AlignmentFile(BAM_FILE, "rb")
    except FileNotFoundError:
        print(f"Error: BAM File {BAM_FILE} not found")
        return motif, chromosome_curr, None
    
    # Process the regions on the chromosome
    query_read_results_df = process_reads_overlapping_regions(
        chr_map=chr_map,
        regions=regions_df,
        bamfile=bam,
        reference_sequence=reference_sequence,
        cluster_distance=cluster_distance
    )
    
    # Close the BAM file
    bam.close()

    # return the chromosome df
    return motif, chromosome_curr, query_read_results_df


def identify_alleles(read_summary_df, minor_thresh=0.20, homozygous_thresh=0.8):
    """Determine alleles based on read support."""

    read_summary_df['RC'] = read_summary_df.groupby('ID')['RS'].transform('sum')
    read_summary_df['PS'] = read_summary_df['RS'] / read_summary_df['RC']

    allele_dfs = []

    for region_id, group in read_summary_df.groupby('ID'):
        max_support = group['PS'].max()
        if max_support >= homozygous_thresh:
            filtered = group[group['PS'] >= homozygous_thresh]
        else:
            filtered = group[(group['PS'] >= minor_thresh)]
        allele_dfs.append(filtered)

    return pd.concat(allele_dfs, ignore_index=True)


def create_tsv_per_sample(sample_allele_df):
    """Create a TSV file for each sample with allele information."""
    
    conv_dict = {
        'CHROM': [],
        'POS': [],
        'ID': [],
        'REF' : [],
        'ALT': [],
        'END': [],
        'MOTIF': [],
        'GT' : [],
        'AL':[],
        'CN' : [],
        'CNM':[],
        'MD' : [],
        'DS_READ': [],
        'DS_REF' : [],
        'RS': [],
        }
    
    uniq_regions = sample_allele_df['ID'].unique()

    for region in uniq_regions:

        region_df = sample_allele_df[sample_allele_df['ID'] == region]

        conv_dict['CHROM'].append(region_df['CHROM'].iloc[0])
        conv_dict['POS'].append(region_df['POS'].iloc[0])
        conv_dict['ID'].append(region)
        conv_dict['REF'].append(region_df['REF'].iloc[0])
        conv_dict['END'].append(region_df['END'].iloc[0])
        conv_dict['MOTIF'].append(region_df['MOTIF'].iloc[0])
        conv_dict['DS_REF'].append(region_df['DS_REF'].iloc[0])
        
        # Add comma-separated values for allele-specific metrics
        for metric in ['AL', 'CN', 'CNM', 'MD', 'DS_READ', 'RS']:
            values = region_df[metric].astype(str).tolist()
            if len(values) == 1:
                values = values * 2
            conv_dict[metric].append(','.join(values))

        # handle GT field
        ref_allele = region_df['REF'].iloc[0]
        observed_alleles = region_df['ALT'].tolist()

        # create mapping index for each allele (will help us build GT field)
        all_possible_alleles = [ref_allele] + [allele for allele in observed_alleles if allele != ref_allele]
        allele_to_idx = {allele: idx for idx, allele in enumerate(all_possible_alleles)}

        # build GT based on which alleles are present and their indices
        gt_indices = [str(allele_to_idx[allele]) for allele in observed_alleles]
        if len(gt_indices) == 1:
            gt_indices = gt_indices * 2
        gt = '/'.join(gt_indices)
        conv_dict['GT'].append(gt)
    
        # handle ALT field 
        alt_alleles = [allele for allele in observed_alleles if allele != ref_allele]
        if len(alt_alleles) == 0:
            conv_dict['ALT'].append('.')
          
        else:
            conv_dict['ALT'].append(','.join(alt_alleles))
            

    return pd.DataFrame(conv_dict)


def run_genotyper(
    sample_id,
    bed_file,
    bam_file,
    fasta,
    cluster_distance,
    minor_threhold,
    major_threshold,
    output_dir,
    num_processes
):
    """
    Main function to run the genotyper.
    """
    # Define constants
    SAMPLE_ID = sample_id
    BAM_FILE = bam_file
    BED_FILE = bed_file
    FASTA_FILE = fasta
    OUTPUT_DIR = output_dir
    NUM_PROCESSES = num_processes
    CLUSTER_DISTANCE = cluster_distance

    # load chromosome name mappings
    with files("vmwhere").joinpath("chr_mapping_simple.txt").open("r") as f:
        chr_names = pd.read_csv(f, header=None, sep=' ')
    chr_map = dict(zip(chr_names[1], chr_names[0]))

    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # fetch the reference sequence
    reference_sequence = fetch_reference_sequence(FASTA_FILE)

    # Load BED file with all regions
    regions_df = pd.read_csv(BED_FILE, sep='\t', header=None)
    chromosomes = regions_df[0].unique()

    # determine all motifs present in the bed file
    all_motifs = regions_df[4].unique()
    all_results = []

    ### set up multiprocessing - collect all tasks first
    all_process_args = []
    
    for motif in all_motifs:
        motif_df = regions_df[regions_df[4] == motif].copy()
        chromosomes = motif_df[0].unique()
        
        # add all chromosome tasks for this motif
        for chrom in chromosomes:
            chrom_regions = motif_df[motif_df[0] == chrom].copy()
            all_process_args.append((motif, chrom, chrom_regions, BAM_FILE, chr_map, reference_sequence, CLUSTER_DISTANCE))

    # use number of total chromosome tasks as process count or to user define processes
    n_processes = min(len(all_process_args), NUM_PROCESSES)

    with mp.Pool(processes=n_processes) as pool:
        # process ALL tasks in one go
        results = pool.map(process_single_chromosome, all_process_args)
        
        # group results by motif and combine
        motif_results = {}
        for motif, chrom, df in results:
            if df is not None and not df.empty:
                if motif not in motif_results:
                    motif_results[motif] = []
                motif_results[motif].append(df)
        
        # combine per-motif results
        for motif, dfs in motif_results.items():
            all_results.append(pd.concat(dfs, ignore_index=True))

    # combine all everything together one sample-wide df
    combined_results = pd.concat(all_results, ignore_index=True)
    
    # perform allele calling 
    sample_allele_df = identify_alleles(combined_results, minor_threhold, major_threshold)

    # convert to be a more appropriate format for tsv (region per row with allele separated values instead of allele per row)
    tsv_df = create_tsv_per_sample(sample_allele_df)
    
    # save TSV
    output_TSV = os.path.join(OUTPUT_DIR, f'{SAMPLE_ID}_vmwhere_results.tsv')
    tsv_df.to_csv(output_TSV, sep='\t', index=False)

    print("Genotyping complete")

