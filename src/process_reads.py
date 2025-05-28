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

def count_total_motifs(sequence, motif):
    """
    Counts the total number of occurrences of a motif in a given sequence (consecutive and non-consecutive).

    Parameters:
        sequence (str): The sequence to search within.
        motif (str): The motif to look fo in the string.

    Returns:
        int: The total number of motif occurrences in the sequence.
    """
    total_repeat_num = sequence.upper().count(motif)
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
            if is_single_nucleotide_variant(seq_string, motif):
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





