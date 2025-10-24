import argparse
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
import os
import logging

# Configure logging once at the start of the script
logging.basicConfig(level=logging.INFO, format="%(message)s")


def find_motif_genomic_coordinates(sequence, motif, buffer_size):
    """Find occurrences of a motif in a genomic sequence."""

    motif = Seq(motif)
    motif_str = str(motif)
    motif_len = len(motif)
    motif_start = 0
    res = []

    while True:
        # find successive occurances of motifs as you move along the contig
        motif_start = sequence.find(motif_str, motif_start)
        if motif_start == -1:
            break

        # Initialize forward scanning variables.
        current_pos = motif_start
        motif_num_next = 1
        perfect_count = 1  # the initial match is perfect

        # Scan forward from the first match to find all repeats.
        while True:
            next_chunk_start = current_pos + motif_len
            next_chunk = sequence[next_chunk_start: next_chunk_start + motif_len]
            if next_chunk == motif_str:
                motif_num_next += 1
                perfect_count += 1
                current_pos = next_chunk_start
            else:
                break

        motif_end = current_pos + motif_len
        buffer_start = max(0, motif_start - buffer_size)
        buffer_end = motif_end + buffer_size

        res.append({
            'Motif': motif_str,
            'Motif_Start': motif_start,
            'Motif_End': motif_end,
            'Total_Repeats': motif_num_next,
            'Perfect_Repeats': perfect_count,
            'Buffer_Start': buffer_start,
            'Buffer_End': buffer_end
        })

        motif_start = motif_end

    return res

def merge_adjacent_motifs(df, max_gap):
    if df.empty:
        return df

    df = df.sort_values(by=['chrom', 'Motif_Start']).reset_index(drop=True)
    df['group'] = (df['Motif_Start'] > (df['Motif_End'].shift() + max_gap)).cumsum()

    # determined merged group bounds
    def merge_group(group):
        contig = group['chrom'].iloc[0]
        motif_val = group['motif'].iloc[0]
        buf_start = group['Buffer_Start'].min()
        buf_end = group['Buffer_End'].max()
        motif_start = group['Motif_Start'].min()
        motif_end = group['Motif_End'].max()
        total_repeats = group['Total_Repeats'].sum()
        perfect_repeats = group['Perfect_Repeats'].max()

        return pd.Series({
            'chrom': contig,
            'Buffer_Start': buf_start,
            'Buffer_End': buf_end,
            'motif': motif_val,
            'Motif_Start': motif_start,
            'Motif_End': motif_end,
            'Total_Repeats': total_repeats,
            'Perfect_Repeats': perfect_repeats,
        })

    merged_df = df.groupby(['chrom', 'group'], group_keys=False).apply(merge_group, include_groups=True).reset_index(drop=True)
    return merged_df

def name_each_region(merged_df):
    """Give each region a unique identifier to be used to track regions during post processing"""
    merged_df = merged_df.copy()
    merged_df['region_name'] = merged_df['chrom'] + '_region_'  + (merged_df.reset_index().index + 1).astype(str)
    return merged_df

def run_find(
        motif,
        repeats,
        gap,
        fasta_file,
        buffer,
        output_dir
        ):

    os.makedirs(output_dir, exist_ok=True)

    results = []
    motif = motif.upper()
    perfect_repeats = repeats
    max_gap = gap
    buffer_size = buffer
    chrom_sequences = {}
    
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        chrom = seq_record.id
        sequence = seq_record.seq.upper()
        chrom_sequences[chrom] = sequence
                
        # Look for motif occurances (always store forward motif)
        forward_results = find_motif_genomic_coordinates(sequence, motif, buffer_size)
        for res in forward_results:
            res['chrom'] = chrom
            res['motif'] = motif
        results.extend(forward_results)
        
        # look for reverse complement occurrances (but store forward motif)
        rev_motif = str(Seq(motif).reverse_complement())
        reverse_results = find_motif_genomic_coordinates(sequence, rev_motif, buffer_size)
        for res in reverse_results:
            res['chrom'] = chrom
            res['motif'] = motif  # Store forward motif, not reverse
        results.extend(reverse_results)
    
    # collect all results into a dataframe
    df_all = pd.DataFrame(results)
    
    if df_all.empty:
        logging.warning("No motif instances found in sample.")
        return

    # order column to match bed format
    df_all = df_all[['chrom', 'Buffer_Start', 'Buffer_End','motif', 'Motif_Start', 'Motif_End', 'Total_Repeats', 'Perfect_Repeats']]
    # filter out regions with less than min number of perfect_repeats
    df_all_filtered = df_all[df_all['Perfect_Repeats'] >= perfect_repeats].copy()
    

    # merge repeat regions per chrom that are in close proximity to be one region instead of two
    merged_all = []
    for contig in df_all_filtered['chrom'].unique():
        df_sub = df_all_filtered[df_all_filtered['chrom'] == contig].copy()
        merged = merge_adjacent_motifs(df_sub, max_gap)
        merged_all.append(merged)

    df_merged = pd.concat(merged_all, ignore_index=True)
    df_merged = name_each_region(df_merged) 
    
    df_final = df_merged[['chrom', 'Buffer_Start', 'Buffer_End', 'motif', 'region_name']]
    df_final.to_csv(os.path.join(output_dir, f'microsatellite_coordinates.bed'), sep='\t', index=False, header=False)
