library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(patchwork)
library(scales)


args <- commandArgs(trailingOnly = TRUE)
genotype_tsv <- args[1]
microsatellite_id <- args[2]
min_AF <- as.integer(args[3])
output_pdf <- args[4]

generate_region_visualization <- function(sample_file, region_id, min_allele_count, output) {
  
  # Check if the file exists
  if (!file.exists(sample_file)) {
    stop(paste("File not found:", sample_file))
  }
  
  # Read the chromosome data
  initial_df <- read.csv(sample_file, sep = '\t')
  
  region_all <- initial_df %>%
    filter(ID == region_id) %>%
    select(CHROM, POS, ID, DS_READ, RS, AL) %>%
    separate_rows(DS_READ, RS, AL, sep = ",", convert = TRUE) %>% 
    transmute(
      CHROM = CHROM,
      POS = POS,
      ID  = ID,
      AL = AL,
      DS_READ = DS_READ,
      RS  = RS
    )
  
  single_alleles <- region_all %>%
    count(DS_READ, name = "allele_count") %>% 
    filter(allele_count < min_allele_count) %>% 
    pull(DS_READ)
  
  region_all_filtered <- region_all %>% 
    filter(!DS_READ %in% single_alleles)
  
  region <- region_all_filtered %>% 
    distinct(DS_READ, .keep_all = TRUE) %>% 
    arrange(AL) %>%
    mutate(ordered_structure_id = row_number())
  
  
  total_alleles <- nrow(region_all_filtered)
  
  read_structure_allele_frequency_df <- region_all_filtered %>%
    count(DS_READ, name = "count") %>%
    mutate(
      frequency = paste0('(', count, '/', total_alleles, ')'),
      frequency_numeric = count / total_alleles
    ) %>%
    select(DS_READ, frequency, frequency_numeric, count)
  
  # Check if any regions match the region id
  if (nrow(region) == 0) {
    stop(paste("No regions found with region id:", region_id))
  }
  
  # Function to parse read structure into components
  parse_read_structure <- function(structure_string) {
    # Handle <NONE> case
    if (structure_string == "<NONE>" || is.na(structure_string) || structure_string == "") {
      result <- data.frame(
        segment = "<NONE>",
        count = 1,
        sequence = "NONE",
        segment_id = 1
      )
      return(result)
    }
    
    segments <- strsplit(structure_string, "_")[[1]]
    result <- data.frame(segment = segments) %>%
      mutate(
        count = as.numeric(str_extract(segment, "^\\d+")),
        count = if_else(is.na(count), 1, count),
        sequence = str_extract(segment, "[A-Z]+$"),
        segment_id = row_number()
      )
    return(result)
  }
  
  # Process each of the unique reads
  read_segments <- data.frame()
  for (i in 1:nrow(region)) {
    segments <- parse_read_structure(region$DS_READ[i])
    segments$read_id <- region$ordered_structure_id[i]
    
    # retrieve frequency using match
    freq_idx <- match(region$DS_READ[i], read_structure_allele_frequency_df$DS_READ)
    freq_value <- if (!is.na(freq_idx)) read_structure_allele_frequency_df$frequency[freq_idx] else NA
    freq_numeric <- if (!is.na(freq_idx)) read_structure_allele_frequency_df$frequency_numeric[freq_idx] else NA
    count_value <- if (!is.na(freq_idx)) read_structure_allele_frequency_df$count[freq_idx] else NA
    
    segments$frequency <- freq_value
    segments$frequency_numeric <- freq_numeric
    segments$allele_count <- count_value
    
    read_segments <- rbind(read_segments, segments)
  }
  
  # Calculate position for each segment block
  read_segments <- read_segments %>%
    group_by(read_id) %>%
    mutate(
      position_start = cumsum(c(0, count[-n()])),
      position_end = cumsum(count)
    ) %>%
    ungroup()
  
  # Add region ID for plot title
  region_id_value <- region$ID[1]
  read_segments$region_id <- region_id_value
  
  read_order = read_segments %>% arrange(read_id) %>% pull(read_id) %>% unique()
  read_segments$read_id <- factor(read_segments$read_id, levels = read_order)
  
  frequency_labels <- read_segments %>% group_by(read_id) %>% slice_head() %>% pull(frequency) 
  
  # Dynamic color assignment based on motifs in data
  unique_sequences <- unique(read_segments$sequence)
  
  # Count frequency of each sequence across all segments (weighted by count)
  sequence_frequency <- read_segments %>%
    group_by(sequence) %>%
    summarise(total_count = sum(count), .groups = 'drop') %>%
    arrange(desc(total_count))
  
  # Find most frequent 4-mer
  four_mers <- sequence_frequency %>% 
    filter(nchar(sequence) == 4 & !sequence %in% c("NONE", "NULL"))
  most_frequent_4mer <- if(nrow(four_mers) > 0) four_mers$sequence[1] else NULL
  
  # Create base color palette for assignment
  available_colors <- c(
    "#66CC99", "bisque4", "orange3", "cadetblue", "mediumpurple", "lavender", "#99FF99", "indianred", 
    "forestgreen", "bisque2", "lightpink", "orangered2", "deepskyblue", "purple3", "darkolivegreen", 
    "#3366CC", "azure3", "#FF9966", "cornsilk3", "gold4", "#660000", "tomato", "skyblue", "#0033FF", "#666699"
  )
  
  # Initialize color assignments
  base_sequence_colors <- c()
  
  # Fixed assignments
  base_sequence_colors["NONE"] <- "white"
  base_sequence_colors["NULL"] <- "cornflowerblue"
  
  # Single bases get different sienna/brown shades
  single_bases <- c("A", "T", "G", "C")
  sienna_colors <- c("A" = "sienna1", "T" = "sienna2", "C" = "sienna3", "G" = "sienna4")
  for(base in single_bases) {
    if(base %in% unique_sequences) {
      base_sequence_colors[base] <- sienna_colors[base]
    }
  }
  
  # Most frequent 4-mer gets special color
  if(!is.null(most_frequent_4mer)) {
    base_sequence_colors[most_frequent_4mer] <- "#FFCC66"
  }
  
  # Assign colors to remaining sequences
  remaining_sequences <- setdiff(unique_sequences, names(base_sequence_colors))
  color_index <- 1
  
  for(seq in remaining_sequences) {
    if(color_index <= length(available_colors)) {
      base_sequence_colors[seq] <- available_colors[color_index]
      color_index <- color_index + 1
    } else {
      # Fallback to a default color if we run out
      base_sequence_colors[seq] <- "lightgray"
    }
  }
  
  # Create main visualization
  main_plot <- ggplot(read_segments, aes(y = read_id)) +
    geom_rect(aes(
      xmin = position_start,
      xmax = position_end,
      ymin = as.numeric(read_id) - 0.4,
      ymax = as.numeric(read_id) + 0.4,
      fill = sequence
    )) +
    geom_text(aes(
      x = position_start + (count/2),
      y = as.numeric(read_id),
      label = ifelse(sequence == "NONE", "None", 
                     ifelse(count >= 1, as.character(count), ""))
    ), size = 3, color = "black") +
    scale_fill_manual(values = base_sequence_colors, na.value = "white") +
    scale_y_discrete(labels = frequency_labels) +
    labs(
      title = region_id_value,
      fill = "Sequence"
    ) +
    theme_void() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 9),
      legend.position = "bottom"
    )
  
  # Create frequency visualization data
  freq_data <- read_segments %>%
    group_by(read_id) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(read_id, frequency_numeric, allele_count)
  
  # Bar width proportional to frequency
  main_plot_with_freq <- ggplot(read_segments, aes(y = read_id)) +
    # Add frequency bars on the left side
    geom_rect(data = freq_data, aes(
      xmin = -max(read_segments$position_end) * 0.4 * frequency_numeric,
      xmax = 0,
      ymin = as.numeric(read_id) - 0.4,
      ymax = as.numeric(read_id) + 0.4
    ), fill = "steelblue4", alpha = 0.6) +
    # motif sequence colored rectangles
    geom_rect(aes(
      xmin = position_start,
      xmax = position_end,
      ymin = as.numeric(read_id) - 0.4,
      ymax = as.numeric(read_id) + 0.4,
      fill = sequence
    )) +
    # motif count text labels
    geom_text(aes(
      x = position_start + (count/2),
      y = as.numeric(read_id),
      label = ifelse(sequence == "NONE", "None", 
                     ifelse(count >= 1, as.character(count), ""))
    ), size = 3, color = "black") +
    # Frequency percentage labels to the left of the bars
    geom_text(data = freq_data, aes(
      x = -max(read_segments$position_end) * 0.3 * frequency_numeric - 1,
      y = as.numeric(read_id),
      label = ifelse(frequency_numeric > 0.1, paste0(round(frequency_numeric * 100, 1), "%"), NA)),
      size = 3.5, color = "black", hjust = 1, na.rm = T) +
    scale_fill_manual(values = base_sequence_colors, na.value = "white") +
    scale_y_discrete(labels = rep("", length(read_order))) +
    labs(
      #title = paste(region_id_value, "- Frequency bars show relative abundance"),
      fill = "Sequence"
    ) +
    theme_void() +
    expand_limits(x = -max(read_segments$position_end) * 0.15) +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "bottom"
    )
  
  ggsave(output, plot = main_plot_with_freq, units = 'in', width = 5, height = 5)

}

generate_region_visualization(genotype_tsv, microsatellite_id, min_AF, output_pdf)

