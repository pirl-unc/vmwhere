library(dplyr)
library(ggplot2)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
sample_file <- args[1]
chrom <- args[2]
start_position <- as.integer(args[3])
output_file <- args[4]


generate_region_visualization <- function(sample_file, chrom, start_position, output_file) {
  
  # Check if the file exists
  if (!file.exists(sample_file)) {
    stop(paste("File not found:", sample_file))
  }
  
  # Read the chromosome data
  sample <- read.csv(sample_file, sep = ',')
  
  # Create region ID column
  sample$region_id <- paste0(paste0(sample$chr, ":", sample$start), "-", sample$end)
  
  # Filter for region of interest based on start position
  region <- sample %>%
    filter(chr == chrom) %>%
    filter(start == start_position) %>%
    arrange(desc(cluster_id))
  
  # Check if any regions match the start position
  if (nrow(region) == 0) {
    stop(paste("No regions found with start position:", start_position))
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
  all_segments <- data.frame()
  for (i in 1:nrow(region)) {
    segments <- parse_read_structure(region$read_structure[i])
    segments$read_id <- region$cluster_id[i]
    segments$read_support <- region$read_support[i]
    
    all_segments <- rbind(all_segments, segments)
  }
  
  # Calculate position for each segment block
  all_segments <- all_segments %>%
    group_by(read_id) %>%
    mutate(
      position_start = cumsum(c(0, count[-n()])),
      position_end = cumsum(count)
    ) %>%
    ungroup()
  
  # Add region ID for plot title
  region_id_value <- region$region_id[1]
  all_segments$region_id <- region_id_value
  
  # Create proper factor levels for read_id (ordered by cluster_id descending)
  all_segments$read_id <- factor(all_segments$read_id, levels = region$cluster_id)
  
  # Get unique read support values for y-axis labels
  read_support_labels <- region %>%
    select(cluster_id, read_support) %>%
    arrange(desc(cluster_id)) %>%
    mutate(label = paste0('(', read_support, ')')) %>%
    pull(label)
  
  # Define base color palette for sequences
  base_sequence_colors <- c(
    "GGAA" = "#FFCC66",
    "GGGA" = "gray80",
    "AGAA" = "#66CC99",
    "GGAG" = "mediumpurple",
    "GGAC" = "lavender",
    "GGAT" = "#99FF99",
    "GGCA" = "bisque2",
    "GGTA" = "forestgreen",
    "GAGG" = "indianred",
    "GAAG" = "lightpink",
    "GAAA" = "orangered2",
    "GCAA" = "#CC0066",
    "CGAA" = "deepskyblue",
    "GTAA" = "darkolivegreen",
    "GAG" = "#3366CC",
    "GGA" = "orange",
    "GAA" = "#FF9966",
    "AGA" = "gold4",
    "GCA" = "yellow",
    "AAA" = "#660000",
    "AAAA" = "tomato",
    "AGTA" = "skyblue",
    "TGAT" = "darkorange2",
    "TGAA" = "gray30",
    "AA" = "#0033FF",
    "AC" = "#666699",
    "TG" = "#000033",
    "T" = 'sienna',
    "G" = 'sienna',
    "C" = 'sienna',
    "A" = 'sienna',
    "NONE" = "lightgray",
    "NULL" = "cornflowerblue"
  )
  
  # Create visualization
  p <- ggplot(all_segments, aes(y = read_id)) +
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
    ), size = 4, color = "black") +
    scale_fill_manual(values = base_sequence_colors, na.value = "white") +
    scale_y_discrete(labels = read_support_labels) +
    labs(
      title = region_id_value,
      fill = "Sequence"
    ) +
    theme_void() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 13),
      legend.position = "bottom"
    )
  
  # Return the plot
  ggsave(output_file, plot = p, units = 'in', width = 5, height = 5)
}


generate_region_visualization(sample_file, chrom, start_position, output_file)

