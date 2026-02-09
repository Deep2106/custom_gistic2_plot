# =============================================================================
# GISTIC2 Visualization Script - BM vs EMD Comparison
# Multiple Myeloma Copy Number Analysis (Siobhan's lab, RCSI)
# Author : Deepak R. Bharti
# date   : 06/02/2026 
# =============================================================================

# Install packages if not already installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  data.table,
  ggplot2,
  gridExtra,
  grid,
  scales,
  dplyr,
  tidyr,
  ggrepel
)

# =============================================================================
# SET FILE PATHS
# =============================================================================

# Set your GISTIC2 output directories
bm_dir <- "~/RCSI_2025_JOB/clinical_bioinformatician/Literature/catherine/WES_SEG/Copy Number SEG files/COHORT/BM_results/"
emd_dir <- "~/RCSI_2025_JOB/clinical_bioinformatician/Literature/catherine/WES_SEG/Copy Number SEG files/COHORT/EMD_results/"

# File paths
bm_scores_file <- paste0(bm_dir, "scores.gistic")
emd_scores_file <- paste0(emd_dir, "scores.gistic")

# Peak files for gene labeling
bm_amp_genes <- paste0(bm_dir, "amp_genes.conf_90.txt")
bm_del_genes <- paste0(bm_dir, "del_genes.conf_90.txt")
emd_amp_genes <- paste0(emd_dir, "amp_genes.conf_90.txt")
emd_del_genes <- paste0(emd_dir, "del_genes.conf_90.txt")

# Output directory
output_dir <- "~/RCSI_2025_JOB/clinical_bioinformatician/Literature/catherine/WES_SEG/Copy Number SEG files/Results/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# FUNCTION 1: Read GISTIC Scores
# =============================================================================

read_gistic_scores <- function(scores_file) {
  # Read scores.gistic file
  scores <- fread(scores_file, header = TRUE, sep = "\t")
  
  # Clean column names
  colnames(scores) <- c("Type", "Chromosome", "Start", "End", 
                        "q_value", "G_score", "avg_amplitude", "frequency")
  
  # Remove "chr" prefix if present and convert to numeric
  scores$Chromosome <- gsub("chr", "", scores$Chromosome)
  scores$Chromosome <- gsub("X", "23", scores$Chromosome)
  scores$Chromosome <- gsub("Y", "24", scores$Chromosome)
  scores$Chromosome <- as.numeric(scores$Chromosome)
  
  # Keep only autosomes and X (chr 1-23)
  scores <- scores[scores$Chromosome %in% 1:23, ]
  
  # Calculate midpoint for plotting
  scores$Position <- (scores$Start + scores$End) / 2
  
  # Separate Amp and Del
  amp_scores <- scores[scores$Type == "Amp", ]
  del_scores <- scores[scores$Type == "Del", ]
  
  # For deletions, make G-scores negative
  del_scores$G_score <- -abs(del_scores$G_score)
  
  return(list(amp = amp_scores, del = del_scores, combined = scores))
}

# =============================================================================
# FUNCTION 2: Create Chromosome Map (hg38/GRCh38)
# =============================================================================

create_chromosome_map <- function() {
  # Official hg38/GRCh38 chromosome lengths (bp)
  # Source: UCSC Genome Browser, GRCh38/hg38 assembly
  chr_lengths <- c(
    248956422,  # chr1
    242193529,  # chr2
    198295559,  # chr3
    190214555,  # chr4
    181538259,  # chr5
    170805979,  # chr6
    159345973,  # chr7
    145138636,  # chr8
    138394717,  # chr9
    133797422,  # chr10
    135086622,  # chr11
    133275309,  # chr12
    114364328,  # chr13
    107043718,  # chr14
    101991189,  # chr15
    90338345,   # chr16
    83257441,   # chr17
    80373285,   # chr18
    58617616,   # chr19
    64444167,   # chr20
    46709983,   # chr21
    50818468,   # chr22
    156040895   # chrX
  )
  
  # Create cumulative positions for plotting
  chr_map <- data.frame(
    Chromosome = 1:23,
    Length = chr_lengths,
    Start = c(0, cumsum(chr_lengths[-23])),
    End = cumsum(chr_lengths)
  )
  
  # Calculate middle position for chromosome labels
  chr_map$Middle <- (chr_map$Start + chr_map$End) / 2
  
  return(chr_map)
}
# =============================================================================
# FUNCTION 3: Add Genomic Position 
# =============================================================================

add_genomic_position_simple <- function(data, chr_map) {
  # Ensure numeric chromosomes
  data$Chromosome <- as.numeric(as.character(data$Chromosome))
  data <- data[!is.na(data$Chromosome) & data$Chromosome %in% 1:23, ]
  
  cat("Valid chromosomes in data:", paste(unique(data$Chromosome), collapse = ", "), "\n")
  
  if (nrow(data) == 0) {
    warning("No valid chromosomes (1-23) found in data")
    return(data.frame())
  }
  
  # Calculate genomic position row by row
  data$GenomicPosition <- sapply(1:nrow(data), function(i) {
    chr <- data$Chromosome[i]
    pos <- data$Position[i]
    chr_start <- chr_map$Start[chr_map$Chromosome == chr]
    
    if (length(chr_start) == 0) {
      warning(paste("No start position found for chromosome", chr))
      return(NA)
    }
    return(chr_start + pos)
  })
  
  # Remove any NAs
  n_before <- nrow(data)
  data <- data[!is.na(data$GenomicPosition), ]
  n_after <- nrow(data)
  
  if (n_before != n_after) {
    cat("Removed", n_before - n_after, "rows with invalid genomic positions\n")
  }
  
  cat("Successfully calculated", n_after, "genomic positions\n")
  
  return(data)
}

# =============================================================================
# FUNCTION 4: Extract Peak Genes (Transposed Format)
# =============================================================================

extract_peak_genes_transposed <- function(amp_genes_file, del_genes_file, chr_map, 
                                          q_threshold = 0.25, top_n = 10) {
  
  peaks_list <- list()
  
  # Function to parse one file
  parse_peak_file <- function(file_path, type) {
    if (!file.exists(file_path)) {
      cat("File not found:", file_path, "\n")
      return(NULL)
    }
    
    # Read the transposed data
    data <- fread(file_path, header = TRUE)
    
    # Extract row names (first column)
    row_names <- data[[1]]
    
    # Find key rows
    q_row_idx <- which(grepl("^q value$", row_names, ignore.case = TRUE))
    boundary_row_idx <- which(grepl("wide peak boundaries", row_names, ignore.case = TRUE))
    genes_row_idx <- which(grepl("genes in wide peak", row_names, ignore.case = TRUE))
    
    if (length(q_row_idx) == 0 || length(boundary_row_idx) == 0) {
      cat("Required rows not found in:", file_path, "\n")
      return(NULL)
    }
    
    # Get peak names (column names, excluding first column) - THESE ARE CYTOBANDS
    peak_names <- colnames(data)[-1]
    
    # Extract q-values, boundaries, and genes
    q_values <- as.numeric(data[q_row_idx, -1, with = FALSE])
    boundaries <- as.character(data[boundary_row_idx, -1, with = FALSE])
    
    if (length(genes_row_idx) > 0) {
      genes <- as.character(data[genes_row_idx, -1, with = FALSE])
    } else {
      genes <- rep("", length(peak_names))
    }
    
    # Create data frame
    peaks_df <- data.frame(
      Cytoband = peak_names,  # Changed from Peak to Cytoband for clarity
      q_value = q_values,
      Boundary = boundaries,
      Genes = genes,
      Type = type,
      stringsAsFactors = FALSE
    )
    
    # Filter by q-value
    peaks_df <- peaks_df[peaks_df$q_value <= q_threshold & !is.na(peaks_df$q_value), ]
    
    if (nrow(peaks_df) == 0) {
      cat("No significant peaks (q <=", q_threshold, ") in:", file_path, "\n")
      return(NULL)
    }
    
    # Sort by q-value and take top N
    peaks_df <- peaks_df[order(peaks_df$q_value), ]
    peaks_df <- head(peaks_df, top_n)
    
    # Parse chromosome and position from boundaries
    # Format: chr1:122327639-143324886
    peaks_df$Chromosome <- gsub("chr", "", gsub(":.*", "", peaks_df$Boundary))
    
    # Extract start and end positions
    positions <- gsub(".*:", "", peaks_df$Boundary)
    peaks_df$Start <- as.numeric(gsub("-.*", "", positions))
    peaks_df$End <- as.numeric(gsub(".*-", "", positions))
    peaks_df$Position <- (peaks_df$Start + peaks_df$End) / 2
    
    # Clean up gene names (remove brackets and pick first gene if multiple)
    peaks_df$Gene_Symbol <- gsub("\\[|\\]", "", peaks_df$Genes)
    peaks_df$Gene_Symbol <- sapply(strsplit(peaks_df$Gene_Symbol, ",|;"), function(x) {
      x <- trimws(x)
      if (length(x) > 0 && x[1] != "") return(x[1]) else return("")
    })
    
    # PRIMARY LABEL: Use cytoband name
    # SECONDARY: If there's a well-known gene, you can optionally show it
    peaks_df$Display_Label <- peaks_df$Cytoband
    
    # Optional: For certain well-known cancer genes, show gene name instead
    # Uncomment if you want to show gene names for specific genes
    # known_genes <- c("FGFR3", "TP53", "MYC", "CCND1", "RB1", "PTEN")
    # peaks_df$Display_Label <- ifelse(peaks_df$Gene_Symbol %in% known_genes, 
    #                                  peaks_df$Gene_Symbol, 
    #                                  peaks_df$Cytoband)
    
    # Convert chromosome to numeric
    peaks_df$Chromosome <- gsub("X", "23", peaks_df$Chromosome)
    peaks_df$Chromosome <- gsub("Y", "24", peaks_df$Chromosome)
    peaks_df$Chromosome <- as.numeric(peaks_df$Chromosome)
    
    # Keep only autosomes and X
    peaks_df <- peaks_df[peaks_df$Chromosome %in% 1:23 & !is.na(peaks_df$Chromosome), ]
    
    return(peaks_df)
  }
  
  # Parse amplifications
  if (file.exists(amp_genes_file)) {
    cat("Parsing amplification peaks...\n")
    amp_peaks <- parse_peak_file(amp_genes_file, "Amp")
    if (!is.null(amp_peaks) && nrow(amp_peaks) > 0) {
      peaks_list[[1]] <- amp_peaks
      cat("Found", nrow(amp_peaks), "significant amplification peaks\n")
    }
  }
  
  # Parse deletions
  if (file.exists(del_genes_file)) {
    cat("Parsing deletion peaks...\n")
    del_peaks <- parse_peak_file(del_genes_file, "Del")
    if (!is.null(del_peaks) && nrow(del_peaks) > 0) {
      peaks_list[[2]] <- del_peaks
      cat("Found", nrow(del_peaks), "significant deletion peaks\n")
    }
  }
  
  # Combine peaks
  if (length(peaks_list) == 0) {
    cat("WARNING: No significant peaks found\n")
    return(data.frame())
  }
  
  peaks <- do.call(rbind, peaks_list)
  
  # Add genomic position using the simple function
  peaks <- add_genomic_position_simple(peaks, chr_map)
  
  if (nrow(peaks) > 0) {
    cat("Total peaks to label:", nrow(peaks), "\n")
    print(peaks[, c("Display_Label", "Cytoband", "q_value", "Type", "Chromosome")])
  }
  
  return(peaks)
}

# =============================================================================
# FUNCTION 5: Plot GISTIC Scores (With Smart Label Positioning)
# =============================================================================

plot_gistic_scores <- function(scores_data, chr_map, peaks = NULL, 
                               title = "GISTIC2 Scores", 
                               q_threshold = 0.25) {
  
  # Combine amp and del
  plot_data <- rbind(scores_data$amp, scores_data$del)
  plot_data <- add_genomic_position_simple(plot_data, chr_map)
  
  # Separate amp and del for different colors
  amp_data <- plot_data[plot_data$G_score > 0, ]
  del_data <- plot_data[plot_data$G_score < 0, ]
  
  # Calculate y-axis range for chromosome bar sizing
  y_max <- max(plot_data$G_score, na.rm = TRUE)
  y_min <- min(plot_data$G_score, na.rm = TRUE)
  y_range <- y_max - y_min
  
  # Create chromosome marking bars - SHORT BARS near x-axis
  chr_shading <- chr_map
  chr_shading$fill <- ifelse(chr_shading$Chromosome %% 2 == 0, "black", "white")
  
  # Bar height as percentage of total plot height
  bar_height <- y_range * 0.08  # 8% of plot height
  
  # Create base plot
  p <- ggplot() +
    # Add SHORT chromosome bars at the bottom (near y=0)
    geom_rect(data = chr_shading, 
              aes(xmin = Start, xmax = End, 
                  ymin = -bar_height/2, ymax = bar_height/2,
                  fill = fill),
              alpha = 1.0, color = NA) +
    scale_fill_identity() +
    
    # Add AMPLIFICATION bars (red)
    geom_segment(data = amp_data,
                 aes(x = GenomicPosition, xend = GenomicPosition,
                     y = 0, yend = G_score),
                 color = "#D6604D", size = 0.4, alpha = 0.8) +
    
    # Add DELETION bars (blue)
    geom_segment(data = del_data,
                 aes(x = GenomicPosition, xend = GenomicPosition,
                     y = 0, yend = G_score),
                 color = "#4393C3", size = 0.4, alpha = 0.8) +
    
    # Add THICK horizontal line at 0
    geom_hline(yintercept = 0, color = "black", size = 1.5) +
    
    # Add chromosome labels on x-axis
    scale_x_continuous(
      breaks = chr_map$Middle,
      labels = c(1:22, "X"),
      expand = c(0.01, 0)
    ) +
    
    # Y-axis with extra space for labels
    scale_y_continuous(
      name = "G-Score",
      expand = expansion(mult = c(0.15, 0.15))  # 15% extra space top and bottom
    ) +
    
    # Theme
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 12, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.margin = margin(10, 10, 10, 10),
      axis.line.x = element_line(color = "black", size = 0.5),
      axis.line.y = element_line(color = "black", size = 0.5)
    ) +
    
    labs(title = title)
  
  # Add cytoband labels - POSITIONED AT EACH PEAK
  if (!is.null(peaks) && nrow(peaks) > 0) {
    
    # For each peak, find the actual G-score at that position
    peaks_with_scores <- peaks
    
    # Match peaks to their actual G-scores from plot_data
    for (i in 1:nrow(peaks_with_scores)) {
      # Find closest genomic position in plot_data
      closest_idx <- which.min(abs(plot_data$GenomicPosition - peaks_with_scores$GenomicPosition[i]))
      peaks_with_scores$Peak_Gscore[i] <- plot_data$G_score[closest_idx]
    }
    
    # Split peaks into amp and del
    amp_peaks <- peaks_with_scores[peaks_with_scores$Type == "Amp", ]
    del_peaks <- peaks_with_scores[peaks_with_scores$Type == "Del", ]
    
    # Add AMPLIFICATION labels (above each peak)
    if (nrow(amp_peaks) > 0) {
      p <- p + geom_text_repel(
        data = amp_peaks,
        aes(x = GenomicPosition, y = Peak_Gscore, label = Display_Label),
        size = 2.5,
        fontface = "italic",
        color = "black",
        direction = "both",        # Allow both x and y movement
        nudge_y = y_range * 0.1,  # Small upward nudge
        segment.size = 0.3,
        segment.color = "gray50",
        segment.alpha = 0.6,
        min.segment.length = 0.1,
        box.padding = 0.35,
        point.padding = 0.3,
        force = 3,
        force_pull = 0.5,
        max.overlaps = Inf,
        max.iter = 10000,
        seed = 42
      )
    }
    
    # Add DELETION labels (below each peak)
    if (nrow(del_peaks) > 0) {
      p <- p + geom_text_repel(
        data = del_peaks,
        aes(x = GenomicPosition, y = Peak_Gscore, label = Display_Label),
        size = 2.5,
        fontface = "italic",
        color = "black",
        direction = "both",        # Allow both x and y movement
        nudge_y = -y_range * 0.1, # Small downward nudge
        segment.size = 0.3,
        segment.color = "gray50",
        segment.alpha = 0.6,
        min.segment.length = 0.1,
        box.padding = 0.35,
        point.padding = 0.3,
        force = 3,
        force_pull = 0.5,
        max.overlaps = Inf,
        max.iter = 10000,
        seed = 42
      )
    }
  }
  
  return(p)
}
# =============================================================================
# MAIN EXECUTION
# =============================================================================

# Create chromosome map
chr_map <- create_chromosome_map()

# Read BM data
cat("========================================\n")
cat("Processing BM (Bone Marrow) cohort\n")
cat("========================================\n")
bm_scores <- read_gistic_scores(bm_scores_file)
bm_peaks <- extract_peak_genes_transposed(bm_amp_genes, bm_del_genes, chr_map, 
                                          q_threshold = 0.25, top_n = 15)

# Read EMD data
cat("\n========================================\n")
cat("Processing EMD (Extramedullary Disease) cohort\n")
cat("========================================\n")
emd_scores <- read_gistic_scores(emd_scores_file)
emd_peaks <- extract_peak_genes_transposed(emd_amp_genes, emd_del_genes, chr_map, 
                                           q_threshold = 0.25, top_n = 15)

# Create plots
cat("\n========================================\n")
cat("Generating plots...\n")
cat("========================================\n")

p_bm <- plot_gistic_scores(bm_scores, chr_map, peaks = bm_peaks, 
                           title = "BM (Bone Marrow)", q_threshold = 0.25)

p_emd <- plot_gistic_scores(emd_scores, chr_map, peaks = emd_peaks, 
                            title = "EMD (Extramedullary Disease)", q_threshold = 0.25)

# Combine plots vertically (EMD on top, BM on bottom - like your reference)
combined_plot <- grid.arrange(p_emd, p_bm, ncol = 1, heights = c(1, 1))

# Save as high-resolution PDF
pdf(paste0(output_dir, "GISTIC2_BM_vs_EMD_comparison.pdf"), 
    width = 12, height = 8)
grid.arrange(p_emd, p_bm, ncol = 1, heights = c(1, 1))
dev.off()

# Save as PNG
png(paste0(output_dir, "GISTIC2_BM_vs_EMD_comparison.png"), 
    width = 3600, height = 2400, res = 300)
grid.arrange(p_emd, p_bm, ncol = 1, heights = c(1, 1))
dev.off()

# Save individual plots
pdf(paste0(output_dir, "GISTIC2_BM_only.pdf"), width = 12, height = 4)
print(p_bm)
dev.off()

pdf(paste0(output_dir, "GISTIC2_EMD_only.pdf"), width = 12, height = 4)
print(p_emd)
dev.off()

cat("\n========================================\n")
cat("✓ Plots saved successfully!\n")
cat("Check:", output_dir, "\n")
cat("========================================\n")