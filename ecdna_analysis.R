# ecDNA Repeat Element Enrichment Analysis Library
# Main analysis functions for regioneR-based permutation testing
# Author: Your Name
# Version: 1.0

# Load required libraries (suppress messages for cleaner output)
suppressMessages({
  library(regioneR)
  library(GenomicRanges)
  library(rtracklayer)
  library(parallel)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
})

#' Complete ecDNA repeat enrichment analysis using regioneR
#'
#' @param ecdna_bed Path to BED file with ecDNA regions
#' @param repeat_bed Path to RepeatMasker BED file  
#' @param genome Genome version (e.g., "hg38", "hg19")
#' @param n_permutations Number of permutations (default: 10000)
#' @param n_cores Number of cores for parallel processing
#' @param randomization_method Method for generating random regions
#' @return List containing results, plots, and statistics
analyze_ecdna_repeats <- function(ecdna_bed, repeat_bed, genome = "hg38", 
                                 n_permutations = 10000, n_cores = 4,
                                 randomization_method = "circularRandomizeRegions") {
  
  cat("=== ecDNA Repeat Element Enrichment Analysis ===\n")
  cat("Using regioneR for robust permutation testing\n\n")
  
  # 1. Load genomic data
  cat("Loading genomic data...\n")
  
  # Load ecDNA regions
  ecdna_gr <- import(ecdna_bed, format = "BED")
  cat("Loaded", length(ecdna_gr), "ecDNA regions\n")
  
  # Load repeat elements
  repeat_gr <- import(repeat_bed, format = "BED")
  cat("Loaded", length(repeat_gr), "repeat elements\n")
  
  # Get genome information
  if (genome == "hg38") {
    genome_gr <- getGenome("hg38")
  } else if (genome == "hg19") {
    genome_gr <- getGenome("hg19")
  } else {
    stop("Supported genomes: hg38, hg19")
  }
  
  # 2. Categorize repeat elements using RepeatMasker family names
  cat("Categorizing repeat elements using RepeatMasker family annotations...\n")
  
  # Your RepeatMasker data format: chr, start, end, repFamily, score, strand
  # The 4th column contains the repeat family name (e.g., AluSp, L1PA10, etc.)
  
  categorize_repeats_by_family <- function(repeat_gr) {
    # Extract repeat family names from the name column (4th column in BED)
    rep_families <- mcols(repeat_gr)$name
    
    # Initialize category vector
    categories <- rep("Other", length(repeat_gr))
    
    # SINE families
    sine_families <- c(
      # Alu subfamilies
      "AluSp", "AluSq", "AluSq2", "AluSx", "AluSx1", "AluSx3", "AluSx4", "AluSz", "AluSz6",
      "AluY", "AluYa5", "AluYa8", "AluYb8", "AluYb9", "AluYc", "AluYc1", "AluYd2", "AluYd6",
      "AluYe5", "AluYf1", "AluYg6", "AluYh3", "AluYh9", "AluYi6", "AluYj4", "AluYk11", "AluYk2",
      # MIR subfamilies
      "MIR", "MIR3", "MIRb", "MIRc", "MIR_Amn",
      # Other SINE families
      "B2", "B4", "ID", "FLAM_A", "FLAM_C"
    )
    sine_mask <- rep_families %in% sine_families | grepl("^Alu", rep_families) | grepl("^MIR", rep_families)
    categories[sine_mask] <- "SINE"
    
    # LINE families
    line_families <- c(
      # L1 subfamilies
      "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6", "L1PA7", "L1PA8", "L1PA10", "L1PA11", 
      "L1PA12", "L1PA13", "L1PA14", "L1PA15", "L1PA16", "L1PA17",
      "L1PB1", "L1PB2", "L1PB3", "L1PB4", "L1PBa",
      "L1M1", "L1M2", "L1M3", "L1M4", "L1M5", "L1MA1", "L1MA2", "L1MA3", "L1MA4", "L1MA5",
      "L1MA6", "L1MA7", "L1MA8", "L1MA9", "L1MA10", "L1MB1", "L1MB2", "L1MB3", "L1MB4", 
      "L1MB5", "L1MB7", "L1MB8", "L1MC", "L1MC1", "L1MC2", "L1MC3", "L1MC4", "L1MC5",
      "L1MD", "L1MD1", "L1MD2", "L1MD3", "L1ME1", "L1ME2", "L1ME3", "L1ME3A", "L1ME4a",
      "L1ME4b", "L1ME5", "L1MF", "L1MF1", "L1MF2", "L1MF3",
      # L2 subfamilies
      "L2", "L2a", "L2b", "L2c",
      # Other LINE families
      "CR1", "RTE1", "RTE-BovB", "Penelope"
    )
    line_mask <- rep_families %in% line_families | grepl("^L1", rep_families) | grepl("^L2", rep_families)
    categories[line_mask] <- "LINE"
    
    # LTR families
    ltr_families <- c(
      # ERV families
      "ERVL", "ERVL-E-int", "ERVL-B4-int", "ERVL-MaLR", "ERVLB4-int", "ERVLE-int",
      "ERV1", "ERV3", "ERV9", "ERVFRD-1", "ERVH", "ERVH48-1", "ERVK", "ERVK3-int",
      "ERVK-int", "ERVK9-int", "ERVK11-int", "ERVK13-int", "ERVKB1-int",
      "ERVW1", "ERVFH21-int",
      # MLT families
      "MLT1A", "MLT1A0", "MLT1A1", "MLT1B", "MLT1C", "MLT1D", "MLT1E", "MLT1E1",
      "MLT1E1A", "MLT1E2", "MLT1E3", "MLT1F", "MLT1F1", "MLT1F2", "MLT1G", "MLT1G1",
      "MLT1G3", "MLT1H", "MLT1H1", "MLT1H2", "MLT1I", "MLT1J", "MLT1J1", "MLT1J2",
      "MLT1K", "MLT1L", "MLT1M", "MLT1N2", "MLT1O", "MLT2A1", "MLT2B1", "MLT2B2",
      "MLT2B3", "MLT2B4", "MLT2B5", "MLT2C1", "MLT2C2", "MLT2D", "MLT2E", "MLT2F",
      # THE families
      "THE1A", "THE1B", "THE1C", "THE1D",
      # Other LTR families
      "HUERS-P1", "HUERS-P2", "HUERS-P3", "LTR1", "LTR2", "LTR3", "LTR4", "LTR5",
      "LTR7", "LTR12", "LTR13", "LTR16", "LTR25", "LTR33", "LTR40", "LTR41"
    )
    ltr_mask <- rep_families %in% ltr_families | grepl("^ERV", rep_families) | 
                grepl("^MLT", rep_families) | grepl("^THE1", rep_families)
    categories[ltr_mask] <- "LTR"
    
    # DNA transposon families
    dna_families <- c(
      # hAT families
      "Charlie1", "Charlie2", "Charlie3", "Charlie4", "Charlie5", "Charlie6", "Charlie7",
      "Charlie8", "Charlie9", "Charlie10", "Charlie11", "Charlie12", "Charlie13", "Charlie14",
      "Charlie15", "Charlie16", "Charlie17", "Charlie18", "Charlie19", "Charlie20",
      # TcMar families
      "Tc1", "Tc2", "Tc4", "Tigger1", "Tigger2", "Tigger3", "Tigger4", "Tigger5", "Tigger6",
      "Tigger7", "Tigger8", "Tigger9", "Tigger10", "Tigger11", "Tigger12", "Tigger13",
      "Tigger14", "Tigger15", "Tigger16", "Tigger17", "Tigger18", "Tigger19", "Tigger20",
      "Tigger21", "Tigger22", "Tigger23", "Tigger24",
      # Merlin/Tip100 families
      "Merlin1", "Tip100", "Tip100a", "Tip100b",
      # MER families (various DNA transposons)
      "MER1", "MER2", "MER4", "MER5", "MER20", "MER21", "MER34", "MER41", "MER44", "MER45",
      "MER46", "MER48", "MER49", "MER50", "MER51", "MER52", "MER53", "MER54", "MER55",
      "MER56", "MER57", "MER58", "MER61", "MER63", "MER65", "MER66", "MER67", "MER68",
      "MER70", "MER71", "MER72", "MER73", "MER74", "MER75", "MER81", "MER82", "MER83",
      "MER84", "MER85", "MER87", "MER91", "MER92", "MER95", "MER96", "MER97", "MER98",
      "MER99", "MER101", "MER103", "MER104", "MER105", "MER106", "MER107", "MER108",
      "MER110", "MER111", "MER112", "MER113", "MER114", "MER115", "MER116", "MER117",
      "MER119", "MER120", "MER121", "MER130", "MER131", "MER132", "MER133", "MER134",
      "MER136", "MER137",
      # Other DNA families
      "HSMAR1", "HSMAR2", "RICKSHA", "DNA6-1_DR", "DNA7-1_DR"
    )
    dna_mask <- rep_families %in% dna_families | grepl("^Charlie", rep_families) | 
                grepl("^Tigger", rep_families) | grepl("^MER[0-9]", rep_families) | 
                grepl("^HSMAR", rep_families)
    categories[dna_mask] <- "DNA"
    
    # Simple repeats and low complexity
    simple_families <- c(
      "(GAATG)n", "(CATTC)n", "(TTTTT)n", "(AAAAT)n", "(A)n", "(T)n", "(AT)n", "(TA)n",
      "(CA)n", "(TG)n", "(GAA)n", "(TTC)n", "(TAAA)n", "(TTTA)n", "(GAAA)n", "(TTTC)n",
      "(GAATG)n", "(CATTC)n", "(ATTCC)n", "(
