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
      "(GAATG)n", "(CATTC)n", "(ATTCC)n", "(GGAAT)n",
      "ALR/Alpha", "(CTTTTT)n", "(AAAAAG)n", "(CAAAAA)n", "(TTTTTT)n"
    )
    simple_mask <- rep_families %in% simple_families | grepl("\\([ATCG]+\\)n$", rep_families) |
                   grepl("^ALR", rep_families)
    categories[simple_mask] <- "Simple"
    
    # Satellite DNA
    satellite_families <- c(
      "HSATII", "SATR1", "SATR2", "BSR/Beta", "GSATII", "ACRO1", "SUBTEL1"
    )
    satellite_mask <- rep_families %in% satellite_families | grepl("^HSAT", rep_families) |
                     grepl("^GSAT", rep_families) | grepl("SAT", rep_families)
    categories[satellite_mask] <- "Satellite"
    
    # RNA genes
    rna_families <- c(
      "U1", "U2", "U3", "U4", "U5", "U6", "U7", "U8", "U11", "U12",
      "7SK", "7SL", "5S", "tRNA", "rRNA", "scRNA", "snRNA", "srpRNA"
    )
    rna_mask <- rep_families %in% rna_families | grepl("^U[0-9]", rep_families) |
                grepl("RNA$", rep_families)
    categories[rna_mask] <- "RNA"
    
    return(categories)
  }
  
  # Add category information to repeat elements
  mcols(repeat_gr)$category <- categorize_repeats_by_family(repeat_gr)
  
  # Create separate GRanges for each category
  repeat_categories <- split(repeat_gr, mcols(repeat_gr)$category)
  
  cat("Repeat categories found:\n")
  for (cat_name in names(repeat_categories)) {
    cat(sprintf("  %s: %d elements\n", cat_name, length(repeat_categories[[cat_name]])))
  }
  
  # 3. Run permutation tests with breakpoint considerations
  cat(sprintf("\nRunning permutation tests (%d iterations) with breakpoint handling...\n", n_permutations))
  
  # Define multiple evaluation functions for different stringency levels
  
  # Function 1: Standard overlap count (any overlap)
  count_any_overlap <- function(A, B) {
    return(numOverlaps(A, B))
  }
  
  # Function 2: Substantial overlap (>50% of repeat length)
  count_substantial_overlap <- function(A, B) {
    overlaps <- findOverlaps(A, B)
    if (length(overlaps) == 0) return(0)
    
    # Calculate overlap lengths
    overlap_ranges <- pintersect(A[queryHits(overlaps)], B[subjectHits(overlaps)])
    overlap_lengths <- width(overlap_ranges)
    repeat_lengths <- width(B[subjectHits(overlaps)])
    
    # Count overlaps where >50% of repeat is covered
    substantial <- overlap_lengths / repeat_lengths > 0.5
    return(sum(substantial))
  }
  
  # Function 3: Near-complete overlap (>80% of repeat length)
  count_complete_overlap <- function(A, B) {
    overlaps <- findOverlaps(A, B)
    if (length(overlaps) == 0) return(0)
    
    overlap_ranges <- pintersect(A[queryHits(overlaps)], B[subjectHits(overlaps)])
    overlap_lengths <- width(overlap_ranges)
    repeat_lengths <- width(B[subjectHits(overlaps)])
    
    # Count overlaps where >80% of repeat is covered
    complete <- overlap_lengths / repeat_lengths > 0.8
    return(sum(complete))
  }
  
  # Function 4: Overlap density (total bp of overlap)
  calculate_overlap_density <- function(A, B) {
    overlaps <- findOverlaps(A, B)
    if (length(overlaps) == 0) return(0)
    
    overlap_ranges <- pintersect(A[queryHits(overlaps)], B[subjectHits(overlaps)])
    total_overlap_bp <- sum(width(overlap_ranges))
    total_ecdna_bp <- sum(width(A))
    
    # Return density per Mb
    return((total_overlap_bp / total_ecdna_bp) * 1e6)
  }
  
  # Set up parallel processing
  if (n_cores > 1) {
    register(MulticoreParam(workers = n_cores))
  }
  
  # Initialize results storage for multiple analyses
  results_list <- list(
    any_overlap = list(),
    substantial_overlap = list(),
    complete_overlap = list(),
    overlap_density = list()
  )
  
  # Run permutation tests for each repeat category with multiple stringency levels
  for (category in names(repeat_categories)) {
    cat(sprintf("Testing %s repeats with multiple stringency levels...\n", category))
    
    # Test 1: Any overlap (standard)
    pt_any <- permTest(
      A = ecdna_gr,
      B = repeat_categories[[category]],
      ntimes = n_permutations,
      randomize.function = get(randomization_method),
      evaluate.function = count_any_overlap,
      genome = genome_gr,
      mc.cores = n_cores
    )
    results_list$any_overlap[[category]] <- pt_any
    
    # Test 2: Substantial overlap (>50%)
    pt_substantial <- permTest(
      A = ecdna_gr,
      B = repeat_categories[[category]],
      ntimes = n_permutations,
      randomize.function = get(randomization_method),
      evaluate.function = count_substantial_overlap,
      genome = genome_gr,
      mc.cores = n_cores
    )
    results_list$substantial_overlap[[category]] <- pt_substantial
    
    # Test 3: Complete overlap (>80%)
    pt_complete <- permTest(
      A = ecdna_gr,
      B = repeat_categories[[category]],
      ntimes = n_permutations,
      randomize.function = get(randomization_method),
      evaluate.function = count_complete_overlap,
      genome = genome_gr,
      mc.cores = n_cores
    )
    results_list$complete_overlap[[category]] <- pt_complete
    
    # Test 4: Overlap density
    pt_density <- permTest(
      A = ecdna_gr,
      B = repeat_categories[[category]],
      ntimes = n_permutations,
      randomize.function = get(randomization_method),
      evaluate.function = calculate_overlap_density,
      genome = genome_gr,
      mc.cores = n_cores
    )
    results_list$overlap_density[[category]] <- pt_density
  }
  
  # 4. Extract and summarize results for all stringency levels
  cat("\nExtracting results for all stringency levels...\n")
  
  # Create comprehensive results table
  comprehensive_results <- data.frame()
  
  stringency_levels <- c("any_overlap", "substantial_overlap", "complete_overlap", "overlap_density")
  stringency_labels <- c("Any Overlap", "Substantial (>50%)", "Complete (>80%)", "Density (bp/Mb)")
  
  for (i in seq_along(stringency_levels)) {
    stringency <- stringency_levels[i]
    label <- stringency_labels[i]
    
    for (category in names(repeat_categories)) {
      pt <- results_list[[stringency]][[category]]
      
      observed <- pt$observed
      permuted <- pt$permuted
      expected_mean <- mean(permuted)
      expected_sd <- sd(permuted)
      fold_change <- observed / expected_mean
      pvalue <- pt$pval
      zscore <- pt$zscore
      
      # Determine enrichment/depletion status
      status <- ifelse(observed > expected_mean, "enriched", "depleted")
      
      # Calculate confidence intervals
      ci <- quantile(permuted, c(0.025, 0.975))
      
      # Add to comprehensive results
      row_data <- data.frame(
        stringency_level = stringency,
        stringency_label = label,
        category = category,
        observed = observed,
        expected_mean = expected_mean,
        expected_sd = expected_sd,
        fold_change = fold_change,
        pvalue = pvalue,
        zscore = zscore,
        status = status,
        ci_lower = ci[1],
        ci_upper = ci[2],
        stringsAsFactors = FALSE
      )
      
      comprehensive_results <- rbind(comprehensive_results, row_data)
    }
  }
  
  # Apply multiple testing correction within each stringency level
  for (stringency in stringency_levels) {
    mask <- comprehensive_results$stringency_level == stringency
    pvals <- comprehensive_results$pvalue[mask]
    comprehensive_results$fdr[mask] <- p.adjust(pvals, method = "fdr")
    comprehensive_results$bonferroni[mask] <- p.adjust(pvals, method = "bonferroni")
    comprehensive_results$significant_fdr[mask] <- comprehensive_results$fdr[mask] < 0.05
  }
  
  # 5. Create enhanced visualizations comparing stringency levels
  cat("Creating comprehensive visualizations...\n")
  
  # Create comparison plots for different stringency levels
  plot_list <- list()
  
  # Plot 1: Fold change comparison across stringency levels
  p1 <- ggplot(comprehensive_results, aes(x = category, y = fold_change, fill = stringency_label)) +
    geom_col(position = "dodge", alpha = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    facet_wrap(~stringency_label, scales = "free_y") +
    labs(title = "Fold Change by Stringency Level", 
         subtitle = "Comparing different overlap thresholds",
         x = "Repeat Category", y = "Fold Change") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # Plot 2: P-value comparison
  p2 <- ggplot(comprehensive_results, aes(x = category, y = -log10(pvalue), fill = significant_fdr)) +
    geom_col(alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    facet_wrap(~stringency_label, scales = "free_y") +
    scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "darkred"), name = "FDR < 0.05") +
    labs(title = "Statistical Significance by Stringency Level", 
         x = "Repeat Category", y = "-log10(p-value)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Plot 3: Observed vs Expected for substantial overlap
  substantial_data <- comprehensive_results[comprehensive_results$stringency_level == "substantial_overlap", ]
  p3 <- ggplot(substantial_data, aes(x = category)) +
    geom_col(aes(y = observed), alpha = 0.7, fill = "red", width = 0.4, position = position_nudge(x = -0.2)) +
    geom_col(aes(y = expected_mean), alpha = 0.7, fill = "blue", width = 0.4, position = position_nudge(x = 0.2)) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1, position = position_nudge(x = 0.2)) +
    labs(title = "Observed vs Expected (Substantial Overlap >50%)", 
         subtitle = "Red: Observed, Blue: Expected from permutations",
         x = "Repeat Category", y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Plot 4: Stringency effect - how results change with stricter thresholds
  stringency_effect <- comprehensive_results %>%
    group_by(category, stringency_label) %>%
    summarise(fold_change = first(fold_change), .groups = "drop")
  
  p4 <- ggplot(stringency_effect, aes(x = stringency_label, y = fold_change, group = category, color = category)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
    labs(title = "Effect of Stringency on Fold Change", 
         subtitle = "How enrichment changes with overlap thresholds",
         x = "Stringency Level", y = "Fold Change", color = "Repeat Category") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  plot_list <- list(fold_change_comparison = p1, pvalue_comparison = p2, 
                   substantial_overlap = p3, stringency_effect = p4)
  
  # Combine plots
  combined_plots <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
  
  # 6. Print summary report
  cat("\n=== RESULTS SUMMARY ===\n")
  cat(sprintf("Analysis completed with %d permutations\n", n_permutations))
  cat(sprintf("Total ecDNA regions: %d\n", length(ecdna_gr)))
  cat(sprintf("Total ecDNA length: %s bp\n", format(sum(width(ecdna_gr)), big.mark = ",")))
  
  # 7. Return comprehensive results with breakpoint analysis
  return(list(
    # Main results table with all stringency levels
    comprehensive_results = comprehensive_results,
    
    # Legacy format for any_overlap (backward compatibility)
    results_table = comprehensive_results[comprehensive_results$stringency_level == "any_overlap", ],
    
    # All permutation test objects organized by stringency
    permtest_objects = results_list,
    
    # Enhanced plots
    plots = plot_list,
    combined_plots = combined_plots,
    
    # Input data
    input_data = list(
      ecdna_regions = ecdna_gr,
      repeat_elements = repeat_gr,
      repeat_categories = repeat_categories
    ),
    
    # Analysis parameters
    parameters = list(
      n_permutations = n_permutations,
      genome = genome,
      randomization_method = randomization_method,
      stringency_levels = stringency_labels
    ),
    
    # Breakpoint analysis summary
    breakpoint_analysis = list(
      description = "Multiple stringency levels to handle breakpoint artifacts",
      stringency_levels = data.frame(
        level = stringency_levels,
        description = c(
          "Any overlap (standard, may include truncated elements)",
          "Substantial overlap (>50% of repeat covered)",
          "Complete overlap (>80% of repeat covered)", 
          "Overlap density (total bp overlap per Mb ecDNA)"
        ),
        recommendation = c(
          "Use for general screening",
          "Recommended for main analysis",
          "Conservative estimate",
          "Best for quantitative analysis"
        )
      )
    )
  ))
}

#' Extract detailed statistics from permutation test results
get_detailed_stats <- function(analysis_results) {
  detailed_stats <- list()
  
  for (category in names(analysis_results$permtest_objects$any_overlap)) {
    pt <- analysis_results$permtest_objects$any_overlap[[category]]
    
    detailed_stats[[category]] <- list(
      observed = pt$observed,
      permuted_values = pt$permuted,
      mean_random = mean(pt$permuted),
      sd_random = sd(pt$permuted),
      min_random = min(pt$permuted),
      max_random = max(pt$permuted),
      percentiles = quantile(pt$permuted, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)),
      zscore = pt$zscore,
      pvalue = pt$pval,
      alternative = pt$alternative
    )
  }
  
  return(detailed_stats)
}
