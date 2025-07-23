#!/usr/bin/env Rscript

# Command-line interface for ecDNA Repeat Element Enrichment Analysis
# Author: Your Name
# Date: 2025

# Load required libraries
suppressMessages({
  library(optparse)
  library(regioneR)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(readr)
})

# Source the main analysis functions
source("/app/ecdna_analysis.R")

# Define command line options
option_list <- list(
  make_option(c("-e", "--ecdna"), type="character", default=NULL,
              help="Path to ecDNA regions BED file [required]", metavar="FILE"),
  
  make_option(c("-r", "--repeats"), type="character", default=NULL,
              help="Path to RepeatMasker BED file [required]", metavar="FILE"),
  
  make_option(c("-g", "--genome"), type="character", default="hg38",
              help="Genome version (hg38 or hg19) [default: %default]", metavar="STRING"),
  
  make_option(c("-n", "--permutations"), type="integer", default=10000,
              help="Number of permutations [default: %default]", metavar="INTEGER"),
  
  make_option(c("-c", "--cores"), type="integer", default=4,
              help="Number of CPU cores to use [default: %default]", metavar="INTEGER"),
  
  make_option(c("-o", "--output"), type="character", default="/app/output",
              help="Output directory [default: %default]", metavar="PATH"),
  
  make_option(c("-p", "--prefix"), type="character", default="ecdna_analysis",
              help="Output file prefix [default: %default]", metavar="STRING"),
  
  make_option(c("-m", "--method"), type="character", default="circularRandomizeRegions",
              help="Randomization method [default: %default]", metavar="STRING"),
  
  make_option(c("-s", "--stringency"), type="character", default="all",
              help="Stringency level (all, any, substantial, complete, density) [default: %default]", metavar="STRING"),
  
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Enable verbose output"),
  
  make_option(c("--save-plots"), action="store_true", default=TRUE,
              help="Save visualization plots [default: TRUE]"),
  
  make_option(c("--save-data"), action="store_true", default=TRUE,
              help="Save intermediate data objects [default: TRUE]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list,
                          description="ecDNA Repeat Element Enrichment Analysis using regioneR",
                          epilogue="Example usage:\n  docker run ecdna-analysis --ecdna /data/ecdna.bed --repeats /data/repeats.bed --genome hg38 --permutations 10000")

opt <- parse_args(opt_parser)

# Validation function
validate_inputs <- function(opt) {
  errors <- c()
  
  # Check required arguments
  if (is.null(opt$ecdna)) {
    errors <- c(errors, "ERROR: --ecdna argument is required")
  } else if (!file.exists(opt$ecdna)) {
    errors <- c(errors, paste("ERROR: ecDNA file does not exist:", opt$ecdna))
  }
  
  if (is.null(opt$repeats)) {
    errors <- c(errors, "ERROR: --repeats argument is required")
  } else if (!file.exists(opt$repeats)) {
    errors <- c(errors, paste("ERROR: RepeatMasker file does not exist:", opt$repeats))
  }
  
  # Check genome version
  if (!opt$genome %in% c("hg38", "hg19")) {
    errors <- c(errors, "ERROR: --genome must be 'hg38' or 'hg19'")
  }
  
  # Check numeric arguments
  if (opt$permutations < 100) {
    errors <- c(errors, "ERROR: --permutations must be at least 100")
  }
  
  if (opt$cores < 1) {
    errors <- c(errors, "ERROR: --cores must be at least 1")
  }
  
  # Check stringency option
  valid_stringency <- c("all", "any", "substantial", "complete", "density")
  if (!opt$stringency %in% valid_stringency) {
    errors <- c(errors, paste("ERROR: --stringency must be one of:", paste(valid_stringency, collapse=", ")))
  }
  
  # Check randomization method
  valid_methods <- c("randomizeRegions", "circularRandomizeRegions", "resampleRegions")
  if (!opt$method %in% valid_methods) {
    errors <- c(errors, paste("ERROR: --method must be one of:", paste(valid_methods, collapse=", ")))
  }
  
  return(errors)
}

# Helper function to print status messages
print_status <- function(message, verbose=TRUE) {
  if (verbose || opt$verbose) {
    cat(sprintf("[%s] %s\n", Sys.time(), message))
  }
}

# Main analysis function
run_ecdna_analysis <- function(opt) {
  
  # Create output directory
  if (!dir.exists(opt$output)) {
    dir.create(opt$output, recursive = TRUE)
  }
  
  print_status("Starting ecDNA Repeat Element Enrichment Analysis...")
  print_status(sprintf("ecDNA file: %s", opt$ecdna))
  print_status(sprintf("RepeatMasker file: %s", opt$repeats))
  print_status(sprintf("Genome: %s", opt$genome))
  print_status(sprintf("Permutations: %d", opt$permutations))
  print_status(sprintf("CPU cores: %d", opt$cores))
  print_status(sprintf("Output directory: %s", opt$output))
  
  # Run the analysis
  tryCatch({
    results <- analyze_ecdna_repeats(
      ecdna_bed = opt$ecdna,
      repeat_bed = opt$repeats,
      genome = opt$genome,
      n_permutations = opt$permutations,
      n_cores = opt$cores,
      randomization_method = opt$method
    )
    
    print_status("Analysis completed successfully!")
    
    # Save results based on stringency selection
    if (opt$stringency == "all") {
      # Save comprehensive results
      results_file <- file.path(opt$output, paste0(opt$prefix, "_comprehensive_results.csv"))
      write_csv(results$comprehensive_results, results_file)
      print_status(sprintf("Comprehensive results saved to: %s", results_file))
      
      # Save summary for each stringency level
      for (level in c("any_overlap", "substantial_overlap", "complete_overlap", "overlap_density")) {
        level_data <- results$comprehensive_results[results$comprehensive_results$stringency_level == level, ]
        level_file <- file.path(opt$output, paste0(opt$prefix, "_", level, "_results.csv"))
        write_csv(level_data, level_file)
      }
      
    } else {
      # Save specific stringency level
      stringency_map <- list(
        "any" = "any_overlap",
        "substantial" = "substantial_overlap", 
        "complete" = "complete_overlap",
        "density" = "overlap_density"
      )
      
      selected_level <- stringency_map[[opt$stringency]]
      level_data <- results$comprehensive_results[results$comprehensive_results$stringency_level == selected_level, ]
      results_file <- file.path(opt$output, paste0(opt$prefix, "_", opt$stringency, "_results.csv"))
      write_csv(level_data, results_file)
      print_status(sprintf("Results (%s stringency) saved to: %s", opt$stringency, results_file))
    }
    
    # Save plots
    if (opt$save_plots) {
      print_status("Generating and saving plots...")
      
      # Save individual plots
      plot_files <- list()
      
      if (opt$stringency == "all") {
        # Save all plots
        for (plot_name in names(results$plots)) {
          plot_file <- file.path(opt$output, paste0(opt$prefix, "_", plot_name, ".png"))
          ggsave(plot_file, results$plots[[plot_name]], width = 12, height = 8, dpi = 300)
          plot_files[[plot_name]] <- plot_file
        }
        
        # Save combined plot
        combined_file <- file.path(opt$output, paste0(opt$prefix, "_combined_plots.png"))
        ggsave(combined_file, results$combined_plots, width = 16, height = 12, dpi = 300)
        plot_files[["combined"]] <- combined_file
        
      } else {
        # Save substantial overlap plot (most informative single plot)
        main_plot_file <- file.path(opt$output, paste0(opt$prefix, "_substantial_overlap.png"))
        ggsave(main_plot_file, results$plots$substantial_overlap, width = 12, height = 8, dpi = 300)
        plot_files[["main"]] <- main_plot_file
      }
      
      print_status("Plots saved successfully!")
      for (name in names(plot_files)) {
        print_status(sprintf("  %s: %s", name, plot_files[[name]]))
      }
    }
    
    # Save data objects
    if (opt$save_data) {
      print_status("Saving analysis data objects...")
      data_file <- file.path(opt$output, paste0(opt$prefix, "_analysis_data.RData"))
      save(results, file = data_file)
      print_status(sprintf("Analysis data saved to: %s", data_file))
    }
    
    # Generate summary report
    print_status("Generating summary report...")
    generate_summary_report(results, opt)
    
    print_status("Analysis pipeline completed successfully!")
    
  }, error = function(e) {
    cat(sprintf("ERROR: Analysis failed with error: %s\n", e$message))
    quit(status = 1)
  })
}

# Function to generate summary report
generate_summary_report <- function(results, opt) {
  report_file <- file.path(opt$output, paste0(opt$prefix, "_summary_report.txt"))
  
  # Get substantial overlap results (most informative)
  substantial_results <- results$comprehensive_results[
    results$comprehensive_results$stringency_level == "substantial_overlap", ]
  
  # Generate report
  report_lines <- c(
    "=" %+% paste(rep("=", 60), collapse="") %+% "=",
    "ecDNA REPEAT ELEMENT ENRICHMENT ANALYSIS SUMMARY",
    "=" %+% paste(rep("=", 60), collapse="") %+% "=",
    "",
    sprintf("Analysis Date: %s", Sys.time()),
    sprintf("ecDNA file: %s", basename(opt$ecdna)),
    sprintf("RepeatMasker file: %s", basename(opt$repeats)),
    sprintf("Genome version: %s", opt$genome),
    sprintf("Permutations: %d", opt$permutations),
    sprintf("CPU cores used: %d", opt$cores),
    sprintf("Total ecDNA regions: %d", length(results$input_data$ecdna_regions)),
    sprintf("Total ecDNA length: %s bp", format(sum(width(results$input_data$ecdna_regions)), big.mark=",")),
    "",
    "RESULTS (Substantial Overlap >50% threshold):",
    paste(rep("-", 50), collapse="")
  )
  
  # Add results for each category
  for (i in 1:nrow(substantial_results)) {
    row <- substantial_results[i, ]
    report_lines <- c(report_lines,
      "",
      sprintf("%s Repeats:", row$category),
      sprintf("  Observed overlaps: %d", row$observed),
      sprintf("  Expected overlaps: %.1f ± %.1f", row$expected_mean, row$expected_sd),
      sprintf("  Fold change: %.2f", row$fold_change),
      sprintf("  P-value: %.4f", row$pvalue),
      sprintf("  FDR-corrected: %.4f", row$fdr),
      sprintf("  Status: %s", row$status),
      sprintf("  Significant (FDR < 0.05): %s", ifelse(row$significant_fdr, "YES", "NO"))
    )
  }
  
  # Add interpretation
  significant_enriched <- substantial_results[substantial_results$significant_fdr & substantial_results$status == "enriched", ]
  significant_depleted <- substantial_results[substantial_results$significant_fdr & substantial_results$status == "depleted", ]
  
  report_lines <- c(report_lines,
    "",
    "INTERPRETATION:",
    paste(rep("-", 20), collapse="")
  )
  
  if (nrow(significant_enriched) > 0) {
    report_lines <- c(report_lines,
      sprintf("Significantly ENRICHED repeat types: %s", paste(significant_enriched$category, collapse=", "))
    )
  }
  
  if (nrow(significant_depleted) > 0) {
    report_lines <- c(report_lines,
      sprintf("Significantly DEPLETED repeat types: %s", paste(significant_depleted$category, collapse=", "))
    )
  }
  
  if (nrow(significant_enriched) == 0 && nrow(significant_depleted) == 0) {
    report_lines <- c(report_lines,
      "No significant enrichment or depletion detected for any repeat category."
    )
  }
  
  report_lines <- c(report_lines,
    "",
    "=" %+% paste(rep("=", 60), collapse="") %+% "="
  )
  
  # Write report
  writeLines(report_lines, report_file)
  print_status(sprintf("Summary report saved to: %s", report_file))
  
  # Also print to console
  if (opt$verbose) {
    cat("\n")
    cat(paste(report_lines, collapse="\n"))
    cat("\n")
  }
}

# String concatenation operator
`%+%` <- function(x, y) paste0(x, y)

# Main execution
main <- function() {
  # Validate inputs
  errors <- validate_inputs(opt)
  
  if (length(errors) > 0) {
    cat("Input validation failed:\n")
    for (error in errors) {
      cat(paste("  ", error, "\n"))
    }
    cat("\nUse --help for usage information.\n")
    quit(status = 1)
  }
  
  # Run analysis
  run_ecdna_analysis(opt)
}

# Execute main function
if (!interactive()) {
  main()
}
