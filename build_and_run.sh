#!/bin/bash
# Docker build and run scripts for ecDNA analysis

# =============================================================================
# BUILD SCRIPT (build_docker.sh)
# =============================================================================

build_docker() {
    echo "Building ecDNA repeat analysis Docker image..."
    
    # Build the Docker image
    docker build -t ecdna-analysis:latest .
    
    if [ $? -eq 0 ]; then
        echo "✅ Docker image built successfully!"
        echo "Image name: ecdna-analysis:latest"
        echo ""
        echo "To run the analysis, use:"
        echo "docker run --rm -v /path/to/your/data:/data -v /path/to/output:/app/output ecdna-analysis:latest [OPTIONS]"
    else
        echo "❌ Docker build failed!"
        exit 1
    fi
}

# =============================================================================
# RUN SCRIPT EXAMPLES
# =============================================================================

# Example 1: Basic usage
run_basic_example() {
    echo "Running basic ecDNA analysis example..."
    
    docker run --rm \
        -v "$(pwd)/test_data:/data" \
        -v "$(pwd)/output:/app/output" \
        ecdna-analysis:latest \
        --ecdna /data/ecdna_regions.bed \
        --repeats /data/repeatmasker.bed \
        --genome hg38 \
        --permutations 10000 \
        --cores 4 \
        --prefix "my_analysis"
}

# Example 2: Quick test with fewer permutations
run_test_example() {
    echo "Running quick test analysis..."
    
    docker run --rm \
        -v "$(pwd)/test_data:/data" \
        -v "$(pwd)/output:/app/output" \
        ecdna-analysis:latest \
        --ecdna /data/ecdna_regions.bed \
        --repeats /data/repeatmasker.bed \
        --genome hg19 \
        --permutations 1000 \
        --cores 2 \
        --stringency substantial \
        --verbose
}

# Example 3: Full analysis with all options
run_full_example() {
    echo "Running comprehensive ecDNA analysis..."
    
    docker run --rm \
        -v "/Users/jeonghuiseog/data:/data" \
        -v "/Users/jeonghuiseog/results:/app/output" \
        ecdna-analysis:latest \
        --ecdna /data/SNU620_roi.bed \
        --repeats /data/repeatmasker_hg19.bed \
        --genome hg19 \
        --permutations 10000 \
        --cores 6 \
        --output /app/output \
        --prefix "SNU620_comprehensive" \
        --method circularRandomizeRegions \
        --stringency all \
        --save-plots \
        --save-data \
        --verbose
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Show help for Docker usage
show_docker_help() {
    echo "ecDNA Repeat Analysis Docker Container"
    echo "======================================"
    echo ""
    echo "BUILD:"
    echo "  docker build -t ecdna-analysis:latest ."
    echo ""
    echo "RUN:"
    echo "  docker run --rm -v /host/data:/data -v /host/output:/app/output ecdna-analysis:latest [OPTIONS]"
    echo ""
    echo "REQUIRED OPTIONS:"
    echo "  --ecdna /data/ecdna.bed       ecDNA regions BED file"
    echo "  --repeats /data/repeats.bed   RepeatMasker BED file"
    echo ""
    echo "OPTIONAL OPTIONS:"
    echo "  --genome hg38|hg19           Genome version [default: hg38]"
    echo "  --permutations N             Number of permutations [default: 10000]"
    echo "  --cores N                    CPU cores [default: 4]"
    echo "  --output /app/output         Output directory [default: /app/output]"
    echo "  --prefix NAME                Output file prefix [default: ecdna_analysis]"
    echo "  --stringency LEVEL           all|any|substantial|complete|density [default: all]"
    echo "  --method METHOD              Randomization method [default: circularRandomizeRegions]"
    echo "  --save-plots                 Save visualization plots [default: TRUE]"
    echo "  --save-data                  Save R data objects [default: TRUE]"
    echo "  --verbose                    Enable verbose output"
    echo ""
    echo "EXAMPLES:"
    echo "  # Basic usage"
    echo "  docker run --rm -v \$(pwd)/data:/data -v \$(pwd)/output:/app/output \\"
    echo "    ecdna-analysis:latest --ecdna /data/ecdna.bed --repeats /data/repeats.bed"
    echo ""
    echo "  # Full analysis with custom settings"
    echo "  docker run --rm -v /path/to/data:/data -v /path/to/results:/app/output \\"
    echo "    ecdna-analysis:latest \\"
    echo "    --ecdna /data/ecdna_regions.bed \\"
    echo "    --repeats /data/repeatmasker.bed \\"
    echo "    --genome hg38 \\"
    echo "    --permutations 10000 \\"
    echo "    --cores 8 \\"
    echo "    --prefix comprehensive_analysis \\"
    echo "    --verbose"
}

# Validate data files exist
validate_data_files() {
    local data_dir="$1"
    
    echo "Validating data files in: $data_dir"
    
    if [ ! -d "$data_dir" ]; then
        echo "❌ Data directory does not exist: $data_dir"
        return 1
    fi
    
    # Check for common file patterns
    ecdna_files=$(find "$data_dir" -name "*ecdna*" -o -name "*roi*" -o -name "*amplicon*" | head -5)
    repeat_files=$(find "$data_dir" -name "*repeat*" -o -name "*rmsk*" | head -5)
    
    echo "Found potential ecDNA files:"
    echo "$ecdna_files"
    echo ""
    echo "Found potential repeat files:"
    echo "$repeat_files"
    echo ""
    
    return 0
}

# Check Docker installation
check_docker() {
    if ! command -v docker &> /dev/null; then
        echo "❌ Docker is not installed or not in PATH"
        echo "Please install Docker: https://docs.docker.com/get-docker/"
        return 1
    fi
    
    if ! docker info &> /dev/null; then
        echo "❌ Docker daemon is not running"
        echo "Please start Docker daemon"
        return 1
    fi
    
    echo "✅ Docker is available"
    return 0
}

# Create example data directory structure
create_example_structure() {
    local base_dir="$1"
    
    echo "Creating example directory structure in: $base_dir"
    
    mkdir -p "$base_dir"/{data,output,scripts}
    
    cat > "$base_dir/data/README.md" << 'EOF'
# Data Directory

Place your input files here:

## Required Files:
- `ecdna_regions.bed`: BED file with ecDNA coordinates
- `repeatmasker.bed`: RepeatMasker annotations in BED format

## File Formats:

### ecDNA regions BED file:
```
chr1    1000000    1500000    ecDNA_1    1000    +
chr2    2000000    2800000    ecDNA_2    800     -
```

### RepeatMasker BED file (6 columns):
```
chr1    16777160    16777470    AluSp    2147    +
chr1    25165800    25166089    AluY     2626    -
chr1    33553606    33554646    L2b      626     +
```
EOF

    cat > "$base_dir/run_analysis.sh" << 'EOF'
#!/bin/bash
# Example analysis script

# Set your data paths
DATA_DIR="$(pwd)/data"
OUTPUT_DIR="$(pwd)/output"

# Run the analysis
docker run --rm \
    -v "$DATA_DIR:/data" \
    -v "$OUTPUT_DIR:/app/output" \
    ecdna-analysis:latest \
    --ecdna /data/ecdna_regions.bed \
    --repeats /data/repeatmasker.bed \
    --genome hg38 \
    --permutations 10000 \
    --cores 4 \
    --prefix "my_analysis" \
    --verbose
EOF

    chmod +x "$base_dir/run_analysis.sh"
    
    echo "✅ Example structure created!"
    echo "Edit $base_dir/run_analysis.sh with your specific file paths"
}

# =============================================================================
# MAIN SCRIPT LOGIC
# =============================================================================

# Parse command line arguments
case "${1:-help}" in
    "build")
        check_docker && build_docker
        ;;
    "test")
        check_docker && run_test_example
        ;;
    "basic")
        check_docker && run_basic_example
        ;;
    "full")
        check_docker && run_full_example
        ;;
    "validate")
        validate_data_files "${2:-./data}"
        ;;
    "setup")
        create_example_structure "${2:-./ecdna_analysis_project}"
        ;;
    "help"|*)
        show_docker_help
        ;;
esac
