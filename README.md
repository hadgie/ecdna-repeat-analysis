# ecDNA Repeat Element Enrichment Analysis - Docker Implementation

A containerized solution for analyzing repeat element enrichment in extrachromosomal DNA (ecDNA) regions using regioneR permutation testing.

## 🚀 Quick Start

### 1. Build the Docker Image

```bash
# Clone/download the repository
git clone <your-repo-url>
cd ecdna-repeat-analysis

# Build the Docker image
docker build -t ecdna-analysis:latest .
```

### 2. Prepare Your Data

Create a data directory with your input files:

```bash
mkdir -p data output

# Copy your files
cp /path/to/your/ecdna_regions.bed data/
cp /path/to/your/repeatmasker.bed data/
```

### 3. Run the Analysis

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/output:/app/output \
  ecdna-analysis:latest \
  --ecdna /data/ecdna_regions.bed \
  --repeats /data/repeatmasker.bed \
  --genome hg38 \
  --permutations 10000 \
  --cores 4
```

## 📋 Requirements

- **Docker** (version 20.0+ recommended)
- **Memory**: At least 4GB RAM (8GB+ recommended for large datasets)
- **CPU**: Multi-core processor (4+ cores recommended)
- **Storage**: ~10GB for Docker image + space for your data and results

## 📁 Input File Formats

### ecDNA Regions BED File
6-column BED format with ecDNA coordinates:
```
chr1	1000000	1500000	ecDNA_1	1000	+
chr2	2000000	2800000	ecDNA_2	800	-
chr3	3000000	3200000	ecDNA_3	200	+
```

### RepeatMasker BED File
6-column format with repeat family names:
```
chr1	16777160	16777470	AluSp	2147	+
chr1	25165800	25166089	AluY	2626	-
chr1	33553606	33554646	L2b	626	+
chr1	50330063	50332153	L1PA10	12545	+
```

## 🔧 Command Line Options

### Required Arguments
- `--ecdna FILE`: Path to ecDNA regions BED file
- `--repeats FILE`: Path to RepeatMasker BED file

### Optional Arguments
- `--genome STRING`: Genome version (`hg38` or `hg19`) [default: hg38]
- `--permutations INTEGER`: Number of permutations [default: 10000]
- `--cores INTEGER`: Number of CPU cores [default: 4]
- `--output PATH`: Output directory [default: /app/output]
- `--prefix STRING`: Output file prefix [default: ecdna_analysis]
- `--stringency LEVEL`: Analysis level (`all`, `any`, `substantial`, `complete`, `density`) [default: all]
- `--method STRING`: Randomization method [default: circularRandomizeRegions]
- `--save-plots`: Save visualization plots [default: TRUE]
- `--save-data`: Save R data objects [default: TRUE]
- `--verbose`: Enable verbose output

## 📊 Output Files

The analysis generates several output files:

### Results Tables
- `*_comprehensive_results.csv`: All stringency levels combined
- `*_substantial_overlap_results.csv`: Recommended main results (>50% overlap)
- `*_any_overlap_results.csv`: Standard overlap results
- `*_complete_overlap_results.csv`: Conservative results (>80% overlap)
- `*_overlap_density_results.csv`: Quantitative density results

### Visualizations
- `*_fold_change_comparison.png`: Fold change across stringency levels
- `*_pvalue_comparison.png`: Statistical significance comparison
- `*_substantial_overlap.png`: Main results visualization
- `*_stringency_effect.png`: Effect of different thresholds
- `*_combined_plots.png`: All plots combined

### Reports
- `*_summary_report.txt`: Human-readable analysis summary
- `*_analysis_data.RData`: Complete R analysis object

## 🐳 Docker Usage Examples

### Basic Analysis
```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/output:/app/output \
  ecdna-analysis:latest \
  --ecdna /data/my_ecdna.bed \
  --repeats /data/my_repeats.bed \
  --genome hg38
```

### Quick Test (Fewer Permutations)
```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/output:/app/output \
  ecdna-analysis:latest \
  --ecdna /data/my_ecdna.bed \
  --repeats /data/my_repeats.bed \
  --genome hg19 \
  --permutations 1000 \
  --cores 2 \
  --stringency substantial
```

### Comprehensive Analysis
```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/output \
  ecdna-analysis:latest \
  --ecdna /data/ecdna_regions.bed \
  --repeats /data/repeatmasker.bed \
  --genome hg38 \
  --permutations 10000 \
  --cores 8 \
  --prefix "comprehensive_analysis" \
  --verbose
```

### Using Your Actual Files
```bash
# For your SNU620 analysis
docker run --rm \
  -v /Users/jeonghuiseog/Library/CloudStorage/Dropbox/Mac/Downloads:/data \
  -v $(pwd)/results:/app/output \
  ecdna-analysis:latest \
  --ecdna /data/SNU620_roi.bed \
  --repeats /data/repeatmasker_hg19.bed \
  --genome hg19 \
  --permutations 10000 \
  --cores 6 \
  --prefix "SNU620_analysis" \
  --verbose
```

## 🔄 Using Docker Compose

For easier management, use Docker Compose:

```bash
# Start the analysis
docker-compose run --rm ecdna-analysis \
  --ecdna /data/ecdna_regions.bed \
  --repeats /data/repeatmasker.bed \
  --genome hg38 \
  --permutations 10000

# Or use the development version
docker-compose run --rm ecdna-analysis-dev \
  --ecdna /data/ecdna_regions.bed \
  --repeats /data/repeatmasker.bed \
  --verbose
```

## 📈 Interpreting Results

### Stringency Levels
The analysis provides multiple stringency levels to handle breakpoint artifacts:

1. **Any Overlap**: Standard approach, counts any overlap ≥1bp
2. **Substantial Overlap** (Recommended): Only counts repeats with >50% overlap
3. **Complete Overlap**: Conservative, only counts repeats with >80% overlap  
4. **Overlap Density**: Quantitative measure (bp overlap per Mb ecDNA)

### Key Metrics
- **Fold Change**: Observed/Expected ratio (>1 = enriched, <1 = depleted)
- **P-value**: Statistical significance from permutation test
- **FDR**: False Discovery Rate corrected p-value
- **Status**: "enriched" or "depleted"

### Recommendations
- Focus on **Substantial Overlap** results for main conclusions
- Elements significant across multiple stringency levels are most robust
- Consider both statistical significance (FDR < 0.05) and biological relevance (fold change)

## 🛠️ Troubleshooting

### Common Issues

**Docker Build Fails**
```bash
# Check Docker is running
docker info

# Build with no cache if needed
docker build --no-cache -t ecdna-analysis:latest .
```

**Out of Memory**
```bash
# Reduce permutations or cores
docker run --memory=4g --cpus=2 ecdna-analysis:latest [options]

# Or increase Docker memory limits in Docker Desktop
```

**File Not Found**
```bash
# Check file paths and permissions
ls -la data/
docker run --rm -v $(pwd)/data:/data alpine ls -la /data
```

**Analysis Takes Too Long**
```bash
# Start with fewer permutations for testing
--permutations 1000

# Use fewer cores if system is overloaded
--cores 2
```

### Performance Tips

- **Memory**: Allocate 8GB+ RAM for large datasets
- **Cores**: Use 4-8 cores for optimal performance
- **Permutations**: Start with 1000 for testing, use 10000+ for final analysis
- **Storage**: Ensure sufficient disk space for results

## 🧬 Biological Interpretation

### Repeat Element Categories
- **SINE**: Short Interspersed Nuclear Elements (mainly Alu)
- **LINE**: Long Interspersed Nuclear Elements (mainly L1)
- **LTR**: Long Terminal Repeat retrotransposons
- **DNA**: DNA transposons
- **Simple**: Simple sequence repeats
- **Satellite**: Satellite DNA

### Expected Results
- ecDNA regions often show enrichment for certain repeat classes
- L1 elements frequently enriched due to recombination hotspots
- Alu elements may be depleted in some contexts
- Results vary by cancer type and ecDNA formation mechanism

## 📚 References

- **regioneR**: Gel et al. (2016) Bioinformatics
- **RepeatMasker**: Smit, AFA, Hubley, R & Green, P
- **ecDNA formation**: Multiple recent studies in Nature, Science, Cell

## 🤝 Support

For issues or questions:
1. Check the troubleshooting section above
2. Verify input file formats
3. Test with smaller datasets first
4. Review Docker logs: `docker logs <container_id>`

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.
