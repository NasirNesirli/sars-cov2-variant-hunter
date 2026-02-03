# Usage Guide

## Getting Started

### 1. Installation

Ensure you have [Pixi](https://pixi.sh) installed on your macOS system:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

Clone and setup the project:

```bash
git clone https://github.com/yourusername/sars-cov2-variant-hunter.git
cd sars-cov2-variant-hunter
pixi install
```

### 2. Data Preparation

#### Download Reference Genome

```bash
mkdir -p data/reference
curl -o data/reference/MN908947.3.fasta \
  https://raw.githubusercontent.com/jonas-fuchs/varVAMP_in_silico_analysis/11bf299d73cc069fbc8fa2be86ed57ada649afcd/reference_seq/SARS-CoV-2/MN908947.3.fasta
```

#### Prepare Your Sequencing Data

Place your paired-end FASTQ files in `data/raw_fastq/`:

```bash
mkdir -p data/raw_fastq
# Copy your files:
# - Sample_R1_001.fastq.gz (forward reads)
# - Sample_R2_001.fastq.gz (reverse reads)
```

**Example data source**: [SARS-CoV-2 FASTQ samples](https://figshare.com/articles/dataset/SARS-CoV-2_fastq_samples/24556747)

### 3. Configure Pipeline

Edit [pipeline.sh](../pipeline.sh) to match your FASTQ filenames:

```bash
# Change these lines:
bwa mem data/reference/MN908947.3.fasta \
    data/raw_fastq/YOUR_SAMPLE_R1_001.fastq.gz \
    data/raw_fastq/YOUR_SAMPLE_R2_001.fastq.gz \
    | samtools view -b -o results/aligned/aligned_reads.bam
```

### 4. Run Pipeline

Execute the complete workflow:

```bash
pixi run bash pipeline.sh
```

Expected runtime: 5-15 minutes (depending on file size and system)

## Advanced Usage

### Check Data Quality

Before running the pipeline:

```bash
# View FASTQ statistics
pixi run seqkit stats data/raw_fastq/*.fastq.gz

# View reference genome info
pixi run seqkit stats data/reference/MN908947.3.fasta
```

### Inspect Results

#### View Alignment Statistics

```bash
# Overall alignment stats
pixi run samtools flagstat results/aligned/sorted_reads.bam

# Coverage statistics
pixi run samtools coverage results/aligned/sorted_reads.bam

# Depth per position
pixi run samtools depth results/aligned/sorted_reads.bam | head -20
```

#### Analyze Variants

```bash
# Count variants
pixi run grep -v "^#" results/vcf/final_variants.vcf | wc -l

# View variant statistics
pixi run bcftools stats results/vcf/final_variants.vcf

# Filter specific positions
pixi run bcftools view -r MN908947.3:14408-14410 results/vcf/final_variants.vcf

# Extract high-confidence variants (QUAL>100, DP>100)
pixi run bcftools view -i 'QUAL>100 && DP>100' results/vcf/final_variants.vcf
```

### Custom Filtering

Adjust quality thresholds in [pipeline.sh](../pipeline.sh):

```bash
# Example: More stringent filtering
bcftools view -i 'QUAL>50 && DP>50' results/vcf/raw_variants.vcf > results/vcf/stringent_variants.vcf
```

## Troubleshooting

### Issue: Pipeline fails at alignment step

**Solution**: Ensure FASTQ filenames match exactly in pipeline.sh

```bash
# List your actual files
ls -lh data/raw_fastq/
# Update pipeline.sh accordingly
```

### Issue: Low variant count

**Possible causes**:
- Low sequencing coverage
- Sample very similar to reference
- Quality thresholds too stringent

**Check**:
```bash
# Check alignment rate
pixi run samtools flagstat results/aligned/sorted_reads.bam

# View raw variants before filtering
pixi run grep -v "^#" results/vcf/raw_variants.vcf | wc -l
```

### Issue: Reference index missing

**Solution**: The pipeline creates indices automatically, but if needed:

```bash
pixi run bwa index data/reference/MN908947.3.fasta
pixi run samtools faidx data/reference/MN908947.3.fasta
```

## Output Files

| File | Description | Size (typical) |
|------|-------------|----------------|
| `results/aligned/aligned_reads.bam` | Raw alignment | 50-500 MB |
| `results/aligned/sorted_reads.bam` | Sorted alignment | 50-500 MB |
| `results/aligned/sorted_reads.bam.bai` | BAM index | 1-10 MB |
| `results/vcf/raw_variants.vcf` | All variants | 10-100 KB |
| `results/vcf/final_variants.vcf` | Filtered variants | 5-50 KB |

## Performance Tips

1. **Use SSD storage** for faster I/O
2. **Increase memory** if processing multiple samples
3. **Parallel processing**: For multiple samples, run pipelines in parallel

```bash
# Example: Process multiple samples
for sample in sample1 sample2 sample3; do
    bash pipeline.sh ${sample} &
done
wait
```

## Next Steps After Variant Calling

1. **Lineage Assignment**: Use [Pangolin](https://pangolin.cog-uk.io/) or [NextClade](https://clades.nextstrain.org/)
2. **Variant Annotation**: Use [SnpEff](https://pcingola.github.io/SnpEff/) or [VEP](https://www.ensembl.org/vep)
3. **Phylogenetic Analysis**: Compare with global sequences
4. **Report Generation**: Create detailed analysis reports

## Additional Resources

- [BWA Manual](http://bio-bwa.sourceforge.net/bwa.shtml)
- [SAMtools Manual](http://www.htslib.org/doc/samtools.html)
- [BCFtools Manual](http://www.htslib.org/doc/bcftools.html)
- [VCF Format Specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)

## Support

For issues or questions:
- Open an issue on [GitHub](https://github.com/yourusername/sars-cov2-variant-hunter/issues)
- Email: nasir.nesirli@gmail.com
