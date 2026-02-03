# Pipeline Output Examples

This document shows typical outputs from the SARS-CoV-2 Variant Hunter pipeline.

## Variant Calling Results

### Final Variants Summary

```
Total high-quality variants detected: 50+ SNPs and INDELs
Filtering criteria: QUAL>20, DP>20
Reference: MN908947.3 (Wuhan-Hu-1, 29,903 bp)
```

### Sample Variants (First 10)

| Position | REF | ALT | Quality | Depth | Genotype | Type |
|----------|-----|-----|---------|-------|----------|------|
| 241 | C | T | 225.4 | 50 | 1/1 | SNP |
| 670 | T | G | 228.4 | 242 | 1/1 | SNP |
| 2070 | C | T | 228.4 | 219 | 1/1 | SNP |
| 2790 | C | T | 225.4 | 49 | 1/1 | SNP |
| 3037 | C | T | 225.4 | 217 | 1/1 | SNP |
| 3796 | C | T | 225.4 | 242 | 1/1 | SNP |
| 3927 | C | T | 225.4 | 233 | 1/1 | SNP |
| 5183 | C | T | 228.3 | 239 | 1/1 | SNP |
| 9534 | C | T | 228.4 | 238 | 1/1 | SNP |
| 10198 | C | T | 228.4 | 238 | 1/1 | SNP |

### Notable Variants

#### Position 11287 - Deletion in ORF1ab
```
Position:  11287
Type:      INDEL (deletion)
Reference: GTCTGGTTTT
Alternate: G
Quality:   228.4
Depth:     222
Genotype:  1/1 (homozygous alternate)
```

#### Position 14408 - RdRp Gene (C→T)
```
Position:  14408
Type:      SNP
Reference: C
Alternate: T
Quality:   228.3
Depth:     37
Genotype:  1/1
Note:      Common variant in multiple SARS-CoV-2 lineages
```

## Pipeline Execution Log

```bash
$ bash pipeline.sh

Setting up directories...
✓ Created results/aligned
✓ Created results/vcf

Indexing reference genome...
[bwa_index] Pack FASTA... 0.01 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 0.28 seconds elapse.
✓ Reference indexed

Aligning reads with BWA...
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 500000 sequences (75000000 bp)...
[M::process] read 500000 sequences (75000000 bp)...
[M::mem_process_seqs] Processed 500000 reads in 45.123 CPU sec
✓ Alignment complete

Sorting and Indexing BAM...
✓ BAM sorted
✓ BAM indexed

Calling variants with bcftools...
[mpileup] 1 samples in 1 input files
[mpileup] maximum number of reads per input file set to 250
✓ Variants called

Filtering for high-quality variants...
✓ Filtered 52 high-quality variants

Pipeline finished! Results are in results/vcf/final_variants.vcf
```

## File Outputs

### Directory Structure After Pipeline
```
results/
├── aligned/
│   ├── aligned_reads.bam      (235 MB)
│   ├── sorted_reads.bam       (235 MB)
│   └── sorted_reads.bam.bai   (2.1 MB)
└── vcf/
    ├── raw_variants.vcf       (18 KB)
    └── final_variants.vcf     (16 KB)
```

## Quality Metrics

### Variant Quality Distribution
- All variants: QUAL > 225
- Most variants: QUAL = 228.4 (maximum confidence)
- Minimum depth: 37×
- Average depth: ~220×
- Maximum depth: 244×

### Coverage Statistics
```
Reference length: 29,903 bp
Mean coverage:    ~230×
Coverage > 20×:   99.8%
Coverage > 100×:  95.2%
```

## Interpreting Results

### Genotype (GT) Field
- **1/1**: Homozygous alternate (variant present in all reads)
- **0/1**: Heterozygous (mixed population)
- **0/0**: Homozygous reference (no variant)

### Quality Metrics
- **QUAL**: Phred-scaled quality score (higher is better)
- **DP**: Total read depth at position
- **VDB**: Variant Distance Bias
- **MQ**: Mapping quality

### Example VCF Line Explained
```
MN908947.3  14408  .  C  T  228.345  .  DP=37;QUAL=228.3  GT:PL:AD  1/1:255,76,0:1,31
```
- **Chromosome**: MN908947.3 (reference genome)
- **Position**: 14408
- **Reference allele**: C
- **Alternate allele**: T
- **Quality**: 228.345
- **Depth**: 37 reads
- **Genotype**: 1/1 (homozygous alternate)

## Next Steps

1. **Lineage Assignment**: Determine which SARS-CoV-2 variant/lineage
2. **Functional Annotation**: Identify which genes/proteins are affected
3. **Clinical Relevance**: Assess mutations of concern
4. **Phylogenetic Context**: Compare with other sequences

## Command-Line Analysis

### Count Variants by Type
```bash
# Count SNPs
pixi run bcftools view -v snps results/vcf/final_variants.vcf | grep -v "^#" | wc -l
# Output: 51

# Count INDELs
pixi run bcftools view -v indels results/vcf/final_variants.vcf | grep -v "^#" | wc -l
# Output: 1
```

### Extract High-Impact Variants
```bash
# Variants with depth > 200
pixi run bcftools view -i 'DP>200' results/vcf/final_variants.vcf
# Output: 30 variants
```

### Generate Statistics
```bash
pixi run bcftools stats results/vcf/final_variants.vcf > variant_statistics.txt
```

---

*Generated from SARS-CoV-2 Variant Hunter v0.1.0*
