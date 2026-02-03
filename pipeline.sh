#!/bin/bash

# ==========================================
# SARS-CoV-2 Variant Calling Pipeline
# Author: Nasir Nesirli
# Description: Automates QC, Alignment, and Variant Calling
# ==========================================

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Configuration
REFERENCE="data/reference/MN908947.3.fasta"
R1_FASTQ="data/raw_fastq/sample_R1.fastq.gz"
R2_FASTQ="data/raw_fastq/sample_R2.fastq.gz"
ALIGNED_DIR="results/aligned"
VCF_DIR="results/vcf"

# Color codes for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Functions
log_info() {
    echo -e "${GREEN}✓${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}⚠${NC} $1"
}

log_error() {
    echo -e "${RED}✗${NC} $1"
    exit 1
}

# Check prerequisites
echo "==========================================
SARS-CoV-2 Variant Calling Pipeline
==========================================="

# Check if reference genome exists
if [ ! -f "$REFERENCE" ]; then
    log_error "Reference genome not found: $REFERENCE"
fi

# Check if FASTQ files exist
if [ ! -f "$R1_FASTQ" ]; then
    log_error "Forward reads not found: $R1_FASTQ"
fi

if [ ! -f "$R2_FASTQ" ]; then
    log_error "Reverse reads not found: $R2_FASTQ"
fi

log_info "All input files found"

# 1. Setup Directories
echo ""
echo "Step 1: Setting up directories..."
mkdir -p "$ALIGNED_DIR" "$VCF_DIR"
log_info "Directories created"

# 2. Index the Reference
echo ""
echo "Step 2: Indexing reference genome..."
if [ -f "${REFERENCE}.bwt" ]; then
    log_warn "Reference already indexed, skipping..."
else
    bwa index "$REFERENCE"
    log_info "Reference indexed successfully"
fi

# 3. Align Reads to Reference
echo ""
echo "Step 3: Aligning reads with BWA-MEM..."
bwa mem "$REFERENCE" "$R1_FASTQ" "$R2_FASTQ" | \
    samtools view -b -o "${ALIGNED_DIR}/aligned_reads.bam"
log_info "Alignment complete: ${ALIGNED_DIR}/aligned_reads.bam"

# 4. Sort and Index BAM
echo ""
echo "Step 4: Sorting and indexing BAM file..."
samtools sort "${ALIGNED_DIR}/aligned_reads.bam" -o "${ALIGNED_DIR}/sorted_reads.bam"
log_info "BAM sorted: ${ALIGNED_DIR}/sorted_reads.bam"

samtools index "${ALIGNED_DIR}/sorted_reads.bam"
log_info "BAM indexed: ${ALIGNED_DIR}/sorted_reads.bam.bai"

# 5. Call Variants
echo ""
echo "Step 5: Calling variants with BCFtools..."
bcftools mpileup -f "$REFERENCE" "${ALIGNED_DIR}/sorted_reads.bam" | \
    bcftools call -mv -o "${VCF_DIR}/raw_variants.vcf"
log_info "Variants called: ${VCF_DIR}/raw_variants.vcf"

# 6. Filter Variants
echo ""
echo "Step 6: Filtering for high-quality variants..."
bcftools view -i 'QUAL>20 && DP>20' "${VCF_DIR}/raw_variants.vcf" > "${VCF_DIR}/final_variants.vcf"

# Count variants
VARIANT_COUNT=$(grep -v "^#" "${VCF_DIR}/final_variants.vcf" | wc -l | tr -d ' ')
log_info "Filtering complete: ${VARIANT_COUNT} high-quality variants detected"

# 7. Summary
echo ""
echo "=========================================="
echo "Pipeline completed successfully!"
echo "=========================================="
echo "Results:"
echo "  - Sorted BAM: ${ALIGNED_DIR}/sorted_reads.bam"
echo "  - Raw variants: ${VCF_DIR}/raw_variants.vcf"
echo "  - Final variants: ${VCF_DIR}/final_variants.vcf"
echo "  - Total variants: ${VARIANT_COUNT}"
echo ""
echo "Next steps:"
echo "  1. View variants: cat ${VCF_DIR}/final_variants.vcf"
echo "  2. Get statistics: bcftools stats ${VCF_DIR}/final_variants.vcf"
echo "  3. Check alignment: samtools flagstat ${ALIGNED_DIR}/sorted_reads.bam"
echo "=========================================="