#!/bin/bash

# ==========================================
# SARS-CoV-2 Variant Analysis Helper
# Quick commands to analyze pipeline results
# ==========================================

set -e

# Color codes
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo "=========================================="
echo "SARS-CoV-2 Variant Analysis Helper"
echo "=========================================="
echo ""

# Check if results exist
if [ ! -f "results/vcf/final_variants.vcf" ]; then
    echo "Error: results/vcf/final_variants.vcf not found"
    echo "Please run pipeline.sh first"
    exit 1
fi

# 1. Variant Summary
echo -e "${BLUE}1. Variant Summary${NC}"
echo "---"
TOTAL_VARIANTS=$(grep -v "^#" results/vcf/final_variants.vcf | wc -l | tr -d ' ')
echo "Total variants: $TOTAL_VARIANTS"

# Count SNPs (without INDEL flag)
SNPS=$(grep -v "^#" results/vcf/final_variants.vcf | grep -v "INDEL" | wc -l | tr -d ' ')
echo "SNPs: $SNPS"

# Count INDELs
INDELS=$(grep -v "^#" results/vcf/final_variants.vcf | grep "INDEL" | wc -l | tr -d ' ')
echo "INDELs: $INDELS"
echo ""

# 2. Quality Statistics
echo -e "${BLUE}2. Quality Metrics${NC}"
echo "---"
if command -v bcftools &> /dev/null; then
    echo "Average Quality:"
    grep -v "^#" results/vcf/final_variants.vcf | awk '{sum+=$6; count++} END {print "  Mean QUAL: " sum/count}'
    echo "Average Depth:"
    grep -v "^#" results/vcf/final_variants.vcf | grep -o "DP=[0-9]*" | cut -d'=' -f2 | awk '{sum+=$1; count++} END {print "  Mean DP: " sum/count}'
else
    echo "bcftools not found, skipping detailed stats"
fi
echo ""

# 3. Top 10 Variants by Position
echo -e "${BLUE}3. Top 10 Variants (by position)${NC}"
echo "---"
echo "Position | Ref | Alt | Quality | Depth | Type"
echo "---------|-----|-----|---------|-------|-----"
grep -v "^#" results/vcf/final_variants.vcf | head -10 | \
    awk '{
        indel = ($8 ~ /INDEL/) ? "INDEL" : "SNP";
        if (match($8, /DP=([0-9]+)/)) {
            dp = substr($8, RSTART+3, RLENGTH-3);
        } else {
            dp = "N/A";
        }
        printf "%-8s | %-3s | %-3s | %-7.1f | %-5s | %s\n", $2, $4, $5, $6, dp, indel
    }'
echo ""

# 4. Notable Positions
echo -e "${BLUE}4. Notable SARS-CoV-2 Mutations${NC}"
echo "---"
# Check for common variant positions
KNOWN_POSITIONS=(241 3037 14408 23403)
for POS in "${KNOWN_POSITIONS[@]}"; do
    VARIANT=$(grep -v "^#" results/vcf/final_variants.vcf | awk -v pos="$POS" '$2 == pos {print $4"->"$5}' || echo "Not found")
    if [ "$VARIANT" != "Not found" ]; then
        echo -e "${GREEN}✓${NC} Position $POS: $VARIANT (detected)"
    else
        echo "  Position $POS: No variant"
    fi
done
echo ""

# 5. Alignment Statistics
echo -e "${BLUE}5. Alignment Statistics${NC}"
echo "---"
if [ -f "results/aligned/sorted_reads.bam" ]; then
    if command -v samtools &> /dev/null; then
        echo "Total reads:"
        samtools view -c results/aligned/sorted_reads.bam || echo "  Error reading BAM file"
        echo ""
        echo "Coverage summary:"
        samtools coverage results/aligned/sorted_reads.bam 2>/dev/null | tail -1 | \
            awk '{printf "  Mean coverage: %.1fx\n  Coverage: %.1f%%\n", $7, $6}' || \
            echo "  Unable to calculate coverage"
    else
        echo "samtools not found"
    fi
else
    echo "BAM file not found"
fi
echo ""

# 6. File Sizes
echo -e "${BLUE}6. Output File Sizes${NC}"
echo "---"
if [ -d "results" ]; then
    du -h results/aligned/*.bam 2>/dev/null | awk '{print "  " $2 ": " $1}' || echo "  No BAM files"
    du -h results/vcf/*.vcf 2>/dev/null | awk '{print "  " $2 ": " $1}' || echo "  No VCF files"
fi
echo ""

echo "=========================================="
echo -e "${GREEN}Analysis Complete!${NC}"
echo "=========================================="
echo ""
echo "For more detailed analysis:"
echo "  - Variant stats: bcftools stats results/vcf/final_variants.vcf"
echo "  - Alignment stats: samtools flagstat results/aligned/sorted_reads.bam"
echo "  - View specific variant: bcftools view -r MN908947.3:14408 results/vcf/final_variants.vcf"
echo ""
