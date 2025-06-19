# post_gwas_SNP_annotation
# Post-GWAS Analysis Pipeline (V.2)

This pipeline processes GWAS results to identify significant SNPs, annotate them with gene information, and find overlaps with fatty acid (FA)-related genes. It uses command-line tools (`awk`, `bedtools`) and Python (`pandas`) for processing.

## Directory Structure

- **Input Files**:
  - GWAS Results: `Y:\Bigdata\computing\Shima\b_carinata\post_gwas_26_05_2025\10_miss\GAPIT.Association.GWAS_Results.SUPER.X22_1.csv`
  - GFF3 File: `/home/mahmoudi/Bigdata/computing/Shima/b_carinata/gff_v2/genomic.gff`
  - Annotated FA Genes: `annotated_genes_plot.txt`
- **Output Files**:
  - `significant_snps.csv`: Filtered significant SNPs (P < 5e-5).
  - `significant_snps.bed`: Significant SNPs in BED format.
  - `genes.bed`: Gene features extracted from GFF3 in BED format.
  - `snp_gene_overlap.txt`: Direct overlaps between SNPs and genes.
  - `snp_gene_nearby.txt`: SNPs with nearby genes (within 1Mb window).
  - `FA_genes_near_SNP_from_blast.csv`: FA-related genes near significant SNPs.

## Prerequisites

- **Tools**:
  - `awk`
  - `bedtools` (v2.30.0 or later)
  - Python 3.x with `pandas`
- **System**: Linux/Unix environment (commands tested on Linux).
- **Input Format**:
  - GWAS CSV: Columns include SNP ID, chromosome, position, and p-value (column 4).
  - GFF3: Standard GFF3 format with gene features.
  - FA Genes: Tab-separated file with columns for chromosome, start, end, and gene ID.

## Pipeline Steps

### 1. Filter Significant SNPs

Filter SNPs with p-value < 5e-5 from GWAS results.

```bash
awk -F',' 'NR==1 || $4 < 5e-5' GAPIT.Association.GWAS_Results.SUPER.X22_1.csv > significant_snps.csv
```

### 2. Convert Significant SNPs to BED Format
Convert filtered SNPs to BED format (0-based, half-open).
```bash
awk -F',' 'NR>1 {OFS="\t"; print $2, $3-1, $3, $1}' significant_snps.csv > significant_snps.bed
```

### 3. Extract Gene Features from GFF3
Extract only gene entries from the GFF3 file and convert to BED format.
```bash

awk '$3 == "gene"' /home/mahmoudi/Bigdata/computing/Shima/b_carinata/gff_v2/genomic.gff | \
awk -F'\t' 'BEGIN{OFS="\t"} { 
  match($9, /ID=([^;]+)/, a); 
  print $1, $4-1, $5, a[1] 
}' > genes.bed
```
### 4. Annotate SNPs with Gene Information
## 4.1 Direct Overlaps
Find SNPs overlapping gene regions.
```bash

bedtools intersect -a significant_snps.bed -b genes.bed -wa -wb > snp_gene_overlap.txt
```
## 4.2 Nearby Genes
Find genes within 1Mb of SNPs.
```bash

bedtools window -a significant_snps.bed -b genes.bed -w 1000000 > snp_gene_nearby.txt
```
### 5. Summarize Nearby Genes
Load snp_gene_nearby.txt and print the number of nearby genes.
python

import pandas as pd

columns = [
    "SNP_chr", "SNP_start", "SNP_end", "SNP_ID",
    "Gene_chr", "Gene_start", "Gene_end", "Gene_ID",
    "Strand", "Annotation"
]

df = pd.read_csv("snp_gene_nearby.txt", sep="\t", header=None, names=columns)
print(f"Nearby genes found: {df.shape[0]}")

### 6. Identify FA-Related Genes Near SNPs
Match nearby genes with annotated FA genes from BLAST output.
python

import pandas as pd

# Load SNP-gene nearby file
nearby_df = pd.read_csv("snp_gene_nearby.txt", sep="\t", header=None,
                        names=["SNP_chr", "SNP_start", "SNP_end", "SNP_ID",
                               "Gene_chr", "Gene_start", "Gene_end", "Gene_ID",
                               "Strand", "Annotation"])

# Load annotated FA genes
fa_genes_df = pd.read_csv("annotated_genes_plot.txt", sep="\t", header=None,
                          names=["FA_chr", "FA_start", "FA_end", "FA_Gene_ID"])

# Standardize gene IDs for matching
nearby_df["Gene_ID_upper"] = nearby_df["Gene_ID"].str.upper()
fa_genes_df["FA_Gene_ID_upper"] = fa_genes_df["FA_Gene_ID"].str.upper()

# Merge to find matches
matched_fa_genes = nearby_df.merge(fa_genes_df, left_on="Gene_ID_upper", right_on="FA_Gene_ID_upper")

# Keep relevant columns
matched_fa_genes = matched_fa_genes[[
    "SNP_ID", "SNP_chr", "SNP_start", "Gene_ID", "Gene_chr", "Gene_start", "Gene_end", "FA_Gene_ID"
]]

# Save results
matched_fa_genes.to_csv("FA_genes_near_SNP_from_blast.csv", index=False)
print(f"{matched_fa_genes.shape[0]} fatty acid-related genes found near SNPs.")

Notes
Adjust file paths as needed for your environment.

The pipeline assumes the GWAS CSV has a specific structure (SNP ID in column 1, chromosome in column 2, position in column 3, p-value in column 4).

The GFF3 file must contain gene features with ID in the attributes field.

The FA genes file (annotated_genes_plot.txt) should be tab-separated with chromosome, start, end, and gene ID columns.

The 1Mb window for nearby genes can be adjusted by changing the -w parameter in bedtools window.

