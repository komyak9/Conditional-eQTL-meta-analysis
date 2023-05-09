# Conditional-eQTL-meta-analysis
Bioinformatics course project

For steps 1-4 we use `tensorqtl` tool.

Required data types:

1. Phenotypes: BED format
2. Genotypes: PLINK format
3. Covariates: a tab-delimited text file (covariates x samples) or dataframe (samples x covariates), with row and column headers

## Steps

### Step 1: Perform normal eQTL analysis for SIGLEC14 gene, testing all common (MAF > 1%) variants in +/- 1Mb window around the promoter of the gene.

Tensorqtl provides 2 analyses: cis-eQTL analysis and trans-eQTL analysis.

cis-eQTL analysis is used to identify genetic variants that affect the expression level of genes located nearby (usually within 1 Mb distance from the transcription start site(TSS)).

So, __cis-eQTL analysis__ is the analysis we use in our project.

#### Normal eQTL analysis

A normal eQTL analysis typically involves testing for genetic variants that are associated with changes in gene expression levels (expression quantitative trait loci or eQTLs) across the genome, without any assumptions about the relationship between the genetic variant and the gene's location. So, no other variables are used.

The term "eQTL" stands for expression quantitative trait loci, and it refers to genetic variants (usually single nucleotide polymorphisms, or SNPs) that are associated with differences in gene expression levels.



##
### Step 2: Identify the most strongly associated variant (lead variant).
##
### Step 3. Perform conditional eQTL analysis for SIGLEC14 by adding the lead variant as a covariate into the model.
##
### Step 4: Repeat this until no significant (p < 1e-5) associations remain.
##
### Step 5: Identify all-but-one conditionally indepedent summary statistics.

If in the step 3 conditionally independent signals are identified, then at this step conditioning has to be performed three times:
1. for Signal 1, condition on Signals 2 and 3 (add both lead SNPs as covariates).
2. for Signal 2, condition on Signals 1 and 3.
3. for Signal 3, condition on Signals 1 and 2.
