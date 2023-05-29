# Conditional-eQTL-meta-analysis
Bioinformatics course project

For steps 1-4 we use `tensorqtl` tool.

Required data types:

1. Phenotypes: BED format
2. Genotypes: PLINK format
3. Covariates: a tab-delimited text file (covariates x samples) or dataframe (samples x covariates), with row and column headers

For step 5 we use R tools.

##
## Step 1: Perform normal eQTL analysis for SIGLEC14 gene, testing all common (MAF > 1%) variants in +/- 1Mb window around the promoter of the gene.

#### Normal eQTL analysis

A normal eQTL analysis typically involves testing for genetic variants that are associated with changes in gene expression levels (expression quantitative trait loci or eQTLs) across the genome, without any assumptions about the relationship between the genetic variant and the gene's location. So, no other variables are used.

The term "eQTL" stands for expression quantitative trait loci, and it refers to genetic variants (usually single nucleotide polymorphisms, or SNPs) that are associated with differences in gene expression levels.

### 1.1. Create phenotypes file.

tensorqtl requires phenotype data be in .bed or .bed.gz format. Phenotype file must contain #Chr, start, end, TargetID columns in the mentioned order.

To create such a file, we processed and merged normalised gene expression data and covariates data from different files: chromosome number, phenotype position (start, and end is start+1) were taken from one file, and TargetIds from another file.

### 1.2. TensorQTL analysis.

Tensorqtl provides 2 analyses: cis-eQTL analysis and trans-eQTL analysis.

cis-eQTL analysis is used to identify genetic variants that affect the expression level of genes located nearby (usually within 1 Mb distance from the transcription start site(TSS)).

So, __cis-eQTL analysis__ is the type of analysis we use in our project. Parameters used: maf_threshold=0.01, window=1000000. Besides, we run the analysis only for 19th chromosome as the gene is located there.

```python
eQTL_result_df = cis.map_cis(genotype_df, variant_df, phenotype_df.loc[phenotype_pos_df['chr']=='19'], phenotype_pos_df, maf_threshold=0.01, window=1000000)
```

##
## Step 2: Identify the most strongly associated variant (lead variant).

At first, we need to extract results for the SIGLEC14 gene:
```python
eQTL_SIGLEC14_df = eQTL_result_df[eQTL_result_df.index=='ENSG00000254415']
```

By default, only the lead variant is displayed after the tensorqtl cis mapping. This, after the first step we get the following dataframe with the lead variant (chr19_51627384_T_A).

phenotype_id | num_var | beta_shape1 | beta_shape2 | true_df | pval_true_df | variant_id | tss_distance | ma_samples | ma_count | af | pval_nominal | slope | slope_se | pval_perm | pval_beta |
--------------- | ----- | ------- | ---------- | ---------- | -------------- | ------------------ | ------ | --- | --- | -------- | -------------- | --------- | -------- | ------ | --------
ENSG00000254415 | 14863 | 1.03784 | 2284.05127 | 397.486603 | 4.515546e-09 | chr19_51627384_T_A | -19442 | 105 | 111 | 0.124719 | 5.973180e-10 | -0.615058 | 0.097147 | 0.0001 | 0.000007

##
## Step 3. Perform conditional eQTL analysis for SIGLEC14 by adding the lead variant as a covariate into the model.
## Step 4: Repeat this until no significant (p < 1e-5) associations remain.

For these steps, we need to create a covariant dataframe. Dataframe row indexes must be sample ids, and each column - variant id. The values are the genotype dosages (0, 1, or 2) for the variant in each individual (sample id).

```python
cov_df = genotype_df.loc[eQTL_SIGLEC14_df['variant_id']].T
```

![image](https://github.com/komyak9/Conditional-eQTL-meta-analysis/assets/42679553/7b90b26c-6d20-45c7-8ebd-7d429f2471f7)

Having such a dataframe, we run the same cis-mapping tensorqtl method with maf_threshold=0.01, window=1000000 parameters and also add the covariant dataframe as a parameter: covariates_df=cov_df.

```python
eQTL_result_df = cis.map_cis(genotype_df, variant_df, phenotype_df.loc[phenotype_pos_df['chr']=='19'],
                                     phenotype_pos_df, covariates_df=cov_df, maf_threshold=0.01, window=1000000)
```

After getting results, we repeat the process: we take the most strongly associated variant, add it to the covariants dataframe, run the conditional eQTL analysis. We stop when the p-value becomes >= 1e5 (that means that no significant associations remain).

```python
    while eQTL_SIGLEC14_df['pval_nominal'].iloc[0] < 1e-5:
        # Performing conditional eQTL analysis
        eQTL_result_df = cis.map_cis(genotype_df, variant_df, phenotype_df.loc[phenotype_pos_df['chr']=='19'],
                                     phenotype_pos_df, covariates_df=cov_df, maf_threshold=0.01, window=1000000)
        
        # Extracting results for the SIGLEC14 gene
        eQTL_SIGLEC14_df = eQTL_result_df[eQTL_result_df.index=='ENSG00000254415']

        # Extracting covariates and updating the covariates dataframe
        cov_df_temp = genotype_df.loc[eQTL_SIGLEC14_df['variant_id']].T
        cov_df = pd.merge(cov_df, cov_df_temp, left_index=True, right_index=True)
```

##
## Step 5: Identify all-but-one conditionally indepedent summary statistics.

If in the step 3 conditionally independent signals are identified, then at this step conditioning has to be performed three times:
1. for Signal 1, condition on Signals 2 and 3 (add both lead SNPs as covariates).
2. for Signal 2, condition on Signals 1 and 3.
3. for Signal 3, condition on Signals 1 and 2.

##
