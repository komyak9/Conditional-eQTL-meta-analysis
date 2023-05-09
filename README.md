# Conditional-eQTL-meta-analysis
Bioinformatics course project

## Steps

### Step 1: Perform normal eQTL analysis for SIGLEC14 gene, testing all common (MAF > 1%) variants in +/- 1Mb window around the promoter of the gene.

### Step 2: Identify the most strongly associated variant (lead variant).

### Step 3. Perform conditional eQTL analysis for SIGLEC14 by adding the lead variant as a covariate into the model.

### Step 4: Repeat this until no significant (p < 1e-5) associations remain.

### Step 5: Identify all-but-one conditionally indepedent summary statistics.

If in the step 3 conditionally independent signals are identified, then at this step conditioning has to be performed three times:
1. for Signal 1, condition on Signals 2 and 3 (add both lead SNPs as covariates).
2. for Signal 2, condition on Signals 1 and 3.
3. for Signal 3, condition on Signals 1 and 2.
