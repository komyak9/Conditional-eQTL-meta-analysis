import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans
import pandas as pd
import os

def load_data():
    file_path = '/gpfs/space/projects/genomic_references/GEUVADIS/metadata/GEUVADIS.tsv'
    df_names = pd.read_csv(file_path, sep='\t')

    pr = genotypeio.PlinkReader("/gpfs/space/projects/genomic_references/GEUVADIS/genotypes/30x-grch38/plink/GEUVADIS")
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed("phenotype_preprocessed.bed.gz")
    for index, row in df_names.iterrows():
        phenotype_df = phenotype_df.rename(columns={row['sample_id']: row['genotype_id']})

    print('Genotypes and phenotypes are loaded.')
    return genotype_df, variant_df, phenotype_df, phenotype_pos_df


def main():
    # Set display options to show all columns without truncation
    pd.set_option('display.max_columns', None)
    pd.set_option('display.expand_frame_repr', False)

    # Load files and covariants
    genotype_df, variant_df, phenotype_df, phenotype_pos_df = load_data()
    covariants_df = pd.read_csv("covariants.tsv", sep='\t')
    covariants_df.set_index("iid", inplace=True)

    ##############################
    ### Signal 1 (covariant 1) ###
    ##############################

    # Extract 2 and 3 variants
    signals_2_3_df = covariants_df[['chr19_51612691_C_G', 'chr19_50966950_G_GT']]
    
    # Run trans eQTL mapping for the SIGLEC14 gene only
    # No p_val threshold to obtain results for all variants
    eQTL_result_df = trans.map_trans(genotype_df, phenotype_df.loc[phenotype_pos_df.index=='ENSG00000254415'],
                                     covariates_df=signals_2_3_df, pval_threshold=1, maf_threshold=0.01)
    
    # Extrant all non-cis regions for the SIGLEC14 gene
    trans_eQTL_result_df = trans.filter_cis(eQTL_result_df,
                                            phenotype_pos_df.loc[phenotype_pos_df.index=='ENSG00000254415'].T.to_dict(), 
                                            variant_df, window=1000000)
    
    # Remove all non-cis regions from the mapping result
    # to get only cis mappings for the SIGLEC14 gene
    merged_df = eQTL_result_df.merge(trans_eQTL_result_df, indicator=True, how='left')
    cis_result_df = merged_df[merged_df['_merge'] == 'left_only'].drop(columns='_merge')

    # Save to a file
    cis_result_df.to_csv("independent_eQTL_signals_2_3_p_val_1.tsv", sep='\t', index=True)

    ##############################
    ### Signal 2 (covariant 2) ###
    ##############################

    signals_1_3_df = covariants_df[['chr19_51627384_T_A', 'chr19_50966950_G_GT']]
    
    eQTL_result_df = trans.map_trans(genotype_df, phenotype_df.loc[phenotype_pos_df.index=='ENSG00000254415'],
                                     covariates_df=signals_1_3_df, pval_threshold=1, maf_threshold=0.01)
    trans_eQTL_result_df = trans.filter_cis(eQTL_result_df,
                                            phenotype_pos_df.loc[phenotype_pos_df.index=='ENSG00000254415'].T.to_dict(), 
                                            variant_df, window=1000000)
    merged_df = eQTL_result_df.merge(trans_eQTL_result_df, indicator=True, how='left')
    cis_result_df = merged_df[merged_df['_merge'] == 'left_only'].drop(columns='_merge')

    cis_result_df.to_csv("independent_eQTL_signals_1_3_p_val_1.tsv", sep='\t', index=True)

    ##############################
    ### Signal 3 (covariant 3) ###
    ##############################

    signals_1_2_df = covariants_df[['chr19_51627384_T_A', 'chr19_51612691_C_G']]
    
    eQTL_result_df = trans.map_trans(genotype_df, phenotype_df.loc[phenotype_pos_df.index=='ENSG00000254415'],
                                     covariates_df=signals_1_2_df, pval_threshold=1, maf_threshold=0.01)
    trans_eQTL_result_df = trans.filter_cis(eQTL_result_df,
                                            phenotype_pos_df.loc[phenotype_pos_df.index=='ENSG00000254415'].T.to_dict(), 
                                            variant_df, window=1000000)
    merged_df = eQTL_result_df.merge(trans_eQTL_result_df, indicator=True, how='left')
    cis_result_df = merged_df[merged_df['_merge'] == 'left_only'].drop(columns='_merge')
    
    cis_result_df.to_csv("independent_eQTL_signals_1_2_p_val_1.tsv", sep='\t', index=True)


if __name__ == "__main__":
    print(tensorqtl.__version__)
    main()