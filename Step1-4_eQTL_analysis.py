import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis
import pandas as pd
import os

def load_data():
    # Reading the file with the sample names
    # that is needed to map sample_ids (ERR) to genotype_ids.
    file_path = '/gpfs/space/projects/genomic_references/GEUVADIS/metadata/GEUVADIS.tsv'
    df_names = pd.read_csv(file_path, sep='\t')

    print('Loading genotypes and variants...')
    pr = genotypeio.PlinkReader("/gpfs/space/projects/genomic_references/GEUVADIS/genotypes/30x-grch38/plink/GEUVADIS")
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    print('Genotypes and variants are loaded.')

    print('Loading phenotypes, mapping ids...')
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed("phenotype_preprocessed.bed.gz")
    for index, row in df_names.iterrows():
        phenotype_df = phenotype_df.rename(columns={row['sample_id']: row['genotype_id']})
    print('Phenotypes are loaded.')

    return genotype_df, variant_df, phenotype_df, phenotype_pos_df


def main():
    # Set display options to show all columns without truncation
    pd.set_option('display.max_columns', None)
    pd.set_option('display.expand_frame_repr', False)

    # Reading data from files
    genotype_df, variant_df, phenotype_df, phenotype_pos_df = load_data()
    

    ##############
    ### Step 1 ###
    ##############

    path = "cis_eQTL_analysis.tsv"
    # Performing normal eQTL analysis
    # If it's beed performed, then read the results from a file
    if os.path.exists(path):
        eQTL_result_df = pd.read_csv(path, sep='\t')
        eQTL_result_df = eQTL_result_df.set_index('phenotype_id')
    else:
        eQTL_result_df = cis.map_cis(genotype_df, variant_df, phenotype_df.loc[phenotype_pos_df.index=='ENSG00000254415'],
                                    phenotype_pos_df, maf_threshold=0.01, window=1000000)
        eQTL_result_df.to_csv(path, sep='\t', index=True)   # Save results to a file


    ##############
    ### Step 2 ###
    ##############

    # Extracting results for the SIGLEC14 gene
    eQTL_SIGLEC14_df = eQTL_result_df
    #eQTL_SIGLEC14_df = eQTL_result_df[eQTL_result_df.index=='ENSG00000254415']
    print(f"SIGLEC QTL analysis results 1:\n{eQTL_SIGLEC14_df.head()}\n\n")


    ##################
    ### Steps 3, 4 ###
    ##################

    # Extracting covariates
    cov_df = genotype_df.loc[eQTL_SIGLEC14_df['variant_id']].T

    counter = 2
    while eQTL_SIGLEC14_df['pval_nominal'].iloc[0] < 1e-5:
        # Performing conditional eQTL analysis
        eQTL_result_df = cis.map_cis(genotype_df, variant_df, phenotype_df.loc[phenotype_pos_df.index=='ENSG00000254415'],
                                     phenotype_pos_df, covariates_df=cov_df, maf_threshold=0.01, window=1000000)
        
        # Extracting results for the SIGLEC14 gene
        eQTL_SIGLEC14_df = eQTL_result_df
        #eQTL_SIGLEC14_df = eQTL_result_df[eQTL_result_df.index=='ENSG00000254415']
        print(f"SIGLEC QTL analysis {counter}:\n{eQTL_SIGLEC14_df.head()}\n\n")
        counter += 1

        # Extracting covariates and updating the covariates dataframe
        cov_df_temp = genotype_df.loc[eQTL_SIGLEC14_df['variant_id']].T
        cov_df = pd.merge(cov_df, cov_df_temp, left_index=True, right_index=True)
    
    cov_df.to_csv("covariants.tsv", sep='\t', index=True)


if __name__ == "__main__":
    print(tensorqtl.__version__)
    main()