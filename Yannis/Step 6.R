library("dplyr")

library("coloc")

library("ggplot2")



#Import pqtl data

pqtl_data = readr::read_tsv("R/Bioinformatics_project_data/step6/QTD000584.lbf_variable.txt.gz")

siglec14_df = dplyr::filter(pqtl_data, molecular_trait_id == "SIGLEC14.8248.222.3..1")



#Import conditional sumstats

signal2 = readr::read_tsv("R/Bioinformatics_project_data/step6/independent_eQTL_signals_1_3_p_val_1.tsv") %>% dplyr::filter(variant_id %in% siglec14_df$variant)

signal1 = readr::read_tsv("R/Bioinformatics_project_data/step6/independent_eQTL_signals_2_3_p_val_1.tsv") %>% dplyr::filter(variant_id %in% siglec14_df$variant)

signal3 = readr::read_tsv("R/Bioinformatics_project_data/step6/independent_eQTL_signals_1_2_p_val_1.tsv")



#Make approximate bayes factors

#Run coloc.abf

signal1_list = list(beta = signal1$b, 
                    
                    varbeta = signal1$b_se^2, 
                    
                    N = rep(445, length(signal1$b)), 
                    
                    MAF = pmin(signal1$af, 1-signal1$af), 
                    
                    snp = signal1$variant_id, 
                    
                    type = "quant")

coloc::check_dataset(signal1_list)



signal2_list = list(beta = signal2$b, 
                    
                    varbeta = signal2$b_se^2, 
                    
                    N = rep(445, length(signal2$b)), 
                    
                    MAF = pmin(signal2$af, 1-signal2$af), 
                    
                    snp = signal2$variant_id, 
                    
                    type = "quant")

coloc::check_dataset(signal2_list)



coloc_df = coloc.abf(signal1_list, signal2_list)

labf_df = dplyr::transmute(coloc_df$results, variant = snp, labf_variable1 = lABF.df1, labf_variable2 = lABF.df2) %>% 
  
  dplyr::as_tibble()



#Make protein lbf matrix

protein_lbf_mat = as.matrix(dplyr::select(siglec14_df, lbf_variable1:lbf_variable10))

row.names(protein_lbf_mat) = siglec14_df$variant

protein_lbf_mat = t(protein_lbf_mat)



#Make gene labf variable matrix

gene_lbf_mat = as.matrix(dplyr::select(labf_df, labf_variable1:labf_variable2))

row.names(gene_lbf_mat) = labf_df$variant

gene_lbf_mat = t(gene_lbf_mat)



#Check matrices

gene_lbf_mat[1:2, 1:10]

protein_lbf_mat[1:10, 1:10]



#Perform colocalisation

lbf_labf_coloc = coloc.bf_bf(gene_lbf_mat, protein_lbf_mat)

dplyr::filter(lbf_labf_coloc$summary, PP.H4.abf > 0.9)



#Visualise LABF vs LBF

full_df = dplyr::left_join(siglec14_df, labf_df, by = "variant")

ggplot(full_df, aes(x = lbf_variable1, y = labf_variable1)) + geom_point()

ggplot(full_df, aes(x = lbf_variable2, y = labf_variable2)) + geom_point()



#Import gene expression lbfs

geuvadis_lbf_data = readr::read_tsv("R/Bioinformatics_project_data/step6/QTD000110.lbf_variable.txt.gz")

siglec14_eqtl_lbf = dplyr::filter(geuvadis_lbf_data, molecular_trait_id == "ENSG00000254415")



#Visualise gene expression LBF vs LABF

full_df2 = dplyr::left_join(siglec14_eqtl_lbf, labf_df, by = "variant")

ggplot(full_df2, aes(x = lbf_variable1, y = labf_variable1)) + geom_point()

ggplot(full_df2, aes(x = lbf_variable2, y = labf_variable2)) + geom_point()

ggplot(full_df2, aes(x = lbf_variable1, y = labf_variable2)) + geom_point()