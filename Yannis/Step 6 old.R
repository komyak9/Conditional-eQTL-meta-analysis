library("dplyr")

library("coloc")

###########################

#import step5 signals ME

signal1 = readr::read_tsv("../Desktop/Yannis_Syncthing/Uni/Spring_23/Bioinformatics/Conditional-eQTL-meta-analysis/Yannis/independent_eQTL_signals_2_3_p_val_1.tsv")
  signal1 = signal1[, -1] ##remove first column
  
signal2 = readr::read_tsv("../Desktop/Yannis_Syncthing/Uni/Spring_23/Bioinformatics/Conditional-eQTL-meta-analysis/Yannis/independent_eQTL_signals_1_3_p_val_1.tsv")
  signal2 = signal2[, -1]
  
signal3 = readr::read_tsv("../Desktop/Yannis_Syncthing/Uni/Spring_23/Bioinformatics/Conditional-eQTL-meta-analysis/Yannis/independent_eQTL_signals_1_2_p_val_1.tsv")
  signal3 = signal3[, -1]
  

#Import marginal sumstats KAUR

marginal_eQTL_sumstats = readr::read_tsv("Bioinformatics_project_data/sumstats/QTD000544.cc.tsv.gz") ###

hal_marginal = dplyr::filter(marginal_eQTL_sumstats, molecular_trait_id == "ENSG00000084110") %>% ##select one trait
  
  dplyr::select(-rsid) %>% ##remove rsid column
  
  dplyr::distinct() ##remove duplicate rows

amdhd1_marginal = dplyr::filter(marginal_eQTL_sumstats, molecular_trait_id == "ENSG00000139344") %>%
  
  dplyr::select(-rsid) %>%
  
  dplyr::distinct()

###########################
###########################

#import QTL LBF ME
lbf_eQTL = readr::read_tsv("../Desktop/Yannis_Syncthing/Uni/Spring_23/Bioinformatics/Conditional-eQTL-meta-analysis/Yannis/QTD000584.lbf_variable.txt.gz")

siglec_lbf = dplyr::filter(lbf_eQTL, molecular_trait_id == "SIGLEC14.8248.222.3..1")

siglec_lbf2 = dplyr::filter(lbf_eQTL, molecular_trait_id == "SIGLEC14.5125.6.3..1")



#Import LBF sumstats KAUR

lbf_eQTL_sumstats = readr::read_tsv("Bioinformatics_project_data/sumstats/QTD000544.lbf_variable.txt.gz") ###

hal_lbf = dplyr::filter(lbf_eQTL_sumstats, molecular_trait_id == "ENSG00000084110")

amdhd1_lbf = dplyr::filter(lbf_eQTL_sumstats, molecular_trait_id == "ENSG00000139344")

###################################
###################################

#Import VitD KAUR

VitD_marginal = readr::read_tsv("Bioinformatics_project_data/sumstats/VitD.coloc3_combined.tsv.gz") %>% ###???
  
  dplyr::filter(region == "chr12_94482127-97481554") %>%
  
  dplyr::filter(maf > 0)

VitD_lbf = readr::read_tsv("Bioinformatics_project_data/sumstats/VitD.coloc5_combined.tsv.gz") %>%  ###
  
  dplyr::filter(region == "chr12_94482127-97481554")


###################################
###################################

#Make LBF matrices for coloc.bf_bf ME

siglec_lbf_mat = as.matrix(dplyr::select(siglec_lbf, lbf_variable1:lbf_variable10))
row.names(siglec_lbf_mat) = siglec_lbf$variant
siglec_lbf_mat = t(siglec_lbf_mat)

siglec_lbf_mat2 = as.matrix(dplyr::select(siglec_lbf2, lbf_variable1:lbf_variable10))
row.names(siglec_lbf_mat2) = siglec_lbf2$variant
siglec_lbf_mat2 = t(siglec_lbf_mat2)


#Make LBF matrices for coloc.bf_bf KAUR

vitd_lbf_mat = as.matrix(dplyr::select(VitD_lbf, lbf_variable1:lbf_variable10))
row.names(vitd_lbf_mat) = VitD_lbf$variant
vitd_lbf_mat = t(vitd_lbf_mat)

hal_lbf_mat = as.matrix(dplyr::select(hal_lbf, lbf_variable1:lbf_variable10))
row.names(hal_lbf_mat) = hal_lbf$variant
hal_lbf_mat = t(hal_lbf_mat)

amdhd1_lbf_mat = as.matrix(dplyr::select(amdhd1_lbf, lbf_variable1:lbf_variable10))
row.names(amdhd1_lbf_mat) = amdhd1_lbf$variant
amdhd1_lbf_mat = t(amdhd1_lbf_mat)

#####################################
#####################################

#Perform coloc ME

siglec_coloc = coloc.bf_bf(siglec_lbf_mat, siglec_lbf_mat) ###????
dplyr::filter(siglec_coloc$summary, PP.H4.abf > 0.9)

siglec_coloc2 = coloc.bf_bf(siglec_lbf_mat2, siglec_lbf_mat2)
dplyr::filter(siglec_coloc2$summary, PP.H4.abf > 0.9)


#Perfrom coloc KAUR

hal_coloc = coloc.bf_bf(vitd_lbf_mat, hal_lbf_mat)
dplyr::filter(hal_coloc$summary, PP.H4.abf > 0.9)


amdhd1_coloc = coloc.bf_bf(vitd_lbf_mat, amdhd1_lbf_mat)
dplyr::filter(amdhd1_coloc$summary, PP.H4.abf > 0.9)

#####################################
#####################################

#Run coloc.abf ME

signal1_list = list(beta = signal1$b, 
                varbeta = signal1$b_se^2, 
                N = length(signal1), 
                MAF = signal1$af, 
                snp = signal1$variant_id, 
                type = "quant")

coloc::check_dataset(signal1_list)


signal2_list = list(beta = signal2$b, 
                    varbeta = signal2$b_se^2, 
                    N = length(signal2), 
                    MAF = signal2$af, 
                    snp = signal2$variant_id, 
                    type = "quant")

coloc::check_dataset(signal2_list)


signal3_list = list(beta = signal3$b, 
                    varbeta = signal3$b_se^2, 
                    N = length(signal3), 
                    MAF = signal3$af, 
                    snp = signal3$variant_id, 
                    type = "quant")

coloc::check_dataset(signal3_list)


signal1_coloc_abf = coloc.abf(signal1_list, signal2_list) ##???



#Run coloc.abf KAUR

hal_list = list(beta = hal_marginal$beta, 
                varbeta = hal_marginal$se^2, 
                N = hal_marginal$an/2, 
                MAF = hal_marginal$maf, 
                snp = hal_marginal$variant, 
                type = "quant")

coloc::check_dataset(hal_list)


amdhd1_list = list(beta = amdhd1_marginal$beta, 
                   varbeta = amdhd1_marginal$se^2, 
                   N = amdhd1_marginal$an/2, 
                   MAF = amdhd1_marginal$maf, 
                   snp = amdhd1_marginal$variant, 
                   type = "quant")

coloc::check_dataset(amdhd1_list)


vitd_list = list(beta = VitD_marginal$beta, 
                 varbeta = VitD_marginal$se^2, 
                 N = VitD_marginal$an/2, 
                 MAF = VitD_marginal$maf, 
                 snp = VitD_marginal$variant, 
                 type = "quant")

coloc::check_dataset(vitd_list)


hal_coloc_abf = coloc.abf(hal_list, vitd_list)

amdhd1_coloc_abf = coloc.abf(amdhd1_list, vitd_list)


#####################################
#####################################

#EXPERIMENTAL: perform coloc between lbf and labf values ME

signal1_labf_mat = dplyr::as_tibble(signal1_coloc_abf$results) %>%
  dplyr::transmute(lbf_variable1 = lABF.df1, lbf_variable2 = 0) %>% as.matrix()
row.names(signal1_labf_mat) = signal1_coloc_abf$results$snp
signal1_labf_mat = t(signal1_labf_mat)

signal1_lbf_labf_coloc = coloc.bf_bf(siglec_lbf_mat, siglec_lbf_mat)
dplyr::filter(signal1_lbf_labf_coloc$summary, PP.H4.abf > 0.9)


#EXPERIMENTAL: perform coloc between lbf and labf values KAUR

HAL_labf_mat = dplyr::as_tibble(hal_coloc_abf$results) %>% 
  dplyr::transmute(lbf_variable1 = lABF.df1, lbf_variable2 = 0) %>% as.matrix()
row.names(HAL_labf_mat) = hal_coloc_abf$results$snp
HAL_labf_mat = t(HAL_labf_mat)


lbf_labf_coloc = coloc.bf_bf(vitd_lbf_mat, HAL_labf_mat)
dplyr::filter(lbf_labf_coloc$summary, PP.H4.abf > 0.9)