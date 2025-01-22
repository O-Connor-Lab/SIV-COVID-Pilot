rm(list = ls())

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ANCOMBC")
#BiocManager::install("mia", force = TRUE)


library(mia)
library(ANCOMBC)


# Stool
tse_stool = importQIIME2(
  featureTableFile = "Data/ancombc2-data/Stool-FT.qza",
  taxonomyTableFile = "Data/taxonomy_2.qza",
  sampleMetaFile = "Data/ancombc2-data/Stool-MD.txt",
  refSeqFile = "Data/ancombc2-data/Stool-RS.qza",
  phyTreeFile = "Data/ancombc2-data/Stool-SEPPFT.qza"
)

tse_stool$time = factor(tse_stool$time, levels = c("D-7", "D3", "D5", "D7", "D10", "D14"))

output_genus_stool = ancombc2(data = tse_stool, assay_name = "counts", tax_level = "Genus",
                              fix_formula = "time",
                              rand_formula = "(1 | patient)",
                              p_adj_method = "BH", pseudo_sens = TRUE,
                              prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                              group = "time", struc_zero = TRUE, neg_lb = TRUE,
                              alpha = 0.05, n_cl = 2, verbose = TRUE,
                              global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                              iter_control = list(tol = 1e-2, max_iter = 20, 
                                                  verbose = TRUE),
                              em_control = list(tol = 1e-5, max_iter = 100),
                              lme_control = lme4::lmerControl(check.nobs.vs.nRE= "ignore"),
                              mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                              trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                          nrow = 2, 
                                                                          byrow = TRUE)),
                                                   node = list(2),
                                                   solver = "ECOS",
                                                   B = 100))

res_prim_stool <- output_genus_stool$res 

# Rectal Swab

tse_rectal = mia::loadFromQIIME2(
  featureTableFile = "Data/ancombc2-data/RectalSwab-FT.qza",
  taxonomyTableFile = "Data/taxonomy_2.qza",
  sampleMetaFile = "Data/ancombc2-data/RectalSwab-MD.txt",
  refSeqFile = "Data/ancombc2-data/RectalSwab-RS.qza",
  phyTreeFile = "Data/ancombc2-data/RectalSwab-SEPPFT.qza"
)

tse_rectal$time = factor(tse_rectal$time, levels = c("D-7", "D3", "D5", "D7", "D10", "D14"))

output_genus_rectal = ancombc2(data = tse_rectal, assay_name = "counts", tax_level = "Genus",
                               fix_formula = "time",
                               rand_formula = "(1 | patient)",
                               p_adj_method = "BH", pseudo_sens = TRUE,
                               prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                               group = "time", struc_zero = TRUE, neg_lb = TRUE,
                               alpha = 0.05, n_cl = 2, verbose = TRUE,
                               global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                               iter_control = list(tol = 1e-2, max_iter = 20, 
                                                   verbose = TRUE),
                               em_control = list(tol = 1e-5, max_iter = 100),
                               lme_control = lme4::lmerControl(check.nobs.vs.nRE= "ignore"),
                               mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                               trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                           nrow = 2, 
                                                                           byrow = TRUE)),
                                                    node = list(2),
                                                    solver = "ECOS",
                                                    B = 100))

res_prim_rectal <- output_genus_rectal$res 

# Throat Swab
tse_throat = mia::loadFromQIIME2(
  featureTableFile = "Data/ancombc2-data/ThroatSwab-FT.qza",
  taxonomyTableFile = "Data/taxonomy_2.qza",
  sampleMetaFile = "Data/ancombc2-data/ThroatSwab-MD.txt",
  refSeqFile = "Data/ancombc2-data/ThroatSwab-RS.qza",
  phyTreeFile = "Data/ancombc2-data/ThroatSwab-SEPPFT.qza"
)

tse_throat$time = factor(tse_throat$time, levels = c("D-7", "D3", "D5", "D7", "D10", "D14"))

output_genus_throat = ancombc2(data = tse_throat, assay_name = "counts", tax_level = "Genus",
                               fix_formula = "time",
                               rand_formula = "(1 | patient)",
                               p_adj_method = "BH", pseudo_sens = TRUE,
                               prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                               group = "time", struc_zero = TRUE, neg_lb = TRUE,
                               alpha = 0.05, n_cl = 2, verbose = TRUE,
                               global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                               iter_control = list(tol = 1e-2, max_iter = 20, 
                                                   verbose = TRUE),
                               em_control = list(tol = 1e-5, max_iter = 100),
                               lme_control = lme4::lmerControl(check.nobs.vs.nRE= "ignore"),
                               mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                               trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                           nrow = 2, 
                                                                           byrow = TRUE)),
                                                    node = list(2),
                                                    solver = "ECOS",
                                                    B = 100))

res_prim_throat <- output_genus_throat$res

# Nasal Swab
tse_nasal = mia::loadFromQIIME2(
  featureTableFile = "Data/ancombc2-data/NasalSwab-FT.qza",
  taxonomyTableFile = "Data/taxonomy_2.qza",
  sampleMetaFile = "Data/ancombc2-data/NasalSwab-MD.txt",
  refSeqFile = "Data/ancombc2-data/NasalSwab-RS.qza",
  phyTreeFile = "Data/ancombc2-data/NasalSwab-SEPPFT.qza"
)

tse_nasal$time = factor(tse_nasal$time, levels = c("D-7", "D3", "D5", "D7", "D10", "D14"))

output_genus_nasal = ancombc2(data = tse_nasal, assay_name = "counts", tax_level = "Genus",
                              fix_formula = "time",
                              rand_formula = "(1 | patient)",
                              p_adj_method = "BH", pseudo_sens = TRUE,
                              prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                              group = "time", struc_zero = TRUE, neg_lb = TRUE,
                              alpha = 0.05, n_cl = 2, verbose = TRUE,
                              global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                              iter_control = list(tol = 1e-2, max_iter = 20, 
                                                  verbose = TRUE),
                              em_control = list(tol = 1e-5, max_iter = 100),
                              lme_control = lme4::lmerControl(check.nobs.vs.nRE= "ignore"),
                              mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                              trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                          nrow = 2, 
                                                                          byrow = TRUE)),
                                                   node = list(2),
                                                   solver = "ECOS",
                                                   B = 100))


res_prim_nasal <- output_genus_nasal$res 


write.csv(output_genus_stool$res,"Outputs/Ancombc2/ancombc2-primary-results-stool.csv")
write.csv(output_genus_rectal$res,"Outputs/Ancombc2/ancombc2-primary-results-rectal.csv")
write.csv(output_genus_throat$res,"Outputs/Ancombc2/ancombc2-primary-results-throat.csv")
write.csv(output_genus_nasal$res,"Outputs/Ancombc2/ancombc2-primary-results-nasal.csv")
