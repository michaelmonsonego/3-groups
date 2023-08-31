

Antigen_processing <- (c("H2-D1", "H2-Q7", "H2-K1", "H2-T23","H2-Q4","Tap1","Tapbp","B2m", "Psmb8", "Psme1",
                         "Psmb9", "Calr", "Psmb10", "Ncf1", "Fcer1g"))


IFNg_signature <- (c("Ccl9", "Cxcl9", "Cxcl10", "Cd274", "Gbp2", "Irf1", "H2-D1", "H2-Q4", "H2-Q7", "H2-K1", 
                     "Gbp7", "Stat1", "Irgm1", "B2m", "Irf9", "Ifitm3", "H2-T23", "Icam1", "Jak2", "Irf2", "Ifngr2"))


Nfkb_signature <- (c("Nfkb1", "Nfkbib", "Nfkbia", "Nfkbie", "Nfkb2", "Nfkbiz", "Rela", "Relb"))


############################################################################
# GOBP_RESPONSE_TO_INTERLEUKIN_2
INTERLEUKIN_2_RESPONSE <- gmtPathways('Genesets\\GOBP_RESPONSE_TO_INTERLEUKIN_2.gmt')[[1]]
INTERLEUKIN_2_RESPONSE <- str_to_title(INTERLEUKIN_2_RESPONSE)

############################################################################
# GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION
ANTIGEN_PROCESSING_AND_PRESENTATION <- gmtPathways('Genesets\\GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION.gmt')[[1]]
ANTIGEN_PROCESSING_AND_PRESENTATION <- str_to_title(ANTIGEN_PROCESSING_AND_PRESENTATION)

############################################################################
# GOMF_INTERLEUKIN_2_RECEPTOR_BINDING
INTERLEUKIN_2_RECEPTOR_BINDING <- gmtPathways('Genesets\\GOMF_INTERLEUKIN_2_RECEPTOR_BINDING.gmt')[[1]]
INTERLEUKIN_2_RECEPTOR_BINDING <- str_to_title(INTERLEUKIN_2_RECEPTOR_BINDING)

############################################################################
# GOBP_I_KAPPAB_KINASE_NF_KAPPAB_SIGNALING
NF_KAPPAB_SIGNALING <- gmtPathways('Genesets\\GOBP_I_KAPPAB_KINASE_NF_KAPPAB_SIGNALING.gmt')[[1]]
NF_KAPPAB_SIGNALING <- str_to_title(NF_KAPPAB_SIGNALING)

############################################################################
M1like <- (c('Ccl11','Ccl2','Ccl3','Ccl4','Ccl5','Ccl8','Cd80','Cd86','Cxcl1','Cxcl10',
             'Cxcl13','Cxcl16','Cxcl3','Cxcl9','Gbp5','Ifit2','Ifit3','Il18','Il1b','Il6',
             'Irf1','Irf5','Isg15','H2-DMb1','H2-Ab1','H2-Eb1','Nfkb1','Stat1','Stat5a',
             'Stat5b','Tlr2','Tlr4','Tnf','Il12rb1','Il23r'))
#'Ccl15','Ccl19','Ccl20','Cxcl11','Cxcl8','Gbp1'

M1like.2 <- c('IL23','TNF','CXCL9','CXCL10','CXCL11','CD86','IL1A','IL1B','IL6','CCL5','IRF5','IRF1','CD40',
              'IDO1' ,'KYNU','CCR7')
M1like.2 <- str_to_title(M1like.2)


M2like <- (c('Arg1','Cbr2','Ccl1','Ccl17','Ccl2','Ccl22','Ccl24','Ccl5','Ccna2',
             'Cd163','Cd69','Cybb','F13a1','Fn1','Folr2','Il10','Il13','Il33','Il4','Irf4',
             'Mrc1','S100a2','S100a9','Sdcbp','Stat6','Tap2','Tgfb1','Tk1','Tlr7'))


M2like.2 <- c('IL4R' ,'CCL4' ,'CCL13' ,'CCL20' ,'CCL17' ,'CCL18' ,'CCL22' ,'CCL24' ,'LYVE1' ,'VEGFA' ,'VEGFB' ,
              'VEGFC' ,'VEGFD' ,'EGF' ,'CTSA' ,'CTSB' ,'CTSC' ,'CTSD' ,'TGFB1' ,'TGFB2' ,'TGFB3' ,'MMP14' ,
              'MMP19' ,'MMP9' ,'CLEC7A' ,'WNT7B' ,'FASL' ,'TNFSF12' ,'TNFSF8' ,'CD276' ,'VTCN1' ,'MSR1' ,
              'FN1' ,'IRF4')
M2like.2 <- str_to_title(M2like.2)
#'Ccl13','Ccl14','Ccl18','Ccl23','Ccl26','Fcca1','Fizz1','Jnjd3','Ppar','Vsig4','Ym1','Ym2'

############################################################################
Tgfb_signature <- c('Abca12','Abcd2','Accsl','Acot12','Adad2','Adra1a','Agpat6','Aldh3b1','Alx1','Arhgdig','Arhgef37','Artn','Ass1','Bloc1s3','Bmp10','Bpifb2','Bub1b','C10orf81','C10orf90','C15orf59','C2orf42','C3','C6orf118','Capn8','Capns2','Ccnb3','Cd163','Cdh18','Cdk20','Cdkn2b','Cdt1','Cep128','Chrm3','Chrm4','Cldn10','Clec3b','Cnn3','Cobl','Cryaa','Cxcl17','Cylc2','Cyp2u1','Cyp4a11','Dclre1c','Dgkb','Dlg5','Dlk1','Dmrt2','Dmrt3','Dnah17','Dpcd','Drd2','Dtna','Dusp8','Ect2','Efcab11','Eif1ad','Epha6','Fabp4','Fam151a','Fam83g','Fbxl13','Fkbp2','Fkbp6','Foxj1','Frem2','Gcm1','Ggt6','Glb1l','Gp2','Gpr174','Grik3','Gsdmc','Gsta5','Hapln4','Hbz','Hcn1','Hrh1','Hspa8','Hspbap1','Ifitm10','Igdcc3','Igsf9','Ikzf1','Il31ra','Iqcf5','Iqsec3','Jph3','Kcnc2','Kcnj15','Kcnj6','Kcnq1','Kcnq4','Kif23','Kif4a','Klk11','Lactbl1','Lamc3','Lif','Lipi','Lipt2','Loc100507656','Lrp2','Mapk10','Mat1a','Mchr1','Mdh1b','Mep1b','Mettl11b','Mgat3','Mis18bp1','Morn3','Mpp1','Mpzl3','Msgn1','Myog','Mzb1','Ndufa4','Necab2','Nkx2-5','Nodal','Nos1','Nos3','Nr2e1','Odz2','Opn1lw','Oprl1','P2ry4','Pabpc5','Pamr1','Papd4','Parp1','Pcdhb6','Pdgfb','Pdilt','Phactr1','Phlda3','Pla2g4e','Plbd1','Plin5','Pnlip','Ppm1e','Pppde2','Prm3','Rab39','Rab3b','Rad51b','Rap2c','Rasgrp4','Rbp3','Reg3a','Rfx4','Rgr','Ribc2','Ripk2','Rnf182','Rtn4rl1','Saa1','Scn2a','Scn4b','Scnn1g','Sdk1','Sel1l3','Sgms2','Sh3rf3','Shox2','Slc1a1','Slc25a12','Slc30a3','Slc30a8','Slitrk5','Smad9','Sncb','Spats1','Sptbn2','Stac','Steap2','Syndig1','Tacr2','Tead1','Tex12','Tlr5','Tmc2','Tmem117','Tmem184c','Trpc5','Ttc18','Ttc22','Ttc36','Ubash3a','Ube2c','Ugt2b15','Vat1l','Vsx1','Wnt7b','Xrcc6','Zdhhc15','Znf541','Znf750','Zyg11a')

############################################################################
Tcells_Co_Inhibitory <- c('Pdcd1','Tigit','Havcr2','Ctla4','Lag3','Entpd1')

Tcells_Co_Stimulatory <- c('Tnfrsf4','Tnfrsf9','Tnfrsf18','Cd28','Icos')

Tcells_Cytotoxicity <- c('Gzma','Gzmb','Gzmk','Gzmh','Prf1','Infg','Tnfa','Il12a','Il12b','Il7r','Klrc1',
                         'Klrd1','Lilrb1','Nectin2')

Tregs_suppression <- c('Il10','Tgfb1','Tgfb3','Ctla4','Il2ra','Nt5e','Fgl2',
                       'Ebi3','Foxp3')
############################################################################
Myeloid_Co_Inhibitory <- c('Cd274','Pvr','Lgals9','Cd80','Adora2a','Adora2b')
Myeloid_Co_Stimulatory <- c('Tnfsf4','Tnfsf9','Tnfsf18','Cd80','Cd86','Icosl')

############################################################################




############################################################################
NKs_Immune_response <- c('Cd160','Cd226','Ceacam1','Crtam','Havcr2','Il12a','Il12a','Nectin2','Pvr')

############################################################################
NK_immature.vs.mature_DN <- gmtPathways('Genesets\\GSE13229_IMM_VS_MATURE_NKCELL_DN.gmt')[[1]]
NK_immature.vs.mature_DN <- str_to_title(NK_immature.vs.mature_DN)

NK_immature.vs.mature_UP <- gmtPathways('Genesets\\GSE13229_IMM_VS_MATURE_NKCELL_UP.gmt')[[1]]
NK_immature.vs.mature_UP <- str_to_title(NK_immature.vs.mature_UP)



############################################################################

Tcell_migration <- gmtPathways('Genesets\\GOBP_T_CELL_MIGRATION.gmt')[[1]]
Tcell_migration <- str_to_title(Tcell_migration)

Tcell_migration_regulation <- gmtPathways('Genesets\\GOBP_REGULATION_OF_T_CELL_MIGRATION.gmt')[[1]]
Tcell_migration_regulation <- str_to_title(Tcell_migration_regulation)

Tcell_chemotaxis <- gmtPathways('Genesets\\GOBP_T_CELL_CHEMOTAXIS.gmt')[[1]]
Tcell_chemotaxis <- str_to_title(Tcell_chemotaxis)

Tcell_chemotaxis_regulation <- gmtPathways('Genesets\\GOBP_REGULATION_OF_T_CELL_CHEMOTAXIS.gmt')[[1]]
Tcell_chemotaxis_regulation <- str_to_title(Tcell_chemotaxis_regulation)


############################################################################

Tcell_Integrin_binding <- gmtPathways('Genesets\\GOMF_INTEGRIN_BINDING.v7.5.1.gmt')[[1]]
Tcell_Integrin_binding <- str_to_title(Tcell_Integrin_binding)

Tcell_cell_adhesion <- gmtPathways('Genesets\\GOBP_CELL_CELL_ADHESION_MEDIATED_BY_INTEGRIN.v7.5.1.gmt')[[1]]
Tcell_cell_adhesion <- str_to_title(Tcell_cell_adhesion)

Tcell_Integrin_activation <- gmtPathways('Genesets\\GOBP_INTEGRIN_ACTIVATION.v7.5.1.gmt')[[1]]
Tcell_Integrin_activation <- str_to_title(Tcell_Integrin_activation)

Tcell_Integrin_signaling <- gmtPathways('Genesets\\GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY.v7.5.1.gmt')[[1]]
Tcell_Integrin_signaling <- str_to_title(Tcell_Integrin_signaling)


