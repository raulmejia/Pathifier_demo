Rscript UptodateKEGG.R ../../Results/
tail -25 ../../Results/KEGG_pathways_in_df_genesymbol.tsv > ../../Results/KEGG_pathways_in_df_genesymbol_demo.tsv
Rscript Pathifier_modified.R TCGA_Control_vs_Basal_indicator_10and10.tsv ../../Results/KEGG_pathways_in_df_genesymbol_demo.tsv ./ ../../Results/Pathifier/ Basal 4 3.75
