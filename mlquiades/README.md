mlquiades

Horvath Lab + Joe Goldfrank

2025 07 11

This package takes in bulk RNA cancer cell line sequencing data and GDSC1 IC50 drug sensitivity scores for palbociclib to build and evaluate 6 machine learning models.
These models include: decision tree, gradient boosted decision tree, neural net (not stable), neural net with hyperband, random forest, and ridge classifier.
There are three options for feature, in this case gene, selection. They include: only CDK4 and CDK6 related genes (the target for palbociclib); only CDK4, CDK6 and cancer genes (COSMIC); and a Pearson correlation method that keeps only the genes that have >=.3 rho value with the IC50 score in the training dataset only.


// example of how to run on cmd

python src/model_ml_drug_sensitivity.py --a sample_data --b output --c CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz --d palbociclib.csv --e 4 --o 3 --s cdk4_6_genes --t cdk4_6_genes.txt --u gencode.v19.genes.v7_model.patched_contigs.gtf.gz --r cdk_4_6_genes
