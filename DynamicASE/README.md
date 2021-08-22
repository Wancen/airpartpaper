# Bulk RNA-seq dataset

* [dynamicASE](https://github.com/Wancen/airpartpaper/blob/main/DynamicASE/dynamicASE.R): the script performing *airpart* analysis with adjusted individual as a batch effect on  top 40 genes with largest variance of weighted allelic ratio mean. 

* [extract_gene_from_snp](https://github.com/Wancen/airpartpaper/blob/main/DynamicASE/extract_gene_from_snp.R): Preprocessing script. Since signals from the same gene are tightly correlated and it would disrupt the inference if we dealt with each SNP independently, so we chose the SNP with largest total counts per gene.

