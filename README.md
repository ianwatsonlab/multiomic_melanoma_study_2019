# multiomic_melanoma_study_2019
Supplemental code and data for the paper entitled "Multi-omic analysis reveals new driver mutations and a sex-specific tumor suppressor in cutaneous melanoma"

### mutation_sex_bias_analysis.R
Reproduces the results presented in Fig. 2b

### mutation_signatures.R
Derives mutational signatures found in our study and reproduces results presented in Supplementary Figs. 2c and 3a.

### smg_alterations_vs_mRNA_groups.R
Reproduces results presented in Fig. 7a

### Somatic variants
An .RDS file containing somatic variants from the 5 melanoma studies can be accessed here: https://drive.google.com/file/d/1Ican79mU_eCva2Q86E3fh257ACAZutXN/view?usp=sharing

Variants are in hg38 coordinates and have been annotated with snpEff. Please refer to the method section of the publication for additional information about processing.
Note that protected variants from TCGA-SKCM are not included in this file, but they can be accessed through the GDC portal: https://portal.gdc.cancer.gov/
