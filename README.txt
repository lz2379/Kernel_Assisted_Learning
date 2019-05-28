the original data: iwpc_data_7_3_09_revised.csv
is public available:  It can be downloaded from
https://www.pharmgkb.org/downloads/ and it is the "IWPC Data"

 the .r file is the data trimming file. The .csv is the
sheet 2 (without comments) of the .xls file. I also attach the NEJM
paper and supp pdf. I am not sure if I trimmed the data in the best
way, so please use it/modify at you will



The covariates I used for constructing my training data are similar to
the S1e Pharmacogenetic dosing algorithm, see page 5 of the NEJM
supplementary pdf  (while I don't use the missing indicator as I use
the complete data only).

I did not impute covariates (except for one SNP in VKORC1, see page 12
of the supp) such that I discard the incomplete data. I agree the
dataset is confusing for the genotype, you may want to read the
meta_data (sheet 1 of the .xls) to get definition of the covariate. In
addition, although there are many SNPs in VKORC1, the S1e only uses
rs9923231, and I think the genotype for CYP2C9 is the variable named
CYP2C9_consensus.



For the Enzyme and Amiodarone (other drugs that can have interaction
with Warfarin), I was quite restrictive of the ascertainment accuracy
such that I discarded the missing data (I keep only the records with
drug status = 0 or 1). This can have a large impact on sample size. We
ended up with 1700+ samples, but it is possible that the NEJM paper
may have let the NA to be 0 when training their algorithm (the paper
did not mention how the missing data was handled). I think if we let
NA to be 0 for drug status, then the sample size (complete samples)
would be about 3500+.  I also discard dose level seems to be outliers
(<6, or >95), which one can argue different cut-off point should be
used.