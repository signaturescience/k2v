
## Kintelligence / k2v SNP count comparison

k2v joins the `<sample>.SnpResult.txt` to an allele table to get
genotypes which are used to convert to VCF.

Building the k2v allele table starts with the
[`hg19_mapping.tsv`](hg19_mapping.tsv) file provided by Verogen, which
contains RSIDs, chromosome, and position information for 10,380 SNPs,
including 90 SNPs on the Y chromosome. 10,230 of these SNPs are
genotyped and included in the `<sample>.SnpResult.txt` output, including
85 SNPs on Y chromosome.

The allele table is created by extracting biallelic SNPs from GnomAD
genomes for 1-22 and X using the GRCh37 GnomAD site VCF. Y chromosome
sites from GnomAD are not available on GRCh37, so the GRCh38 Y
chromosome site VCF from GnomAD is lifted over to GRCh37 prior to
extracting Y chromosome Kintelligence sites. Reference and alternate
alleles are extracted from the GnomAD site VCF(s) and used for
converting 0/0, 0/1, 1/1 genotypes in the `<sample>.SnpResult.txt`.

**48** variants are not captured by k2v. These include 47 Y chromosome
SNPs, and 1 non-Y variant. This variant is rs796296176 (called “N29insA”
in the Kintelligence output): a reference C to CA insertion at
16:89985750. This insertion cannot be captured by k2v. The 47 Y SNPs not
captured by k2v are shown in the table below.

| chr |      pos | rsid       |
|:----|---------:|:-----------|
| Y   |  2661306 | rs35840667 |
| Y   |  6854765 | rs7067290  |
| Y   |  6952061 | rs13303791 |
| Y   |  7035520 | rs3109537  |
| Y   |  7134120 | rs3116341  |
| Y   |  7570678 | rs13304869 |
| Y   |  7571122 | rs13303599 |
| Y   |  8306262 | rs13304587 |
| Y   |  9060781 | rs34064949 |
| Y   |  9098400 | rs35180016 |
| Y   |  9417393 | rs4618902  |
| Y   | 14001232 | rs17316007 |
| Y   | 14028148 | rs35617575 |
| Y   | 14074877 | rs13304612 |
| Y   | 14096577 | rs9785941  |
| Y   | 14814540 | rs13304396 |
| Y   | 14864191 | rs7893052  |
| Y   | 15015396 | rs13304992 |
| Y   | 15042327 | rs13305209 |
| Y   | 15156388 | rs13305939 |
| Y   | 15350616 | rs13304350 |
| Y   | 15401877 | rs28651585 |
| Y   | 15526695 | rs2032664  |
| Y   | 15585558 | rs13304733 |
| Y   | 15816113 | rs13305905 |
| Y   | 16767517 | rs13304414 |
| Y   | 16859582 | rs3898902  |
| Y   | 16898428 | rs13305838 |
| Y   | 16903692 | rs9724561  |
| Y   | 16904849 | rs13305362 |
| Y   | 17594966 | rs9785996  |
| Y   | 17805633 | rs13305443 |
| Y   | 18099341 | rs13304877 |
| Y   | 18871216 | rs13303804 |
| Y   | 19037067 | rs13305519 |
| Y   | 19259569 | rs7067448  |
| Y   | 21157531 | rs7067330  |
| Y   | 21402296 | rs35974449 |
| Y   | 22047290 | rs9943029  |
| Y   | 23059591 | rs9785991  |
| Y   | 23148103 | rs34053023 |
| Y   | 23342709 | rs9785894  |
| Y   | 23655112 | rs9786064  |
| Y   | 23759180 | rs13305670 |
| Y   | 23959349 | rs4144073  |
| Y   | 28650343 | rs9786795  |
| Y   | 28727063 | rs7474433  |

Y chromosome SNPs not captured by k2v. These SNPs are not in the GnomAD
GRCh38 site VCF after liftover to GRCh37.
