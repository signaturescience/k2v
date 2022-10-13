## Kintelligence site VCF

# Create site TSV for extracting sites from GnomAD. The hg19_mapping.tsv file was supplied by Verogen.
cat hg19_mapping.tsv \
  | awk -v OFS='\t' 'NR>1 {print $3,$2,$1}' \
  | sed 's/^23/X/g' \
  | sed 's/^24/Y/g' \
  > kintelligence.sites.tsv
cut -f3 kintelligence.sites.tsv > kintelligence.rsids.txt
bgzip kintelligence.sites.tsv
tabix -f -s1 -b2 -e2 kintelligence.sites.tsv.gz
zcat kintelligence.sites.tsv | cut -f3 > kintelligence.rsids.txt

# Get the Gnomad GRCh37 site VCF (461GB - put in .gitignore)
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi

# Extract kintelligence sites from GnomAD VCF.
bcftools view -R kintelligence.sites.tsv.gz gnomad.genomes.r2.1.1.sites.vcf.bgz \
  | bcftools view -v snps -m2 -M2 \
  | bcftools view -i'ID=@kintelligence.rsids.txt' \
  | bcftools annotate -x ^INFO/AF \
  | bcftools sort -Oz -o kintelligence.autosomal.vcf.gz
tabix -f kintelligence.autosomal.vcf.gz

# Get GnomAD Y chromosome sites on GRCh38 and crossmap
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chrY.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chrY.vcf.bgz.tbi
bcftools annotate -x ^INFO/AF gnomad.genomes.v3.1.2.sites.chrY.vcf.bgz | bcftools view -m2 -M2 -v snps | bcftools sort -Oz -o tmp1.y.vcf.gz && tabix -f tmp1.y.vcf.gz
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
conda activate crossmap
CrossMap.py vcf --chromid s hg38ToHg19.over.chain.gz tmp1.y.vcf.gz ../assets/human_g1k_v37.fasta.gz tmp2.y.vcf
bcftools view -m2 -M2 -v snps tmp2.y.vcf | bcftools sort -Oz -o tmp3.y.vcf.gz && tabix -f tmp3.y.vcf.gz
bcftools view -R kintelligence.sites.tsv.gz tmp3.y.vcf.gz \
  | bcftools view -v snps -m2 -M2 \
  | bcftools view -i'ID=@kintelligence.rsids.txt' \
  | bcftools sort -Oz -o kintelligence.Y.vcf.gz
tabix -f kintelligence.Y.vcf.gz
bcftools index -s kintelligence.Y.vcf.gz
rm -f tmp*

bcftools concat kintelligence.autosomal.vcf.gz kintelligence.Y.vcf.gz | bcftools sort -Oz -o kintelligence.sites.vcf.gz
tabix -f kintelligence.sites.vcf.gz
bcftools index -s kintelligence.sites.vcf.gz
bcftools index -n kintelligence.sites.vcf.gz

# Create temporaty table that needs deduplication
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' kintelligence.sites.vcf.gz > tmp

# Deduplicate keeping the more common alt allele
Rscript dedupe.R && rm -f tmp

# Move to assets
cp kintelligence.alleles.csv ../assets

# # Get the dbSNP 155 VCF
# wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
# wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi

# # Reannotate chromosomes. Contents of GCF_000001405.25.gz-contigmap.tsv:
# # NC_000001.10    1
# # NC_000002.11    2
# # NC_000003.11    3
# # NC_000004.11    4
# # NC_000005.9     5
# # NC_000006.11    6
# # NC_000007.13    7
# # NC_000008.10    8
# # NC_000009.11    9
# # NC_000010.10    10
# # NC_000011.9     11
# # NC_000012.11    12
# # NC_000013.10    13
# # NC_000014.8     14
# # NC_000015.9     15
# # NC_000016.9     16
# # NC_000017.10    17
# # NC_000018.9     18
# # NC_000019.9     19
# # NC_000020.10    20
# # NC_000021.8     21
# # NC_000022.10    22
# # NC_000023.10    X
# # NC_000024.9     Y
# bcftools annotate --rename-chrs GCF_000001405.25.gz-contigmap.tsv -Oz -o GCF_000001405.25.reanno.vcf.gz GCF_000001405.25.gz && tabix -f GCF_000001405.25.gz
#
# # Extract kintelligence sites from dbSNP VCF.
# time bcftools view -R kintelligence.sites.tsv.gz GCF_000001405.25.reanno.vcf.gz \
#   | bcftools view -v snps -m2 -M2 \
#   | bcftools view -i'ID=@kintelligence.rsids.txt' \
#   | bcftools sort -Oz -o kintelligence.sites.vcf.gz
# tabix -f kintelligence.sites.vcf.gz
#
# # Create temporaty table that needs deduplication
# bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' kintelligence.sites.vcf.gz > tmp
#
# # Deduplicate keeping the more common alt allele
# Rscript dedupe.R && rm -f tmp
#
# # Move to assets
# cp kintelligence.alleles.csv ../assets



## HG004

wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed

bcftools merge -R kintelligence.sites.vcf.gz kintelligence.sites.vcf.gz HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz \
  | perl -pe 's/\t\.\t\.$/\tGT\t.\/./g' \
  | bcftools +setGT -- -t. -n0 \
  | bcftools view -T kintelligence.sites.vcf.gz \
  | bcftools annotate -x INFO,^FORMAT/GT \
  | bcftools view -T HG004_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed \
  | bcftools sort -Oz -o HG004.giab.vcf.gz
tabix -f HG004.giab.vcf.gz
bcftools index -s HG004.giab.vcf.gz
bcftools index -n HG004.giab.vcf.gz
mv HG004.giab.vcf.gz* ../testdata


## HG002

wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed

bcftools merge -R kintelligence.sites.vcf.gz kintelligence.sites.vcf.gz HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz \
  | perl -pe 's/\t\.\t\.$/\tGT\t.\/./g' \
  | bcftools +setGT -- -t. -n0 \
  | bcftools view -T kintelligence.sites.vcf.gz \
  | bcftools annotate -x INFO,^FORMAT/GT \
  | bcftools view -T HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed \
  | bcftools sort -Oz -o HG002.giab.vcf.gz
tabix -f HG002.giab.vcf.gz
bcftools index -s HG002.giab.vcf.gz
bcftools index -n HG002.giab.vcf.gz
mv HG002.giab.vcf.gz* ../testdata

