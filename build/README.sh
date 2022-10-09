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
  | bcftools sort -Oz -o kintelligence.sites.vcf.gz
tabix -f kintelligence.sites.vcf.gz

# Create temporaty table that needs deduplication
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' kintelligence.sites.vcf.gz > tmp

# Deduplicate keeping the more common alt allele
Rscript dedupe.R && rm -f tmp

# Move to assets
cp kintelligence.alleles.csv ../assets
