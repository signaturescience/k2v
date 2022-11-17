#!/bin/sh

# Usage: k2v.sh <KintelligenceInput.SnpResult.txt>

# Get input filename, get sample ID and create output filename
INFILE=$1
OUTDIR=$(dirname $INFILE)
SAMPLEID=$(basename $INFILE | sed "s/\.[Ss]np[Rr]esult\.txt$//g" | sed "s/\.txt$//g")
OUTVCF=${OUTDIR}/${SAMPLEID}.vcf.gz

if [[ -f "$INFILE" ]]; then
  echo
  echo "Sample:    $SAMPLEID"
  echo "Input:     $INFILE"
  echo "Output:    $OUTVCF"
  echo
else
  echo "File not found: $INFILE"
  exit 1
fi

# Convert to tabular format for bcftools convert
Rscript /scripts/k2tsv.R $INFILE

# Convert to VCF, sort, index
{ bcftools convert --tsv2vcf /tmp/tmp.23.txt --fasta-ref /assets/human_g1k_v37.fasta.gz --samples $SAMPLEID | bcftools sort -Oz -o $OUTVCF && tabix -f $OUTVCF; } > /dev/null 2>&1
