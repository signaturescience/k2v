# k2v: Kintelligence reports to VCF

## Quickstart

Get some example data:

```sh
wget https://github.com/signaturescience/k2v/raw/main/testdata/HG002.SnpResult.txt
```

Pull the image from Docker Hub, and run the k2v container:

```sh
docker pull sigsci/k2v
docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) sigsci/k2v HG002.SnpResult.txt
```

## Output

Take a look at the VCF produced:

```sh
bcftools view HG002.vcf.gz | grep -v "^##" | head
```

```
#CHROM  POS     ID              REF     ALT     QUAL    FILTER  INFO    FORMAT  HG002
1       1447325 rs6690515       G       .       .       .       .       GT      0/0
1       1560103 rs28635343      C       .       .       .       .       GT      0/0
1       1900232 rs16824588      T       .       .       .       .       GT      0/0
1       2071340 rs424079        C       A       .       .       .       GT      0/1
1       2280661 rs2055204       G       .       .       .       .       GT      0/0
1       2756621 rs4648384       C       A       .       .       .       GT      0/1
1       2853759 rs16823228      A       .       .       .       .       GT      0/0
1       3047252 rs2817178       T       C       .       .       .       GT      1/1
1       3072692 rs731031        T       .       .       .       .       GT      0/0
```

Output data in a tabular format similar to 23andMe output:

```sh
bcftools query -f '%ID\t%CHROM\t%POS\t[%TGT]\n' HG002.vcf.gz | sed -s 's/[\/\|]//g' | head
```

```
rs6690515  1 1447325	GG
rs28635343 1 1560103	CC
rs16824588 1 1900232	TT
rs424079   1 2071340	CA
rs2055204  1 2280661	GG
rs4648384  1 2756621	CA
rs16823228 1 2853759	AA
rs2817178  1 3047252	CC
rs731031   1 3072692	TT
rs2455107  1 3180158	TT

```

Note that only 10,229 of the 10,230 sites on the Kintelligence are genotyped. `k2v` is currently limited to SNPs only, and cannot capture the one insertion that Kintelligence assays (rs796296176).

```sh
bcftools index --nrecords HG002.vcf.gz
```

```
10229
```

Number of sites genotyped on each chromosome:

```sh
bcftools index --stats HG002.vcf.gz
```


```
1	  249250621	811
2	  243199373	759
3	  198022430	663
4	  191154276	604
5	  180915260	600
6	  171115067	587
7	  159138663	529
8	  146364022	470
9	  141213431	462
10	135534747	501
11	135006516	469
12	133851895	496
13	115169878	375
14	107349540	351
15	102531392	344
16	90354753	376
17	81195210	351
18	78077248	351
19	59128983	274
20	63025520	307
21	48129895	178
22	51304566	180
X	  155270560	106
Y	  59373566	85
```

## Details

### Data location

Input to k2v is the `.SnpResult.txt` file produced by the Kintelligence assay. On the UAS server, navigate to the `Reports` folder within the relevant analysis folder for the sample of interest. E.g., `C:\verogen\fuas\runs\1\analyses\1\1\Reports` represents `C:\verogen\fuas\runs\<run-number>\analyses\<sample-number>\<analysis-number>\Reports`. The file is named `<sample-number>_S1_S.SnpResult.txt`.

### Singularity

If you're in an environment where you can't use Docker, use Singularity instead. Build a singularity image, and run the k2v container:

```sh
singularity build k2v.sif docker://sigsci/k2v:latest
singularity run k2v.sif HG002.SnpResult.txt
```

### Build

```sh
git clone https://github.com/signaturescience/k2v
cd k2v
docker build --no-cache -t k2v .
```
