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

Take a look at the results:

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

```sh
bcftools query -f '%ID\t%CHROM\t%POS\t[%TGT]\n' HG002.vcf.gz | head
```

```
rs6690515	  1	1447325	G/G
rs28635343	1	1560103	C/C
rs16824588	1	1900232	T/T
rs424079	  1	2071340	C/A
rs2055204	  1	2280661	G/G
rs4648384	  1	2756621	C/A
rs16823228	1	2853759	A/A
rs2817178	  1	3047252	C/C
rs731031	  1	3072692	T/T
rs2455107	  1	3180158	T/T
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

### Count comparison

Please see [build/comparison.md](build/comparison.md) for details on sites not captured by k2v.
