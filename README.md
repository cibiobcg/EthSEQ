# EthSEQ


Whole exome sequencing (WES) is widely utilized both in translational cancer genomics studies and in the setting of precision medicine. Stratification of individual’s ethnicity is fundamental for the correct interpretation of personal genomic variation impact. We implemented EthSEQ to provide reliable and rapid ethnicity annotation from whole exome sequencing individual’s data. EthSEQ can be integrated into any WES based processing pipeline and exploits multi-core capabilities.

EthSEQ requires genotype data at SNPs positions for a set of individuals with known ethnicity (the reference model) and either a list of BAM files or genotype data (in VCF format) of individuals with unknown ethnicity. EthSEQ annotates the ethnicity of each individual using an automated procedure and returns detailed information
about individual’s inferred ethnicity, including aggregated visual reports. 

***

## Perform ethnicity analysis with individuals genotype data from VCF file

Analysis of 6 individuals from 1,000 Genome Project using a reference model built from 1,000 Genome Project individual's genotype data. Genotype data for 10,000 SNPs included in Agilent Sure Select v2 captured regions are provided in input to EthSEQ in VCF format while reference model is provided in GDS format and describes genotype data for 1,000 Genome Project individuls for the same SNPs set. 

```{r}
library(EthSEQ)

## Run the analysis
ethseq.Analysis(
  target.vcf = system.file("extdata", "Samples_SS2_10000SNPs.vcf",
	package="EthSEQ"),
  out.dir = "/tmp/EthSEQ_Analysis/",
  model.gds = system.file("extdata", "Reference_SS2_10000SNPs.gds",
	package="EthSEQ"),
  verbose=TRUE,
  composite.model.call.rate = 1)

## Load and display computed ethnicity annotations
ethseq.annotations = read.delim("/tmp/EthSEQ_Analysis/Report.txt",
	sep="\t",as.is=TRUE,header=TRUE)
head(ethseq.annotations)

## Delete analysis folder
unlink("/tmp/EthSEQ_Analysis/",recursive=TRUE)
```
