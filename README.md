# EthSEQ: ethnicity annotations from whole exome sequencing data

Whole exome sequencing (WES) is widely utilized both in translational cancer genomics studies and in the setting of precision medicine. Stratification of individual’s ethnicity is fundamental for the correct interpretation of personal genomic variation impact. We implemented EthSEQ to provide reliable and rapid ethnicity annotation from whole exome sequencing individual’s data. EthSEQ can be integrated into any WES based processing pipeline and exploits multi-core capabilities.

EthSEQ requires genotype data at SNPs positions for a set of individuals with known ethnicity (the reference model) and either a list of BAM files or genotype data (in VCF format) of individuals with unknown ethnicity. EthSEQ annotates the ethnicity of each individual using an automated procedure and returns detailed information
about individual’s inferred ethnicity, including aggregated visual reports. 

***

## Installation

You can either install EthSEQ from github repository using devtools package or directly from CRAN repository.

## Perform ethnicity analysis with individuals genotype data from VCF file

Analysis of 6 individuals from 1,000 Genome Project using a reference model built from 1,000 Genome Project individual's genotype data. Genotype data for 10,000 SNPs included in Agilent Sure Select v2 captured regions are provided in input to EthSEQ in VCF format while reference model is provided in GDS format and describes genotype data for 1,000 Genome Project individuls for the same SNPs set. 

```
library(EthSEQ)

out.dir = file.path(tempdir(),"EthSEQ_Analysis/")

## Run the analysis
ethseq.Analysis(
  target.vcf = system.file("extdata", "Samples_SS2_10000SNPs.vcf",
	package="EthSEQ"),
  out.dir = out.dir,
  model.gds = system.file("extdata", "Reference_SS2_10000SNPs.gds",
	package="EthSEQ"),
  verbose=TRUE,
  composite.model.call.rate = 1,
  space = "3D") # Default space is 2D

## Load and display computed ethnicity annotations
ethseq.annotations = read.delim(file.path(out.dir,"Report.txt"),
	sep="\t",as.is=TRUE,header=TRUE)
head(ethseq.annotations)

## Delete analysis folder
unlink(out.dir,recursive=TRUE)
```

Current version of EthSEQ manages only VCF files with the following format:
- FORMAT column should be only "GT"
- Only genotypes 0/0, 0/1, 1/1 and ./. are admitted
- Only positions with single reference and single alternative base are admitted
- No duplicate IDs are admitted (so no multiple variants with ID equal to ".") 
- No duplicated sample names are admitted
- No duplicated positions are admitted

## Perform ethnicity analysis using pre-computed reference model

Analysis of 6 individuals from 1,000 Genome Project using a reference model built from 1,000 Genome Project individual's genotype data. Genotype data for 123,292 SNPs included in Agilent Sure Select v2 captured regions are provided in input to EthSEQ in VCF format while reference model selected among the set of pre-computed reference model. Reference model SS2.Major refers to the reference model built from 1550 individuals (from AFR, EUR, SAS and EAS major populations) from 1,000 Genome Project and considering 123,292 SNPs included in Agilent Sure Select v2 captured regions. Note that a reference model version called SS2 considering gentoype data for more than 2,000 individuals from 1,000 Genome Project is also available.

```
library(EthSEQ)

data.dir = file.path(tempdir(),"EthSEQ_Data/")
out.dir = file.path(tempdir(),"EthSEQ_Analysis/")

## Download genotype data in VCF format
dir.create(data.dir)
download.file("https://github.com/aromanel/EthSEQ_Data/raw/master/Sample_SS2.vcf",
  destfile = file.path(data.dir,"Sample_SS2.vcf"))

## Run the analysis
ethseq.Analysis(
  target.vcf =  file.path(data.dir,"Sample_SS2.vcf"),
  out.dir = out.dir,
  model.available = "SS2.Major",
  model.folder = data.dir,
  verbose=TRUE,
  composite.model.call.rate = 1,
  space = "3D") # Default space is 2D

## Delete analysis folder
unlink(data.dir,recursive=TRUE)
unlink(out.dir,recursive=TRUE)
```

## Perform ethnicity analysis from BAM files list

Analysis of individual NA07357 from 1,000 Genome Project using a reference model built from 1,000 Genome Project individual's genotype data. Genotype data for 10,000 SNPs included in Agilent Sure Select v2 captured regions are provided in input to EthSEQ with a BAM file. reference model is provided in GDS format and describes genotype data for 1,000 Genome Project individuls for the same SNPs set. Note than the BAM given in input to EthSEQ is a toy BAM file containing only reads overlapping the positions of the 10,000 SNPs considered in the analysis.

```
library(EthSEQ)

data.dir = file.path(tempdir(),"EthSEQ_Data")
out.dir = file.path(tempdir(),"EthSEQ_Analysis")

## Download BAM file used in the analysis
dir.create(data.dir)
download.file(
 "https://github.com/aromanel/EthSEQ_Data/raw/master/NA07357_only10000SNPs.bam",
 destfile = file.path(data.dir,"Sample.bam"))
download.file(
 "https://github.com/aromanel/EthSEQ_Data/raw/master/NA07357_only10000SNPs.bam.bai",
 destfile = file.path(data.dir,"Sample.bam.bai"))

## Create BAM files list 
write(file.path(data.dir,"Sample.bam"),file.path(data.dir,"BAMs_List.txt"))

## Run the analysis
ethseq.Analysis(
  bam.list = file.path(data.dir,"BAMs_List.txt"),
  out.dir = out.dir,
  model.gds = system.file("extdata","Reference_SS2_10000SNPs.gds",
     package="EthSEQ"),
  verbose=TRUE,
  aseq.path = out.dir,
  mbq=20,
  mrq=20,
  mdc=10,
  run.genotype = TRUE,
  composite.model.call.rate = 1,
  cores=1,
  bam.chr.encoding = FALSE) # chromosome names encoded without "chr" prefix in BAM files

## Delete analysis folder
unlink(data.dir,recursive=TRUE)
unlink(out.dir,recursive=TRUE)
```

## Perform ethnicity analysis using multi-step refinement

Multi-step refinement Analysis of individuals from 1,000 Genome Project using a reference model built from 1,000 Genome Project individual's genotype data (set of analysed individuals and individuals used for the reference model are disjoint). Genotype data for 10,000 SNPs included in Agilent Sure Select v2 captured regions are provided in input to EthSEQ in GDS format while reference model is provided in GDS format and describes genotype data for 1,000 Genome Project reference individuls for the same SNPs set. Multi-step refinement tree is constructed as matrix. Non-empty cells in columns i contains parent nodes for non-empty cells in columns i+1. Ethnic groups in child nodes should be included in parent nodes, while siblings node ethnic groups should be disjoint. Consult EthSEQ paper supplementary material for more complicated examples.  

```
library(EthSEQ)

out.dir = file.path(tempdir(),"EthSEQ_Analysis")
data.dir = file.path(tempdir(),"EthSEQ_Data")

## Download genotype data in VCF format
dir.create(data.dir)
download.file("https://github.com/aromanel/EthSEQ_Data/raw/master/Target_SS2_10000SNPs.gds",
  destfile = file.path(data.dir,"Target_SS2_10000SNPs.gds"))

## Create multi-step refinement matrix
m = matrix("",ncol=2,nrow=2)
m[1,1] = "SAS|EUR|EAS"
m[2,2] = "SAS|EUR"

## Run the analysis on a toy example with only 10000 SNPs
ethseq.Analysis(
  target.gds = file.path(data.dir,"Target_SS2_10000SNPs.gds"),
  out.dir = out.dir,
  model.gds = system.file("extdata","Reference_SS2_10000SNPs.gds",
	package="EthSEQ"),
  verbose=TRUE,
  composite.model.call.rate = 1,
  refinement.analysis = m,
  space="3D")

## Delete analysis folder
unlink(out.dir,recursive=TRUE)
```

## Create a reference model from multiple VCF genotype data files

Construction of a reference model from two genotype data files in VCF format and a corresponding annotation files which described ethnicity and sex of each sample contained in the genotype data files.

```
library(EthSEQ)

out.dir = tempdir()
dir.create(out.dir)

### Load list of VCF files paths
vcf.files = 
  c(system.file("extdata", "VCF_Test_1.vcf", package="EthSEQ"),
    system.file("extdata", "VCF_Test_2.vcf", package="EthSEQ"))

### Load samples annotations
annot.samples = read.delim(system.file("extdata", "Annotations_Test.txt",
	package="EthSEQ"))

### Create reference model
ethseq.RM(
  vcf.fn = vcf.files,
  annotations = annot.samples,
  out.dir = out.dir,
  model.name = "Reference.Model",
  bed.fn = NA,
  call.rate = 1,
  cores = 1)

## Delete example file
unlink(out.dir,recursive=TRUE)
```

## Reference paper

EthSEQ: ethnicity annotation from whole exome sequencing data.

Alessandro Romanel#, Tuo Zhang, Olivier Elemento, Francesca Demichelis#.

Bioinformatics, 33(15):2402-2404, 2017.

## Additional papers

Comprehensive Analysis of Genetic Ancestry and Its Molecular Correlates in Cancer.

Jian Carrot-Zhang*, Nyasha Chambwe*, Jeffrey S. Damrauer*, Theo A. Knijnenburg*, A. Gordon Robertson*, Christina Yau*, Wanding Zhou*, Ashton C. Berger*, Kuan-lin Huang*, Justin Y. Newberg*, R. Jay Mashl**, Alessandro Romanel**,  Rosalyn W. Sayaman**, Francesca Demichelis, Ina Felau, Garrett M. Frampton, Seunghun Han, Katherine A. Hoadley, Anab Kemal, Peter W. Laird, Alexander J. Lazar, Xiuning Le, Ninad Oak, Hui Shen, Christopher K. Wong, Jean C. Zenklusen, Elad Ziv, Cancer Genome Atlas Analysis Network, Andrew D. Cherniack#, Rameen Beroukhim#. 

Cancer Cell, 37(5):639-654, 2020.

