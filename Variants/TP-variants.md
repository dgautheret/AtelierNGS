---
title: 'Lab session: Exome Sequencing Pipeline'
author: Daniel Gautheret - Universit<c3><a9> Paris-Sud, I2BC & Gustave Roussy Bioinformatics
  Core
date: "5/1/2019"
output:
  html_document:
    keep_md: true
  pdf_document: default
---



---
# Lab session: Exome Sequencing Pipeline
# Author: Daniel Gautheret - I2BC
# 4/1/2019
---

## 0/ Pre-requisites

This lab session requires:

- Knowledge of the Unix command line and simple scripts ("for" loops)
- Access to a Unix workstation or Virtual Machine (VM), hearafter called "VM" 
- Bioconda must be installed on the VM. all required software can be installed when needed using conda install.

## 1/ Data source

We will use the Open Access pancreatic cancer dataset from Baylor College of Medicine, Houston, published in this paper: 

- Becnel et al. (2015) An open access pilot freely sharing cancer genomic data from participants in Texas. Scientific Data, 3:160010. DOI: 10.1038/sdata.2016.10

Tumor and adjacent normal tissues from 7 patients were submitted to whole exome sequencing. Libraries were prepared from purified DNA and processed with WES capture kit NimbleGen SeqCap EZ Exome Library SR. Paired-end sequencing (2x100nt) was performed on Illumina HiSeq 2000, to a depth of about 2x60M reads (tumor & normal).

Here we'll work on data from patient ID TCRBOA7. To reduce file size, sequences mapping to Chromosome 16 have been extracted from the original BAM file and converted into fastq.gz files to simulate the actual output from a sequencing platform. 

## 2/ Downloading fastq files

Login into your Unix VM. Create a working directory that will host all files and subdirectories for this session. 

From within this directory, retrieve all fastq files using wget from URL:

https://transfert.u-psud.fr/om9m3ke/download

Unarchive then check file sizes and estimate read numbers for each sample.

## 3/ Creating BWA index

Create an index directory. Using wget, download into this directory the human chromosome 16 in fasta format from the UCSC web site:


```bash
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr16.fa.gz
```

Créez un dossier d'index, copiez le génome dans ce dossier, puis créez l'index BWA avec:

## 4/ Quality control + trimming 

Launch fastqc for an R1 & R2 file. Check quality drop at extremities. Use trimmomatic to eliminate low quality sequences at read extermities. There is no need to remove adapters with this dataset. 


```bash
trimmomatic PE <path/to/R1.gz> <path/to/R2.gz> -baseout <output-name.fastq>  LEADING:20 TRAILING:20 MINLEN:50
```

With these parameters, Trimmomatic removes bases from 5' et 3' ends as long as their quality is below 20. If the final read is shorter than 50, it is deleted. Four files are created: two R1+R2 files with trimmed reads (suffix 1P,2P) and two R1+R2 files with deleted reads (suffix 1U,2U). Argument -baseout specifies prefix.suffix of created files. Trimmomatic can work with .gz files. 

For the next 3 sections, we recommend to run the procedure for one single sample (R1+R2). When you are sure it works, create a shell script to run automatically the procedure for normal and tumor samples.


```bash

bwa index -a bwtsw <index-dir>/chr16.fa.gz

```

## 5/ Mapping with BWA

Create a directory for BWA output, then launch BWA with the following parameters:


```bash

bwa mem -M -t 2 -A 2 -E 1 <index-dir>/chr16.fa.gz <path/to/R1> <path/to/R2> > </path/to/outputsamfile> 

# Options used:
# -t INT        Number of threads 
# -A INT        Matching score. 
# -E INT Gap extension penalty. A gap of length k costs O + k*E (i.e. -O
#   is for opening a zero-length gap). 
# -M  Mark shorter split hits as secondary (for Picard compatibility).

```

Make sure the SAM file was created. It should be about 10 times bigger than the initial R1/R2 fastq.gz files.  

## 6/ Processing SAM files

- Use samtools view to convert SAM to BAM
- Use samtools sort to sort the BAM file.
- Use samtools index to index the BAM file
- Use samtools flagstats to check mapping stats


```bash
#Sam 2 Bam
samtools view -b </path/to/samfile> -o </path/to/bamfile> 

#Sort Bam
samtools sort </path/to/bamfile> -o </path/to/sortedbamfile> 

#Index bam file
samtools index </path/to/bamfile>

# flagstats
samtools flagstat </path/to/sorted-indexed-bamfile> > <path-to-flagstat-file>

```

Check the flagstat file. What fraction of reads were mapped? What fraction of pairs were properly mapped?

Now use samtools mpileup to convert each bam file to the pileup format. Beware: samtools mpileup requires a genome fasta file not in .gz format.  


```bash
#Convert to Mpileup
samtools mpileup -B -A -f <index-dir>/chr16.fa  <path/to/bam-file> > <path/to/mpileup-file>

```

## 8/ Shell Script for Normal + tumor

Now create a shell script to perform the previous analysis on normal and tumor samples. To loop over normal and tumor samples, you may use something like:


```bash
for F in normal tumor; do
  <trimmomatic> bla bla ... $F.fastq.gz ... bla bla
  <bwa> bla bla
  <samtools> bla bla
done
```

Your script should check whether directories exist and create them if necessary:


```bash
if [ ! -f BWAoutput ];then
  mkdir BWAoutput
fi
```

Do not include software installations and genome indexing in your main script, as it is something that is done only once. You may save all these preliminary commands in a separate file.

Your full pipeline should run as a single command line. Test it.  

## 9/ Calling somatic variants with Varscan

Now we use varscan to create variant lists from Normal and Tumor Pileup files:


```bash
varscan somatic <path/to/normal-pileup-file> <path/to/tumor-pileup-file> <path/to/vcf-file> \
  --variants --p-value 0.001 --min-avg-qual 15 --output-vcf 1
```

Check the VCF file is correct. Should contain somatic variants. (grep "SOMATIC").

## Extract somatic mutations
Using a combination of grep and the following awk command, extract all somatic mutations from VCF files (SNP and INDEL) and convert to BED format.


```bash
grep 'SOMATIC' <path/to/vcf-file> > <path/to/filtered-vcf>
awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' \
   <path/to/filtered-vcf> > <path/to/bed-file>
```

Using "bedtools intersect" and the Gencode GTF annotation (available in prof/Databases), extract all Gencode annotations for the somatic mutations.
Using a combination of grep and awk, extract only GTF lines corresponding to gene and the gene names:


```bash
bedtools intersect -a <path/to/genome-annot.gtf> -b <path/to/bed-file> > <path/to/interset-file> 
grep '\sgene\s' <path/to/interset-file> | awk '{print " " $1 " " $4 " " $5 " " $16}'
```

## 10/ Annotating mutations

CNV annotation is often performed with the Annovar software. However installation of all Annovar databases requires a lot of time and resources. As a simple alternative, you may use the Broad institute's Oncotator web site: https://portals.broadinstitute.org/oncotator/
Oncotator needs tab delimited file produced as follows:


```bash
awk '{print $1, "\t", $2, "\t", $2, "\t", $4, "\t",$5}' <path/to/vcf-file> > <path/to/output.tsv-file>
```

## 11/ Visualizing alignements with IGV

(Requires installing IGV on your local station. Then IGV needs Internet access to download reference genomes)

Download normal and tumor BAM files (+indexes) on local computer and view with IGV.
Locate and visualize a few somatic mutations identified by Varscan.

