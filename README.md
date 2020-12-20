# Taiga_Metagenome_Project

This project is a training project in Bioinformatic institute. In this case, we will try to describe two different types of soils - black taiga soil and regional common soil. We analyse Oxford Nanopore reads and Illumina 16s amplicones and try to extract some useful information about communities of microorganisms.

## 1. Aims and Tasks

Current task was to decribe composition of soil microbiome in different soils, in particular:

* perform contigs binning and try to obtain MAGs;
* extract 16s rDNA from contigs, compare with amplicone sequences;
* annotate assembled metagenomes
* extract information about glycoside hydrolases

## 2. Sources and Methods

In current work we compared two spots,

* N1 - black taiga soil samples
* N2 - regional blank soil (control)

Sampling, DNA extraction and sequencing was performed previously, raw ONP reads was assembled by Flye and described by QUAST

In this analysis we worked with assembled metagenomes (called as N1 and N2), and 16s amplicone sequences (s6 - s10, s26-s30)

## 3. System Requement

We perform followed analysis in server computer (20 cores, 500Gb memory, Ubintu 16 LTS), but with several limits it can be performed in a PC up to 4 cores and 8Gb memory. Also we use followed software:

* [Anaconda](https://www.anaconda.com/) 4.6.11
* [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/) 2:2.15
* [CheckM](https://github.com/Ecogenomics/CheckM/wiki) 1.0.18
* [barrnap](https://github.com/tseemann/barrnap) 0.9
* [QIIME](http://qiime.org/) 1.9.1
* [PROKKA](https://github.com/tseemann/prokka) 1.14.6
* [R](https://www.r-project.org/) 4.0.3 and libraries:
    * Biostrings 2.54.0
    * dada2 1.14.1
    * dplyr 1.0.2
    * plyr 1.8.6
    * ggplot2 3.3.2
    * phyloseq 1.30
* [SILVA](https://www.arb-silva.de/) 132 database


## 4. Processing

### 4.1 Quality information

Information about assembly obtained from QUAST software. All required information extracted from info files, provided alonside assemples. Here is information about assembles:

|   | N1 | N2 |
|---|----|----|
|N50   | 61 838 | 34 361 |
|contigs >1000  |15k|15k|
|largest contig   | 3.6Mbp | 1.2Mbp|

We can see, that our N1 assembly is better in contigs' length, but not in contigs' numbers. Large contig info also interesting - as we can remember *E.coli* genome size as 4.6Mbp, we can expect almost-full-MAG in one contig.

## 4.2 Binning

We perform binning using `MetaBAT2` by followed command:

```
$metabat2 -i source/N1_N2_taiga_nanopore/N1/consensus.fasta -o binned2/N1/ --minContig 2000 --numThreads 20
    50 bins (617270216 bases in total) formed.

$metabat2 -i source/N1_N2_taiga_nanopore/N2/consensus.fasta -o binned2/N2/ --minContig 2000 --numThreads 20
    17 bins (316161957 bases in total) formed.
```

In this case, `--minContig` is minimum length for contig to be binned, and `--threads` is number o cores, used for calculation.

Quite expected, we find more bins in assembly with more elongate contigs. Let's see, how correct they are.

## 4.3 Check quality of bins

We check quality of bins using `checkM` program. Commands were:

```
$checkm lineage_wf -x fa --tab_table -f checkm/N1_MAGs_checkm.tab -t 20 binned2/N1/ checkm/N1
$checkm bin_qa_plot checkm/N1/ binned2/N1 checkm/plots_N1 -x fa

$checkm lineage_wf -x fa --tab_table -f checkm/N2_MAGs_checkm.tab -t 20 binned2/N2/ checkm/N2
$checkm bin_qa_plot checkm/N2/ binned2/N2 checkm/plots_N2 -x fa
```

Flags are `-x` is a file format for bins (in our case is .fa); `-f` for printing results in file; `--tab-table` for tab-separated values table. Commands were `lineage_wf` for lineage, analyze and quality analysis, and `bin_qa_plot` for plotting quality data

N1:
![N1](N1_bin_qa_plot.jpg)

N2:
![N2](N2_bin_qa_plot.jpg)

So, as we can see, our bins are overcontaminated, and we can't use them.


## 4.5 Extract rRNA gene

For rRNA genes extraction, we use Barrnap by this command:
```
$barrnap source/N1/consensus.fasta >16s_rdna/rrna_N1.gff -o 16s_rdna/rrna_N1.fna --threads 20 --lencutoff 95
$grep -A 1 ">16S" 16s_rdna/rrna_N1.fna >16s_rdna/16s_rdna_N1.fna
$grep ">" 16s_rdna/16s_rdna_N1.fna | wc -l
    95

$barrnap source/N2/consensus.fasta >16s_rdna/rrna_N2.gff -o 16s_rdna/rrna_N2.fna --threads 20 --lencutoff 95
$grep -A 1 ">16S" 16s_rdna/rrna_N2.fna >16s_rdna/16s_rdna_N2.fna
$grep ">" 16s_rdna/16s_rdna_N2.fna | wc -l
    44
```

In this case, flags `--treads` for number of cores, `--lencutoff` for percent of match for RNA gene sequence. Flag fod `grep` `-A` determain, how many strings should be shown after match.

Here we extract 95 and 44 16s rRNA gene sequences from metagenomes N1 and N2 respectively.

## 4.6 Annotation if 16s rRNA gene sequences

Make an OTU table for extracted 16s rRNA genes, using QIIME

```
$pick_de_novo_otus.py -r /home/alexey/tax_n_refs/silva_132_97_16S.fna -t /home/alexey/tax_n_refs/taxonomy_7_level.txt -i all.fasta -o otus
```
Downstream analysis performed in R. All code from this step you can find in `processing.R` file

Here we can show taxonomical structure in phyla level for amplicone and metagenomic data

16s rRNA gene from amplicons:
![amplicone](amplicons.jpg)

16s rDNA from metagenomes:
![metagenome](metagenome.jpg)

## 4.7 PROKKA annotations

Assembles were annotate by PROKKA:

```
$prokka source/N1/consensus.fasta --outdir prokka/N1/ --metagenome --cpus 20
$prokka source/N2/consensus.fasta --outdir prokka/N2/ --metagenome --cpus 20

$cat N1/PROKKA_12102020.gff | grep "product=hypothetical" | wc -l
$cat N1/PROKKA_12102020.gff | grep "product=" | wc -l
$cat N1/PROKKA_12102020.gff | grep "eC_number=3.2.1" | wc -l

$cat N2/PROKKA_12102020.gff | grep "product=hypothetical" | wc -l
$cat N2/PROKKA_12102020.gff | grep "product=" | wc -l
$cat N2/PROKKA_12102020.gff | grep "eC_number=3.2.1" | wc -l
```
Results are in followed table:

|   | N1 | N2 |
|---|----|----|
|Assigned protein   | 176k | 89k |
|Hypotetical protein  |750k | 402k|
| Glycoside Hydrolases | 1238 | 629 |

## 5. Results

In our work, we performed basic analysis of metagenomic assembles, and compare with amplicone sequences. 
1. Our assembles have same number of contigs, but different length of them. By this reason, we can't compare number of genes, we should normalise them first
2. Binning via MetaBAT2 is not sucsessful due high contamination level. It can be linked to high diversity of soils, and high variety of matrix. But, perhaps, we should use another binning package (for example, Maxbin)
3. 16s rRNA genes, extracted from metagenome assembles, all are discovered in amplicone reads on genus level. In common view, on Phyla level, taxonomical structure of amplicone and metagenome communities are similar
4. PROKKA allowed to annotate hypotetical and annotated proteins (and glycoside hydrolase genes in particular), but their numbers should be corrected









