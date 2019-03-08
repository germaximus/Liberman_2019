# Transcriptome and Ribosome profiling analysis

**Prerequisites:**  
[cutadapt 1.18](https://cutadapt.readthedocs.io/en/stable/index.html)  
[STAR-2.6.1d](https://github.com/alexdobin/STAR)  
[gffread utility](http://ccb.jhu.edu/software/stringtie/gff.shtml)  
Transcriptome samples were sequenced in paired-end 150 nt mode on Illumina sequencer.
Ribosome profiling samples were prepared with Illumina Small RNA TruSeq kit and sequenced in single-end 50 nt mode on Illumina sequencer.
Raw sequencing files are available from [GEO]().


### Preparing genome annotation and index files
C. elegans genomic sequences and annotation files (WS268) were downloaded from the [Wormbase](https://wormbase.org/).

| files                                       | MD5 check sum (unzipped)         | Description              |
| ------------------------------------------- |:--------------------------------:| -------------------------|
| c_elegans.PRJNA13758.WS268.annotations.gff3 | 2b353175bf6e8410815aede3a77a8a62 | annotation               |
| c_elegans.PRJNA13758.WS268.genomic.fa       | d570defcdc006a7c2859fc92dbb21bc4 | Genome sequence          |


<details><summary><b>Building necessary index files</b></summary>
  
```bash  
bowtie-build Elegans_rRNA.fa ./Elegans_indices/Elegans_rRNA  
```
</details>



### Ribo-seq Sequencing reads filtering and mapping   
<details><summary><b>Illumina adapter trimming.</b></summary>

```bash
cutadapt -j 20 -m 23 -a TGGAATTCTCGGGTGCCAAGG -o out.fastq input.fq.gz 
# -j      - number of threads
# -m      - discard reads shorter than 23 nucleotides after adapter trimming
```
</details>

<details><summary><b>Discard reads mapping to rRNA and PhiX</b></summary>
  
```bash
bowtie -p 36 --un filtered.fastq ./bowtie-1.2.1.1/Elegans_indices/Elegans_rRNA trimmed.fastq >/dev/null
```
</details>

<details><summary><b>Ribo-seq Sequencing Summary Statistics</b></summary>

|   sample   | total number of reads  |   rRNA + PhiX [%]  | footprints (first sequencing) | footprints (second sequencing) |
|:---------: |:----------------------:|:------------------:|:-----------------------------:|:------------------------------:|
|1WT20       |  29931496              |   25.51            |  22296445                     |  0                             |
|2WT20       |  26428417              |   41.50            |  15460677                     |  0                             |
|4WT20       |  22083711              |   23.88            |  16809424                     |  0                             |
|1WT37       |  18893098              |   91.37            |  1630864                      |  7563889                       |
|2WT37       |  27558119              |   50.91            |  13528926                     |  0                             |
|4WT37       |  18328652              |   50.20            |  9126816                      |  0                             |
|1CD20       |  31528904              |   29.90            |  22101589                     |  0                             |
|2CD20       |  26507787              |   46.05            |  14300563                     |  0                             |
|4CD20       |  18448868              |   25.54            |  13737481                     |  0                             |
|1CD37       |  15408038              |   60.19            |  6134583                      |  0                             |
|2CD37       |  21317358              |   59.35            |  8664500                      |  0                             |
|4CD37       |  16801168              |   57.87            |  7078945                      |  0                             |

<img src="Figures/RiboSeq_Summary_statistics1.png" width="600">

First round of sequencing revealed that sample 1WT37 yeilded low number of ribosomal footprints. Therefore, we re-sequenced this sample to increase the coverage.   

<img src="Figures/RiboSeq_Summary_statistics2.png" width="600">

</details>



