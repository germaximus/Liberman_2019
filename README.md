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


<details><summary><b>Building necessary index files</b></summary>
  
```bash  
bowtie-build Elegans_rRNA.fa ./Elegans_indices/Elegans_rRNA  
```
</details>



### Sequencing reads filtering and mapping   
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
