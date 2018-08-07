About
==========
Tn5 transposase is now utilized for fragmentation of the DNA in ATAC-seq (Buenrostro et al., 2013), scRNA-seq (Islam et al., 2014) and WGBS (Wang et al, 2013). Compared with other fragmentases and physical methods, the Tn5 is efficient for generating sequencing libraries by fragmenting the DNA and appending a 19-bp ME sequence and adapters in a single step.  This program is designed to identify ME sequence in the paired-end reads and trim the 'polluted' sequences. I recommend you to do this step before mapping reads since it will improve the performance of alignment and probably influence the downstream analysis based on insert size.

Install
=====
```
## clone the updated program from GitHub
git clone https://github.com/shiquan/TN5dyncut
## get into the source codes directory
cd TN5dyncut
## build the source code
make
```

Usage
=====
```
## bwa-mem
dyncut -t 4 -m 2 -tail 3 -d SRR891269_1.fastq.gz SRR891269_2.fastq.gz | bwa mem -pt 12 ref.fa | samtools view -Sb -o aln.bam

## bowtie
dyncut -t 4 -m 2 -tail 3 -d SRR891269_1.fastq.gz SRR891269_2.fastq.gz -1 read1.fq.gz -2 read2.fq.gz
bowtie2 -x ref.fa -1 read1.fq.gz -2 read2.fq.gz | samtools view -Sb -o aln.bam
```

How to cite
====
This program is free to use, distribute and modify in research and commercial purposes.
Please cite https://github.com/shiquan/TN5dyncut if you use it.

