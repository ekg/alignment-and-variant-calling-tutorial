# NGS alignment and variant calling

## Part 0: Setup

We're going to use a bunch of fun tools for working with genomic data:

1. [bwa](https://github.com/lh3/bwa)
2. [samtools](https://github.com/samtools/samtools)
3. [htslib](https://github.com/samtools/htslib)
4. [vt](https://github.com/atks/vt)
5. [freebayes](https://github.com/ekg/freebayes)
6. [vcflib](https://github.com/ekg/vcflib/)
7. [seqtk](https://github.com/lh3/seqtk)

In most cases, you can download and build these using this kind of pattern:

```bash
git clone https://github.com/lh3/bwa
cd bwa && make
```

or

```bash
git clone --recursive https://github.com/ekg/vcflib
cd vcflib && make
```

Let's assume you're in an environment where you've already got them available.

## Part 1: Alignment walk-through

[E. Coli K12](https://en.wikipedia.org/wiki/Escherichia_coli#Model_organism) is a common laboratory strain that has lost its ability to live in the human intestine, but is ideal for manipulation in a controlled setting.
The genome is relatively short, and so it's a good place to start learning about alignment and variant calling.

### E. Coli K12 reference

We'll get some test data to play with. First, [the E. Coli K12 reference](http://www.ncbi.nlm.nih.gov/nuccore/556503834), from NCBI. It's a bit of a pain to pull out of the web interface so [you can also download it here](http://hypervolu.me/~erik/genomes/E.coli_K12_MG1655.fa).

```bash
# the start of the genome, which is circular but must be represented linearly in FASTA
curl -s http://hypervolu.me/%7Eerik/genomes/E.coli_K12_MG1655.fa | head
# >NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
# AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
# TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA
# ...
```

### E. Coli K12 Illumina 2x300bp MiSeq sequencing results

For testing alignment, let's get some data from a [recently-submitted sequencing run on a K12 strain from the University of Exeter](http://www.ncbi.nlm.nih.gov/sra/?term=SRR1770413). We can us the [sratoolkit](https://github.com/ncbi/sratoolkit) to directly pull the sequence data (in paired FASTQ format) from the archive:

```bash
fastq-dump --split-files SRR1770413
```

`fastq-dump` is in the SRA toolkit. It allows directly downloading data from a particular sequencing run ID. SRA stores data in a particular compressed format (SRA!) that isn't directly compatible with any downsteam tools, so it's necessary to put things into [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) for further processing. The `--split-files` part of the command ensures we get two files, one for the first and second mate in each pair. We'll use them in this format when aligning.

```bash
# alternatively, you may want to first download, and then dump
# but this seems to fail sometimes for me
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR177/SRR1770413/SRR1770413.sra
sra-dump --split-files SRR1770413.sra
```

These appear to be paired 300bp reads from a modern MiSeq.

### E. Coli O104:H4 HiSeq 2000 2x100bp

As a point of comparison, let's also pick up a [sequencing data set from a different E. Coli strain](http://www.ncbi.nlm.nih.gov/sra/SRX095630[accn]). This one is [famous for its role in foodborne illness](https://en.wikipedia.org/wiki/Escherichia_coli_O104%3AH4#Infection) and is of medical interest.

```bash
fastq-dump --split-files SRR341549
```

## Setting up our reference indexes

### FASTA file index

First, we'll want to allow tools (such as our variant caller) to quickly access certain regions in the reference. This is done using the samtools `.fai` FASTA index format, which records the lengths of the various sequences in the reference and their offsets from the beginning of the file.

```bash
samtools faidx E.coli_K12_MG1655.fa
```

Now it's possible to quickly obtain any part of the E. Coli K12 reference sequence. For instance, we can get the 200bp from position 1000000 to 1000200. We'll use a special format to describe the target region `[chr]:[start]-[end]`.

```bash
samtools faidx E.coli_K12_MG1655.fa NC_000913.3:1000000-1000200
```

We get back a small FASTA-format file describing the region:

```text
>NC_000913.3:1000000-1000200
GTGTCAGCTTTCGTGGTGTGCAGCTGGCGTCAGATGACAACATGCTGCCAGACAGCCTGA
AAGGGTTTGCGCCTGTGGTGCGTGGTATCGCCAAAAGCAATGCCCAGATAACGATTAAGC
AAAATGGTTACACCATTTACCAAACTTATGTATCGCCTGGTGCTTTTGAAATTAGTGATC
TCTATTCCACGTCGTCGAGCG
```

### BWA's FM-index

BWA uses the [FM-index](https://en.wikipedia.org/wiki/FM-index), which a compressed full-text substring index based around the [Burrows-Wheeler transform](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform).
To use this index, we first need to build it:

```bash
bwa index E.coli_K12_MG1655.fa
```

You should see `bwa` generate some information about the build process.

```bash
[bwa_index] Pack FASTA... 0.04 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 2.26 seconds elapse.
[bwa_index] Update BWT... 0.04 sec
[bwa_index] Pack forward-only FASTA... 0.03 sec
[bwa_index] Construct SA from BWT and Occ... 0.72 sec
[main] Version: 0.7.8-r455
[main] CMD: bwa index E.coli_K12_MG1655.fa
[main] Real time: 3.204 sec; CPU: 3.121 sec
```

And, you should notice a new index file which has been made using the FASTA file name as prefix:

```bash
ls -rt1 E.coli_K12_MG1655.fa*
# -->
E.coli_K12_MG1655.fa
E.coli_K12_MG1655.fa.fai
E.coli_K12_MG1655.fa.bwt
E.coli_K12_MG1655.fa.pac
E.coli_K12_MG1655.fa.ann
E.coli_K12_MG1655.fa.amb
E.coli_K12_MG1655.fa.sa
```
