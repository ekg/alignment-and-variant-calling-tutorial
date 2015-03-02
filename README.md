# NGS alignment and variant calling

## Part 0: Setup

We're going to use a bunch of fun tools for working with genomic data:

1. [bwa](https://github.com/lh3/bwa)
2. [samtools](https://github.com/samtools/samtools)
3. [htslib](https://github.com/samtools/htslib)
4. [vt](https://github.com/atks/vt)
5. [freebayes](https://github.com/ekg/freebayes)
6. [vcflib](https://github.com/ekg/vcflib/)

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

### E. Coli data

We'll get some test data to play with. First, [the E. Coli K12 reference](http://www.ncbi.nlm.nih.gov/nuccore/556503834), from NCBI. It's a bit of a pain to pull out of the web interface so [you can also download it here](http://hypervolu.me/~erik/genomes/E.coli_K12_MG1655.fa).

```bash
# the start of the genome, which is circular but must be represented linearly in FASTA
curl -s http://hypervolu.me/%7Eerik/genomes/E.coli_K12_MG1655.fa | head
# >NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
# AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
# TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA
# ...
``` bash

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

