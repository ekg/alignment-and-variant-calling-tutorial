# Alignment and variant calling

This tutorial steps through some basic tasks in alignment and variant calling using a handful of Illumina sequencing data sets. For theoretical background, please refer to the included [presentation on alignment and variant calling](https://docs.google.com/presentation/d/1t921ccF66N0_oyn09gbM0w8nzADzWF20rfZkeMv3Sy8/edit?usp=sharing), or the [included PDF from a previous year](https://github.com/ekg/alignment-and-variant-calling-tutorial/raw/master/presentations/Alignment%2C%20Variant%20Calling%2C%20and%20Filtering%20(WGC%20NGS%20Bioinformatics).pdf).
    
## Part 0: Setup

We're going to use a bunch of fun tools for working with genomic data:

1. [bwa](https://github.com/lh3/bwa)
2. [samtools](https://github.com/samtools/samtools)
3. [htslib](https://github.com/samtools/htslib)
4. [vt](https://github.com/atks/vt)
5. [freebayes](https://github.com/ekg/freebayes)
6. [vcflib](https://github.com/ekg/vcflib/)
7. [sambamba](https://github.com/lomereiter/sambamba)
8. [seqtk](https://github.com/lh3/seqtk)
9. [sra-tools](https://github.com/ncbi/sra-tools/wiki/HowTo:-Binary-Installation)
10. [rtg-tools](https://www.realtimegenomics.com/products/rtg-tools)

In most cases, you can download and build these using this kind of pattern:

```bash
git clone https://github.com/lh3/bwa
cd bwa && make
```

or, in the case of several packages (vcflib, sambamba, freebayes, and vg), submodules are used to control the dependencies of the project, and so the whole source tree must be cloned using the `--recursive` flag to git. For example, here is how we'd clone and build freebayes:

```bash
git clone --recursive https://github.com/ekg/freebayes
cd freebayes && make
```

In some cases you can download precompiled binaries. For instance, you can head
to the [sambamba releases page](https://github.com/lomereiter/sambamba/releases) to find binaries that
should work on any modern Linux or OSX distribution. The same applies to the [sra-toolkit, which is probably easier to install from available binaries](https://github.com/ncbi/sra-tools/wiki/HowTo:-Binary-Installation).

Otherwise, let's assume you're in an environment where you've already got them available.

## Part 1: Aligning E. Coli data with `bwa mem`

[E. Coli K12](https://en.wikipedia.org/wiki/Escherichia_coli#Model_organism) is a common laboratory strain that has lost its ability to live in the human intestine, but is ideal for manipulation in a controlled setting.
The genome is relatively short, and so it's a good place to start learning about alignment and variant calling.

### E. Coli K12 reference (SKIP this step for Evomics2024 as we have already downloaded the files for you!)

We'll get some test data to play with. First, [the E. Coli K12 reference](http://www.ncbi.nlm.nih.gov/nuccore/556503834), from NCBI. It's a bit of a pain to pull out of the web interface so [you can also download it here](http://hypervolu.me/~erik/genomes/E.coli_K12_MG1655.fa).

```bash
# the start of the genome, which is circular but must be represented linearly in FASTA
curl -s http://hypervolu.me/%7Eerik/genomes/E.coli_K12_MG1655.fa | head
# >NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
# AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
# TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA
# ...
# now download and save the genome
curl -s http://hypervolu.me/%7Eerik/genomes/E.coli_K12_MG1655.fa > E.coli_K12_MG1655.fa
```

### E. Coli K12 Illumina 2x300bp MiSeq sequencing results (SKIP this step for Evomics2024 as we have already downloaded the files for you!)

For testing alignment, let's get some data from a [recently-submitted sequencing run on a K12 strain from the University of Exeter](http://www.ncbi.nlm.nih.gov/sra/?term=SRR1770413). We can use the [sratoolkit](https://github.com/ncbi/sratoolkit) to directly pull the sequence data (in paired FASTQ format) from the archive:

```bash
fastq-dump --split-files SRR1770413
```

`fastq-dump` is in the SRA toolkit. It allows directly downloading data from a particular sequencing run ID. SRA stores data in a particular compressed format (SRA!) that isn't directly compatible with any downstream tools, so it's necessary to put things into [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) for further processing. The `--split-files` part of the command ensures we get two files, one for the first and second mate in each pair. We'll use them in this format when aligning.

```bash
# alternatively, you may want to first download, and then dump
# but this seems to fail sometimes for me
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR177/SRR1770413/SRR1770413.sra
sra-dump --split-files SRR1770413.sra
```

These appear to be paired 300bp reads from a modern MiSeq.

### E. Coli O104:H4 HiSeq 2000 2x100bp (SKIP this step for Evomics2024 as we have already downloaded the files for you!)

As a point of comparison, let's also pick up a [sequencing data set from a different E. Coli strain](http://www.ncbi.nlm.nih.gov/sra/SRX095630[accn]). This one is [famous for its role in foodborne illness](https://en.wikipedia.org/wiki/Escherichia_coli_O104%3AH4#Infection) and is of medical interest.

```bash
fastq-dump --split-files SRR341549
```

#### Downsampling for faster processing

To go faster in our tutorial, we'll downsample the FASTQs using `seqtk`:

```
seqkit sample -p 0.2 SRR1770413_1.fastq >SRR1770413_1.subsampled.fastq
seqkit sample -p 0.2 SRR1770413_2.fastq >SRR1770413_2.subsampled.fastq
seqkit sample -p 0.2 SRR341549_1.fastq >SRR341549_1.subsampled.fastq
seqkit sample -p 0.2 SRR341549_2.fastq >SRR341549_2.subsampled.fastq
```

This works because the same seed for the PRNG is used in both runs (`-s 11` is the default).
So we sample the same read pairs in both files even though the commands are independent.
We get files that have pretty low depth, but for our work this won't matter too much and it might make things more interesting.

### Setting up our reference indexes (START here for Evomics 2024!)

#### FASTA file index

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

#### BWA's FM-index

BWA uses the [FM-index](https://en.wikipedia.org/wiki/FM-index), which a compressed full-text substring index based around the [Burrows-Wheeler transform](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform).
To use this index, we first need to build it:

```bash
bwa index E.coli_K12_MG1655.fa
```

You should see `bwa` generate some information about the build process:

```text
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

### Aligning our data against the E. Coli K12 reference

Here's an outline of the steps we'll follow to align our K12 strain against the K12 reference:

1. use bwa to generate SAM records for each read
2. convert the output to BAM
3. sort the output
4. mark PCR duplicates that result from exact duplication of a template during amplification

We could do the steps one-by-one, generating an intermediate file for each step.
However, this isn't really necessary unless we want to debug the process, and it will make a lot of excess files which will do nothing but confuse us when we come to work with the data later.
Thankfully, it's easy to use [unix pipes](https://en.wikiepdia.org/wiki/Pipeline_%28Unix%29) to stream most of these tools together (see this [nice thread about piping bwa and samtools together on biostar](https://www.biostars.org/p/43677/) for a discussion of the benefits and possible drawbacks of this).

You can now run the alignment using a piped approach. _Replace `$threads` with the number of CPUs you would like to use for alignment._ Not all steps in `bwa` run in parallel, but the alignment, which is the most time-consuming step, does. You'll need to set this given the available resources you have.

```bash
bwa mem -t 8 -R '@RG\tID:K12\tSM:K12' \
    E.coli_K12_MG1655.fa SRR1770413_1.subsampled.fastq SRR1770413_2.subsampled.fastq \
    | samtools view -b - >SRR1770413.raw.bam
sambamba sort SRR1770413.raw.bam
sambamba markdup SRR1770413.raw.sorted.bam SRR1770413.bam
```

Breaking it down by line:

- *alignment with bwa*: `bwa mem -t $threads -R '@RG\tID:K12\tSM:K12'` --- this says "align using so many threads" and also "give the reads the read group K12 and the sample name K12"
- *reference and FASTQs* `E.coli_K12_MG1655.fa SRR1770413_1.fastq SRR1770413_2.fastq` --- this just specifies the base reference file name (`bwa` finds the indexes using this) and the input alignment files. The first file should contain the first mate, the second file the second mate.
- *conversion to BAM*: `samtools view -b -` --- this reads SAM from stdin (the `-` specifier in place of the file name indicates this) and converts to BAM.
- *sorting the BAM file*: `sambamba sort SRR1770413.raw.bam` --- sort the BAM file, writing it to `.sorted.bam`.
- *marking PCR duplicates*: `sambamba markdup SRR1770413.raw.sorted.bam SRR1770413.bam` --- this marks reads which appear to be redundant PCR duplicates based on their read mapping position. It [uses the same criteria for marking duplicates as picard](http://lomereiter.github.io/sambamba/docs/sambamba-markdup.html).

Now, run the same alignment process for the O104:H4 strain's data. Make sure to specify a different sample name via the `-R '@RG...` flag incantation to specify the identity of the data in the BAM file header and in the alignment records themselves:

```bash
bwa mem -t 8 -R '@RG\tID:O104_H4\tSM:O104_H4' \
    E.coli_K12_MG1655.fa SRR341549_1.subsampled.fastq  SRR341549_2.subsampled.fastq \
    | samtools view -b - >SRR341549.raw.bam
sambamba sort SRR341549.raw.bam
sambamba markdup SRR341549.raw.sorted.bam SRR341549.bam
```

As a standard post-processing step, it's helpful to add a BAM index to the files. This let's us jump around in them quickly using BAM compatible tools that can read the index. `sambamba` does this for us by default, but if it hadn't or we had used a different process to generate the BAM files, we could use samtools to achieve exactly the same index.

```bash
samtools index SRR1770413.bam
samtools index SRR341549.bam
```

### Using minimap2

It's easy to use `minimap2` instead of `bwa mem`. This may help in some contexts, as it can be several fold faster with minimal reduction in alignment quality. In the case of these short reads, we'd use it as follows. The only major change from bwa mem is that we'll tell it we're working with short read data using `-ax sr`:

```bash
minimap2 -ax sr -t 8 -R '@RG\tID:O104_H4\tSM:O104_H4' \
    E.coli_K12_MG1655.fa SRR341549_1.subsampled.fastq  SRR341549_2.subsampled.fastq \
    | samtools view -b - >SRR341549.raw.minimap2.bam
sambamba sort SRR341549.raw.minimap2.bam
sambamba markdup SRR341549.raw.minimap2.sorted.bam SRR1770413.minimap2.bam
```

```bash
minimap2 -ax sr -t 8 -R '@RG\tID:K12\tSM:K12' \
    E.coli_K12_MG1655.fa SRR1770413_1.subsampled.fastq  SRR1770413_2.subsampled.fastq \
    | samtools view -b - >SRR1770413.raw.minimap2.bam
sambamba sort SRR1770413.raw.minimap2.bam
sambamba markdup SRR1770413.raw.minimap2.sorted.bam SRR1770413.minimap2.bam
```

However, it is widely acknowledged that `minimap2` is not optimal for short read alignment, so use this with caution and in contexts where you can accept lower accuracy in your alignments. Really, you should be applying `minimap2` to long read alignments.

## Part 2: Calling variants

Now that we have our alignments sorted, we can quickly determine variation against the reference by scanning through them using a variant caller.
There are many options, including [samtools mpileup](http://samtools.sourceforge.net/samtools.shtml), [platypus](http://www.well.ox.ac.uk/platypus), and the [GATK](https://www.broadinstitute.org/gatk/).

For this tutorial, we'll keep things simple and use [freebayes](https://github.com/ekg/freebayes). It has a number of advantages in this context (bacterial genomes), such as long-term support for haploid (and polyploid) genomes. However, the best reason to use it is that it's very easy to set up and run, and it produces a very well-annotated VCF output that is suitable for immediate downstream filtering.

### Variant calls with `freebayes`

It's quite easy to use `freebayes` provided you have your BAM file completed. We use `--ploidy 1` to indicate that the sample should be genotyped as haploid.

```bash
freebayes -f E.coli_K12_MG1655.fa --ploidy 1 SRR1770413.bam >SRR1770413.vcf
```

We can run on the other sample too:

```bash
freebayes -f E.coli_K12_MG1655.fa --ploidy 1 SRR341549.bam >SRR341549.vcf
```

### Joint calling

We can put the samples together if we want to find differences between them. Calling them jointly can help if we have a population of samples to use to help remove calls from paralogous regions. The Bayesian model in freebayes combines the data likelihoods from sequencing data with an estimate of the probability of observing a given set of genotypes under assumptions of neutral evolution and a [panmictic](https://en.wikipedia.org/wiki/Panmixia) population. For instance, [it would be very unusual to find a locus at which all the samples are heterozygous](https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle). It also helps improve statistics about observational biases (like strand bias, read placement bias, and allele balance in heterozygotes) by bringing more data into the algorithm.

However, in this context, we only have two samples and the best reason to call them jointly is to make sure we have a genotype for each one at every locus where a non-reference allele passes the caller's thresholds in either sample.

We would run a joint call by dropping in both BAMs on the command line to freebayes:

```bash
freebayes -f E.coli_K12_MG1655.fa --ploidy 1 SRR1770413.bam SRR341549.bam >e_colis.vcf
```

We can do joint calling as long as we've added the read group (@RG) flags when we aligned (or did so after with [bamaddrg](https://github.com/ekg/bamaddrg). This lest the variant caller know which read fits with which sample, and allows it to model the variation in both samples simultaneously.

### `bgzip` and `tabix`

We can speed up random access to VCF files by compressing them with `bgzip`, in the [htslib](https://github.com/samtools/htslib) package.
`bgzip` is a "block-based GZIP", which compresses files in chunks of lines. This chunking let's us quickly seek to a particular part of the file, and support indexes to do so. The default one to use is tabix. It generates indexes of the file with the default name `.tbi`.

```bash
bgzip SRR1770413.vcf # makes SRR1770413.vcf.gz
tabix -p vcf SRR1770413.vcf.gz
```

Now you can pick up a single part of the file. For instance, we could count the variants in a particular region:

```bash
tabix SRR1770413.vcf.gz NC_000913.3:1000000-1500000 | wc -l
```

If we want to pipe the output into a tool that reads VCF, we'll need to add the `-h` flag, to output the header as well.

```bash
# tabix -h SRR1770413.vcf.gz NC_000913.3:1000000-1500000 | vcffilter ...
# example vcf filter
tabix -h SRR1770413.vcf.gz NC_000913.3:1000000-1500000 | vcffilter -f 'DP > 20' | wc -l
```

The `bgzip` format is very similar to that used in BAM, and the indexing scheme is also similar (blocks of compressed data which we build a chromosome position index on top of).

### Take a peek with `vt`

[vt](https://github.com/atks/vt) is a toolkit for variant annotation and manipulation. In addition to other methods, it provides a nice method, `vt peek`, to determine basic statistics about the variants in a VCF file.

We can get a summary like so:

```bash
vt peek SRR1770413.vcf.gz
```

### Filtering using the transition/transversion ratio (ts/tv) as a rough guide

`vt` produces a nice summary with the transition/transversion ratio. Transitions are mutations that switch between DNA bases that have the same base structure (either a [purine](https://en.wikipedia.org/wiki/Purine) or [pyrimidine](https://en.wikipedia.org/wiki/Pyrimidine) ring).

In most biological systems, [transitions (A<->G, C<->T) are far more likely than transversions](https://upload.wikimedia.org/wikipedia/commons/3/35/Transitions-transversions-v3.png), so we expect the ts/tv ratio to be pretty far from 0.5, which is what it would be if all mutations between DNA bases were random. In practice, we tend to see something that's at least 1 in most organisms, and ~2 in some, such as human. In some biological contexts, such as in mitochondria, we see an even higher ratio, perhaps as much as 20.

As we don't have validation information for our sample, we can use this as a simple guide for our first filtering attempts. An easy way is to try different filterings using `vcffilter` and check the ratio of the resulting set with `vt peek`:

```bash
# a basic filter to remove low-quality sites
vcffilter -f 'QUAL > 10' SRR1770413.vcf.gz | vt peek -

# scaling quality by depth is like requiring that the additional log-unit contribution
# of each read is at least N
vcffilter -f 'QUAL / AO > 10' SRR1770413.vcf.gz | vt peek -
```

Note that the second filtering removes a large region near the beginning of the reference where there appears to be some paralogy, which could be caused by a duplication of this region in the sequenced sample relative to the reference genome. The read counts for reference and alternate are each around half of the total depth, which is unusual for a sequenced haploid clone and may indicate some structural differences between the sample and the original reference.

## Part 3: When you know the truth

For serious applications, it's not sufficient to simply filter on the basis of bulk metrics like the ts/tv ratio. Some external validation information should be used to guide the development of pipelines for processing genomic data. In our case, we're just using free data from the web, and unless we find some validation data associated with the strains that were sequenced, we can only filter on intuition, bulk metrics like ts/tv, and with an eye for the particular question we're interested in. What we want is to know the truth for a particular context, so as to understand if our filtering criteria make sense.

### The NIST Genome in a Bottle truth set for NA12878

Luckily, a group at the [National Institute of Standards and Technology](https://en.wikipedia.org/wiki/National_Institute_of_Standards_and_Technology) (NIST) has developed a kind of truth set based on the [HapMap CEU cell line NA12878](https://catalog.coriell.org/0/Sections/Search/Sample_Detail.aspx?Ref=GM12878). It's called the [Genome in a Bottle](https://sites.stanford.edu/abms/giab). In addition to characterizing the genome of this cell line using extensive sequencing and manual curation of inconsistencies found between sequencing protocols, the group actually distributes reference material from the cell line for use in validating sequencing pipelines.

To download the truth set, head over to the [Genome in a Bottle ftp site](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/) and pick up the latest release. As of writing this, we're at [GiAB v3.3.2](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/). Download the highly confident calls and the callable region targets:

```bash
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz{,.tbi}
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed
```

### The human reference (SKIP this step for Evomics2024 as we have already downloaded the files for you!)

For variant calling work, we can use the [1000 Genomes Project's version of the GRCh37 reference](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz). We could also use the [version of the reference that doesn't include dummy sequences](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz), as we're just doing variant calling.

```bash
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
samtools faidx hs37d5.fa
```

### Calling variants in [20p12.1](http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr20%3A12100001-17900000&hgsid=220600397_Vs2XvVv0rRPE9lPwepHAL4Iq3ndi)

To keep things quick enough for the tutorial, let's grab a little chunk of an NA12878 dataset. Let's use [20p12.1](http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr20%3A12100001-17900000&hgsid=220600397_Vs2XvVv0rRPE9lPwepHAL4Iq3ndi). We'll use a downsampled alignment made from Illumina HiSeq sequencing of NA12878 (HG001) that was used as an input to the [NIST Genome in a Bottle](https://github.com/genome-in-a-bottle) truth set for this sample. (Other relevant data can be found in the [GiAB alignment indexes](https://github.com/genome-in-a-bottle/giab_data_indexes).)

We don't need to download the entire BAM file to do this. `samtools` can download the BAM index (`.bai`) provided it hosted alongside the file on the HTTP/FTP server and then use this to jump to a particular target in the remote file.

```bash
samtools view -b ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam 20:12100000-17900000 >NA12878.20p12.1.30x.bam
samtools index NA12878.20p12.1.30x.bam
```

We can call variants as before. Note that we drop the `--ploidy 1` flag. `freebayes` assumes its input is diploid by default. We can use bgzip in-line here to save the extra command for compression.

```bash
freebayes -f hs37d5.fa NA12878.20p12.1.30x.bam | bgzip >NA12878.20p12.1.30x.vcf.gz
tabix -p vcf NA12878.20p12.1.30x.vcf.gz
```

### Comparing our results to the GiAB truth set

We'll need to download the [GiAB truth set](ftp://ftp-trace.ncbi.nih.gov/giab/ftp/release/NA12878_HG001/NISTv2.19/). Its core consists of a VCF file defining "true" variants and a BED file defining callable regions.

In order to compare, we need to exclude things in our output that are outside the callable region, and then intersect with the truth set. That which we don't see in the truth set, and is also in the callable region should be considered a false positive.

First, we'll prepare a reduced representation of this dataset to match 20p12.1:

```bash
# subset the callable regions to chr20 (makes intersection much faster)
cat HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed | grep ^20 >giab_callable.chr20.bed

# subset the high confidence calls to 20p12.1 and rename the sample to match the BAM
tabix -h HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz 20:12100000-17900000 \
    | sed s/HG001/NA12878/ | bgzip >NIST_NA12878_20p12.1.vcf.gz
tabix -p vcf NIST_NA12878_20p12.1.vcf.gz
```

Now, we can compare our results to the calls to get a list of potentially failed sites.

```bash
vcfintersect -r hs37d5.fa -v -i NIST_NA12878_20p12.1.vcf.gz NA12878.20p12.1.30x.vcf.gz \
    | vcfintersect -b giab_callable.chr20.bed \
    | bgzip >NA12878.20p12.1.30x.giab_failed.vcf.gz
tabix -p vcf NA12878.20p12.1.30x.giab_failed.vcf.gz
```

We can now examine these using `vt peek` and `vcfstats`, or manually by inspecting them either serially:

```bash
zcat NA12878.20p12.1.30x.giab_failed.vcf.gz | less -S
```

... or by looking at loci which fail in `samtools tview`.


### Variant normalization

Many of the failed variants are unusual in their normalization. For instance:

```text
20      17432903        .       TATC    TGATA   490.237
```

To ensure that comparisons work correctly, we should "normalize" the variants so that they are represented solely as short indels and SNPs.

There are two main problems:

1. freebayes represents short haplotypes in the VCF output
2. indels may not be completely left-aligned, there could be additional bases on the call that should be removed so that it can be represented in most-normalized form

Finally, the variants in the GiAB set have been normalized using a similar process, and doing so will ensure there are not any discrepancies when we compare.

```bash
vcfallelicprimitives -kg NA12878.20p12.1.30x.vcf.gz \
    | vt normalize -r hs37d5.fa - \
    | bgzip >NA12878.20p12.1.30x.norm.vcf.gz
tabix -p vcf NA12878.20p12.1.30x.norm.vcf.gz
```

Here, `vcfallelicprimitives -kg` decomposes any haplotype calls from `freebayes`, keeping the genotype and site level annotation. (This isn't done by default because in some contexts doing so is inappropriate.) Then `vt normalize` ensures the variants are left-aligned. This isn't important for the comparison, as `vcfintersect` is haplotype-based, so it isn't affected by small differences in the positioning or description of single alleles, but it is good practice.

We can now compare the results again:

```bash
vcfintersect -r hs37d5.fa -v -i NIST_NA12878_20p12.1.vcf.gz NA12878.20p12.1.30x.norm.vcf.gz \
    | vcfintersect -b giab_callable.chr20.bed \
    | bgzip >NA12878.20p12.1.30x.norm.giab_failed.vcf.gz
tabix -p vcf NA12878.20p12.1.30x.norm.giab_failed.vcf.gz
```

Here we observe why normalization is important when comparing VCF files. Fortunately, the best package available for comparing variant calls to truth sets, [rtg-tools](https://www.realtimegenomics.com/products/rtg-tools), addresses exactly this concern, and also breaks comparisons into three parts matching the three types of information provided by the VCF file--- positional, allele, and genotype. We'll get into that later in this section.

### Hard filtering strategies

The failed list provides a means to examine ways to reduce our false positive rate using post-call filtering. We can look at the failed list to get some idea of what might be going on with the failures.

For example, we can test how many of the failed SNPs are removed by applying a simple quality filter and checking the output file's statistics.

```bash
vcffilter -f "QUAL > 10" NA12878.20p12.1.30x.norm.giab_failed.vcf.gz \
    | vt peek -
```

We might also want to measure our sensitivity from different strategies. To do this, just invert the call to `vcfintersect` by removing the `-v` flag (which tells it to invert):

```bash
vcfintersect -r hs37d5.fa -i NIST_NA12878_20p12.1.vcf.gz NA12878.20p12.1.30x.norm.vcf.gz \
    | vcfintersect -b giab_callable.chr20.bed \
    | bgzip >NA12878.20p12.1.30x.norm.giab_passed.vcf.gz
tabix -p vcf NA12878.20p12.1.30x.norm.giab_passed.vcf.gz
```

Now we can test how many variants remain after using the same filters on both:

```bash
vcffilter -f "QUAL / AO > 10 & SAF > 0 & SAR > 0" NA12878.20p12.1.30x.norm.giab_passed.vcf.gz | vt peek -
vcffilter -f "QUAL / AO > 10 & SAF > 0 & SAR > 0" NA12878.20p12.1.30x.norm.giab_failed.vcf.gz | vt peek -
```

### Using RTG-eval for comparison to the truth

The best variant calling comparison and evaluation framework in current use was developed by Real Time Genomics: [rtg-tools](https://www.realtimegenomics.com/products/rtg-tools). This was subsequently used for the basis of comparison in the PrecisionFDA challenges in 2016.

We can easily apply `rtg vcfeval` to our results, but we will need to prepare the reference in RTG's "SDF" format first.

```bash
rtg format -o hs37d5.sdf hs37d5.fa  # done already at evomics2024 --- takes a while so don't re-do!
```

Now we can proceed and test the performance of our previous freebayes run against the truth set.

```bash
rtg vcfeval -f QUAL -o eval1 -t hs37d5.sdf -e giab_callable.chr20.bed \
    -b NIST_NA12878_20p12.1.vcf.gz -c NA12878.20p12.1.30x.vcf.gz
```

The output of `rtg vcfeval` is a set of reports and files tallying true and false positives.

We also get a nice report in the shell that describes the optimal QUAL cutoff (here equivalent to `vcffilter -f 'QUAL > 3.587'`) and the result if we have no cutoff (`None`).

```
Selected score threshold using: maximized F-measure
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
    3.587               8861           8686          7         18     0.9992       0.9980     0.9986
     None               8868           8693        110         11     0.9875       0.9988     0.9931
```

In this case, we can get a quick overview by looking in the files and directories prefixed by `eval1`. It is also quick to clean up with `rm -rf eval1.*`. _Make sure you clean up before re-running on a new file, or use a different prefix!_

### Bonus: new freebayes configurations, ROC plots, and false positive investigations

Take the output of `rtg vcfeval` in `eval1/` and complete some of the following tasks:

1. Plot the ROC curve based on the QUAL field that's given in `eval1/snp_roc.tsv.gz`.
2. Look at "false positives" that are over the optimal QUAL threshold. These are in `eval1/fp.vcf.gz`. Take a few and use `samtools tview` or `igv` to look at the alignments around the putative error. What's going on? If anyone completes this we can discuss as a group.
3. Look at "false negatives" in the `fn.vcf.gz` file and see what's happening in the alignments around each locus. We can also discuss some of these together.
4. Re-run `freebayes` with different configurations and use `rtg vcfeval` to see if you can get a better F-measure than the default settings!
