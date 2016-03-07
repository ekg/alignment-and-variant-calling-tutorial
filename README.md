# NGS alignment and variant calling

This tutorial steps through some basic tasks in alignment and variant calling using a handful of Illumina sequencing data sets. For theoretical background, please refer to the included [presentation on alignment and variant calling](https://github.com/ekg/alignment-and-variant-calling-tutorial/raw/master/presentations/Alignment%2C%20Variant%20Calling%2C%20and%20Filtering%20(WGC%20NGS%20Bioinformatics).pdf).
    
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
9. [mutatrix](https://github.com/ekg/mutatrix)
8. [sra-tools](https://github.com/ncbi/sra-tools/wiki/HowTo:-Binary-Installation)

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

In some cases you can download precompiled binaries. For instance, you can head
to the [sambamba releases page](https://github.com/lomereiter/sambamba/releases) to find binaries that
should work on any modern Linux or OSX distribution. The same applies to the [sra-toolkit, which is probably easier to install from available binaries](https://github.com/ncbi/sra-tools/wiki/HowTo:-Binary-Installation).

Otherwise, let's assume you're in an environment where you've already got them available.

## Part 1: Aligning E. Coli data with `bwa mem`

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

For testing alignment, let's get some data from a [recently-submitted sequencing run on a K12 strain from the University of Exeter](http://www.ncbi.nlm.nih.gov/sra/?term=SRR1770413). We can use the [sratoolkit](https://github.com/ncbi/sratoolkit) to directly pull the sequence data (in paired FASTQ format) from the archive:

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

### Setting up our reference indexes

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

We could the steps one-by-one, generating an intermediate file for each step.
However, this isn't really necessary unless we want to debug the process, and it will make a lot of excess files which will do nothing but confuse us when we come to work with the data later.
Thankfully, it's easy to use [unix pipes](https://en.wikiepdia.org/wiki/Pipeline_%28Unix%29) to stream most of these tools together (see this [nice thread about piping bwa and samtools together on biostar](https://www.biostars.org/p/43677/) for a discussion of the benefits and possible drawbacks of this).

You can now run the alignment using a piped approach. _Replace `$threads` with the number of CPUs you would like to use for alignment._ Not all steps in `bwa` run in parallel, but the alignment, which is the most time-consuming step, does. You'll need to set this given the available resources you have.

```bash
bwa mem -t $threads -R '@RG\tID:K12\tSM:K12' \
    E.coli_K12_MG1655.fa SRR1770413_1.fastq.gz SRR1770413_2.fastq.gz \
    | samtools view -Shu - \
    | sambamba sort /dev/stdin -o SRR1770413.raw.bam
sambamba markdup SRR1770413.raw.bam SRR1770413.bam
```

Breaking it down by line:

- *alignment with bwa*: `bwa mem -t $threads -R '@RG\tID:K12\tSM:K12'` --- this says "align using so many threads" and also "give the reads the read group K12 and the sample name K12"
- *reference and FASTQs* `E.coli_K12_MG1655.fa SRR1770413_1.fastq.gz SRR1770413_2.fastq.gz` --- this just specifies the base reference file name (`bwa` finds the indexes using this) and the input alignment files. The first file should contain the first mate, the second file the second mate.
- *conversion to BAM*: `samtools view -Shu -` --- this reads SAM from stdin (`-S` and the `-` specifier in place of the file name indicate this) and converts to uncompressed BAM (there isn't need to compress, as it's just going to be parsed by the next program in the pipeline.
- *sorting the BAM file*: `sambamba sort /dev/stdin -o /dev/stdout` --- sort the BAM file, reading from stdin and writing to stdout. We then use shell redirection to write it to a file, but you could put tools that can work on streams directly after this step.
- *marking PCR duplicates*: `sambamba markdup SRR1770413.raw.bam SRR1770413.bam` --- this marks reads which appear to be redundant PCR duplicates based on their read mapping position. It [uses the same criteria for marking duplicates as picard](http://lomereiter.github.io/sambamba/docs/sambamba-markdup.html).

Now, run the same alignment process for the O104:H4 strain's data. Make sure to specify a different sample name via the `-R '@RG...` flag incantation to specify the identity of the data in the BAM file header and in the alignment records themselves:

```bash
bwa mem -t $threads -R '@RG\tID:O104_H4\tSM:O104_H4' \
    E.coli_K12_MG1655.fa SRR341549_1.fastq.gz  SRR341549_2.fastq.gz \
    | samtools view -Shu - \
    | sambamba sort /dev/stdin -o /dev/stdout >SRR341549.raw.bam
sambamba markdup SRR341549.raw.bam SRR341549.bam
```

As a standard post-processing step, it's helpful to add a BAM index to the files. This let's us jump around in them quickly using BAM compatible tools that can read the index. `sambamba` does this for us by default, but if it hadn't we could use samtools to achieve exactly the same index.

```bash
samtools index SRR1770413.bam
samtools index SRR341549.bam
```

## Part 2: Calling variants

Now that we have our alignments sorted, we can quickly determine variation against the reference by scanning through them using a variant caller.
There are many options, including [samtools mpileup](http://samtools.sourceforge.net/samtools.shtml), [platypus](http://www.well.ox.ac.uk/platypus), and the [GATK](https://www.broadinstitute.org/gatk/).

For this tutorial, we'll keep things simple and use [freebayes](https://github.com/ekg/freebayes). It has a number of advantages in this context (bacterial genomes), such as long-term support for haploid (and polyploid) genomes. However, the best reason to use it is that it's very easy to set up and run, and it produces a very well-annotated VCF output that is suitable for immediate downstream filtering.

### Variant calls with `freebayes`

It's quite easy to use `freebayes` provided you have your BAM file completed. We use `--ploidy 1` to indicate that the sample should be genotyped as haploid.

```bash
freebayes -f E.coli_K12_MG1655.fa --ploidy 1 SRR1770413.bam >SRR1770413.vcf
```

### Joint calling

We can put the samples together if we want to find differences between them. Calling them jointly can help if we have a population of samples to use to help remove calls from paralagous regions. The Bayesian model in freebayes combines the data likelihoods from sequencing data with an estimate of the probability of observing a given set of genotypes under assumptions of neutral evolution and a [panmictic](https://en.wikipedia.org/wiki/Panmixia) population. For instance, [it would be very unusual to find a locus at which all the samples are heterozygous](https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle). It also helps improve statistics about observational biases (like strand bias, read placement bias, and allele balance in heterozygotes) by bringing more data into the algorithm.

However, in this context, we only have two samples and the best reason to call them jointly is to make sure we have a genotype for each one at every locus where a non-reference allele passes the caller's thresholds in either sample.

We would run a joint call by dropping in both BAMs on the command line to freebayes:

```bash
freebayes -f E.coli_K12_MG1655.fa --ploidy 1 SRR1770413.bam SRR341549.bam >e_colis.vcf
```

As long as we've added the read group (@RG) flags when we aligned (or did so after with [bamaddrg](https://github.com/ekg/bamaddrg), that's all we need to do to run the joint calling. (NB: due to the amount of data in SRR341549, this should take about 20 minutes.)

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
tabix -h SRR1770413.vcf.gz NC_000913.3:1000000-1500000 | vcffilter ...
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

![DNA mutations](https://upload.wikimedia.org/wikipedia/commons/3/35/Transitions-transversions-v3.png)

In most biological systems, transitions (A<->G, C<->T) are far more likely than transversions, so we expect the ts/tv ratio to be pretty far from 0.5, which is what it would be if all mutations between DNA bases were random. In practice, we tend to see something that's at least 1 in most organisms, and ~2 in some, such as human. In some biological contexts, such as in mitochondria, we see an even higher ratio, perhaps as much as 20.

As we don't have validation information for our sample, we can use this as a simple guide for our first filtering attempts. An easy way is to try different filterings using `vcffilter` and check the ratio of the resulting set with `vt peek`:

```bash
# a basic filter to remove low-quality sites
vcffilter -f 'QUAL > 10' SRR1770413.vcf.gz | vt peek -

# scaling quality by depth is like requiring that the additional log-unit contribution
# of each read is at least N
vcffilter -f 'QUAL / AO > 10' SRR1770413.vcf.gz | vt peek -
```

Note that the second filtering removes a large region near the beginning of the reference where there appears to be some paralogy. The read counts for reference and alternate aare each around half of the total depth, which is unusual for a sequenced clone and may indicate some structural differences between the sample and the original reference.

## Part 3: When you know the truth

For serious applications, it's not sufficient to simply filter on the basis of bulk metrics like the ts/tv ratio. Some external validation information should be used to guide the development of pipelines for processing genomic data. In our case, we're just using free data from the web, and unless we find some validation data associated with the strains that were sequenced, we can only filter on intuition, bulk metrics like ts/tv, and with an eye for the particular question we're interested in. What we want is to know the trut for a particular context, so as to understand if our filtering criteria make sense.

### The NIST Genome in a Bottle truth set for NA12878

Luckily, a group at the [National Institute of Standards and Technology](https://en.wikipedia.org/wiki/National_Institute_of_Standards_and_Technology) (NIST) has developed a kind of truth set based on the [HapMap CEU cell line NA12878](https://catalog.coriell.org/0/Sections/Search/Sample_Detail.aspx?Ref=GM12878). It's called the [Genome in a Bottle](https://sites.stanford.edu/abms/giab). In addition to characterizing the genome of this cell line using extensive sequencing and manual curation of inconsistencies found between sequencing protocols, the group actually distributes reference material from the cell line for use in validating sequencing pipelines.

To download the truth set, head over to the [Genome in a Bottle ftp site](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/) and pick up the latest release. As of writing this, we're at [GiAB v2.19](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv2.19/). Download the highly confident calls and the callable region targets:

```bash
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv2.19/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv2.19/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio.bed.gz
```

### The human reference

For variant calling work, we can use the [1000 Genomes Project's version of the GRCh37 reference](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz). We could also use the [version of the reference that doesn't include dummy sequences](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz), as we're just doing variant calling.

```bash
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
samtools faidx hs37d5.fa
```

### Calling variants in [20p12.2](http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr20%3A9200001-12100000)

To keep things quick enough for the tutorial, let's grab a little chunk of an NA12878 dataset. Let's use [20p12.2](http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr20%3A9200001-12100000).

We don't need to download the entire BAM file to do this. `samtools` can download the BAM index (`.bai`) provided it hosted alongside the file on the HTTP/FTP server and then use this to jump to a particular target in the remote file.

```bash
samtools view -b ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NA12878_data_other_projects/alignment/XPrize_Illumina_WG.bam 20:9200001-12100000 >NA12878.20p12.2.XPrize.bam
```

We can call variants as before. Note that we drop the `--ploidy 1` flag. `freebayes` assumes its input is diploid by default. We can use bgzip in-line here to save the extra command for compression.

```bash
freebayes -f hs37d5.fa NA12878.20p12.2.XPrize.bam | bgzip >NA12878.20p12.2.XPrize.vcf.gz
tabix -p vcf NA12878.20p12.2.XPrize.vcf.gz
```

### Comparing our results to the GiAB truth set

We'll need to download the [GiAB truth set](ftp://ftp-trace.ncbi.nih.gov/giab/ftp/release/NA12878_HG001/NISTv2.19/). Its core consists of a VCF file defining "true" variants and a BED file defining callable regions.

In order to compare, we need to exclude things in our output that are outside the callable region, and then intersect with the truth set. That which we don't see in the truth set, and is also in the callable region should be considered a false positive.

First, we'll prepare a reduced representation of this dataset to match 20p12.2:

```bash
# subset the callable regions to chr20 (makes intersection much faster)
zcat union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio.bed.gz | grep ^20 >giab_callable.chr20.bed

# index the high-confidence calls
tabix -p vcf NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz

# and subset them to 20p12.2
tabix -h NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz 20:9200001-12100000 \
    | bgzip >NIST_NA12878_20p12.2.vcf.gz
tabix -p vcf NIST_NA12878_20p12.2.vcf.gz
```

Now, we can compare our results to the calls to get a list of potentially failed sites.

```bash
vcfintersect -r hs37d5.fa -v -i NIST_NA12878_20p12.2.vcf.gz NA12878.20p12.2.XPrize.vcf.gz \
    | vcfintersect -b giab_callable.chr20.bed \
    | bgzip >NA12878.20p12.2.XPrize.giab_failed.vcf.gz
tabix -p vcf NA12878.20p12.2.XPrize.giab_failed.vcf.gz
```

We can now examine these using `vt peek` and `vcfstats`, or manually by inspecting them either serially:

```bash
zcat NA12878.20p12.2.XPrize.giab_failed.vcf.gz | less -S
```

... or by looking at loci which fail in `samtools tview`.

### Variant normalization

Many of the failed variants are unusual in their normalization. For instance:

```text
20   9575773 .  GAGAG  TATAT  1172.52
```

To ensure that comparisons work correctly, we should "normalize" the variants so that they are represented solely as short indels and SNPs.

There are two main problems:

1. freebayes represents short haplotypes in the VCF output
2. indels may not be completely left-aligned, there could be additional bases on the call that should be removed so that it can be represented in most-normalized form

Finally, the variants in the GiAB set have been normalized using a similar process, and doing so will ensure there are not any discrepancies when we compare.

```bash
vcfallelicprimitives -kg NA12878.20p12.2.XPrize.vcf.gz \
    | vt normalize -r hs37d5.fa - \
    | bgzip >NA12878.20p12.2.XPrize.norm.vcf.gz
```

Here, `vcfallelicprimitives -kp` decomposes any haplotype calls from `freebayes`, keeping the genotype and site level annotation. (This isn't done by default because in some contexts doing so is inappropriate.) Then `vt normalize` ensures the variants are left-aligned. This isn't important for the comparison, as `vcfintersect` is haplotype-based, so it isn't affected by small differences in the positioning or descripition of single alleles, but it is good practice.

We can now compare the results again:

```bash
vcfintersect -r hs37d5.fa -v -i NIST_NA12878_20p12.2.vcf.gz NA12878.20p12.2.XPrize.norm.vcf.gz \
    | vcfintersect -b giab_callable.chr20.bed \
    | bgzip >NA12878.20p12.2.XPrize.norm.giab_failed.vcf.gz
tabix -p vcf NA12878.20p12.2.XPrize.norm.giab_failed.vcf.gz
```

### Hard filtering strategies

The failed list provides a means to examine ways to reduce our false positive rate using post-call filtering. We can look at the failed list to get some idea of what might be going on with the failures.

For example, we can test how many of the failed SNPs are removed by applying a simple quality filter and checking the output file's statistics.

```bash
vcffilter -f "QUAL > 10" NA12878.20p12.2.XPrize.norm.giab_failed.vcf.gz \
    | vt peek -
```

We might also want to measure our sensitivity from different strategies. To do this, just invert the call to `vcfintersect` by removing the `-v` flag (which tells it to invert):

```bash
vcfintersect -r hs37d5.fa -i NIST_NA12878_20p12.2.vcf.gz NA12878.20p12.2.XPrize.norm.vcf.gz \
    | vcfintersect -b giab_callable.chr20.bed \
    | bgzip >NA12878.20p12.2.XPrize.norm.giab_passed.vcf.gz
tabix -p vcf NA12878.20p12.2.XPrize.norm.giab_passed.vcf.gz
```

Now we can test how many variants remain after using the same filters on both:

```bash
vcffilter -f "QUAL / AO > 10 & SAF > 0 & SAR > 0" NA12878.20p12.2.XPrize.norm.giab_passed.vcf.gz | wc -l
vcffilter -f "QUAL / AO > 10 & SAF > 0 & SAR > 0" NA12878.20p12.2.XPrize.norm.giab_failed.vcf.gz | wc -l
```

## Part 4: Exploring alternative genomic contexts

`freebayes` is capable of handling a wide array of genomic contexts, for instance, pooled experiments and polyploid genomes.

To evaluate these, few truth sets exist. One option is to use simulation to generate a dataset we can work from.

### Simulation with mutatrix and wgsim

First, we'll want to pick up a reference genome that's pretty small so that we can quickly iterate the tests.

```bash
# makes a 100kb fasta sequence from E. coli K12
samtools faidx E.coli_K12_MG1655.fa NC_000913.3:1000000-1100000 >NC_000913.3:1000000-1100000.fa
# index it
samtools faidx NC_000913.3:1000000-1100000.fa
```

Now we can simulate genomes using mutatrix.

```bash
mutatrix -S sample -p 2 -n 10 NC_000913.3:1000000-1100000.fa | bgzip >answers.vcf.gz
tabix -p vcf answers.vcf.gz
```

[This example may be helpful to see how the whole thing works](https://github.com/ekg/mutatrix/blob/master/freebayes_mosaik_simulation.sh).

### A simulated population

### Polyploids and pooled sequencing contexts

## errata

If you're part of the Biology for Adaptation genomics course, [here is a shared document describing system-specific information about available data sets and binaries](http://goo.gl/JvNIRv).
