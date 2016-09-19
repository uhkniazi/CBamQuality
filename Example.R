# Name: Example.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 19/9/2016
# Desc: example usage of bam quality checks

source('CBamQuality.R')

## from the vigenette of GenomicAlignments
csFile = system.file('extdata', 'ex1.bam', package='Rsamtools')

# load the bam file as GAlignment object
oGABam = readGAlignments(csFile)

# load the bam file using BamFile function or BamFileList
bf = BamFile(csFile)
# get names of sequences
sn = seqinfo(bf)
gr = as(sn, 'GRanges')
seqlengths(bf)

## open a larger bamfile 
csFile = file.choose()
bf = BamFile(csFile)
sn = seqinfo(bf)
gr = as(sn, 'GRanges')

# we loaded a human rna-seq alignment file
# looking at chr1 i.e. one scaffold using CBamScaffold class

ob = CBamScaffold(bf, 'chr1')

# create some sample plots
plot.coverage.distribution(ob, 'chr1 sample')
plot.coverage.vector(ob, 'chr1 not logged')
plot.coverage.vector(ob, 'chr1 logged', bLog = T)

# plot coverage in a different manner
# get matrix/vector of coverage
m = CBamScaffold.getCoverage(ob)
m2 = lowess(m, f=2/10)$y

plot(m2, type='l', main='lowess fit to average binned coverage')

# get the gamma shape parameter for the binned data
getCoverageGammaParam(ob)

# average read width
boxplot(CBamScaffold.getReadWidthSample(ob), main='chr1 read width distribution')

# total number of reads in scaffold in millions
length(CBamScaffold.getReadWidth(ob))/1e+6

## for more usage examples over a whole set of scaffolds see example
## https://github.com/uhkniazi/BRC_Keloid/blob/master/Keloid_main/bam_files_qa.R