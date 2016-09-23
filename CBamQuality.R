# Name: CBamQuality.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 15/9/2016
# Desc: class to create a CBamQuality and CBamScaffold Objects 


library(methods)
if (!require(Rsamtools) || !require(GenomicAlignments)) stop('BamQuality.R: libraries missing Rsamtools, GenomicAlignments')

##### Class CBamScaffold
# Name: Class CBamScaffold
# Desc: Uses the functionality in the Rsamtools and GenomicAlignments libraries to create an object to hold information 
#       about the sequence/scaffold

# ivCoverage = binned vector with coverage over each bin
# ivReadWidth = integer vector of the width of each bin
# cSeqname = character string with sequence name
# cBamFile = file name and path
# lWhat = list, extendible for future use
setClass('CBamScaffold', slots=list(ivCoverage='numeric', ivReadWidth='numeric', cSeqname='character', cBamFile='character', lWhat='list'))


## constructor
## ARGS:
## oBamFile = object of class BamFile
## cSeqname = character string with sequence/scaffold name
## iBinSize = the number of bins to bin the scaffold in
## flag, what = parameters for the ScanBamParam(flag=flag, what=what, which=which) function
## fViewSummary = the function to call to get the summary, other options are viewMeans, viewMins, viewSums, viewMaxs
CBamScaffold = function(oBamFile, cSeqname, iBinSize=1000, flag = scanBamFlag(), what = scanBamWhat(), fViewSummary=viewMeans,...){
  ######### internal functions
  # Function: f_bin_vector
  # DESC: Takes the start and end values of the vector; and the number of
  #       bins, and returns a data.frame with the start and ends of the 
  #       bins
  # ARGS: start coordinate, end coordinate, and number of bins
  # RETS: a data.frame object with starts and ends of each bin
  f_bin_vector = function(start, end, bins){
    s = floor(seq(start, end, length.out=bins+1))
    e = s-1
    e[length(e)] = s[length(s)]
    length(s) = length(s)-1
    e = e[2:length(e)]
    return(data.frame(start=s, end=e))
  }# f_bin_vector
  ############################ internal functions
  
  ## read the name of the bamfile with path
  if (class(oBamFile) != 'BamFile') stop('CBamScaffold: oBamFile is not object of class BamFile')
  cBamFile = path(oBamFile)
  # check if sequence is present in bam file
  gr = as(seqinfo(oBamFile), 'GRanges')
  if(!any(as.character(seqnames(gr)) == cSeqname)) stop('CBamScaffold: cSeqname is not present in Sequence names')
  which = gr[cSeqname]
  param = ScanBamParam(flag=flag, what=what, which=which)
  # read the GAlignments object
  oGA = readGAlignments(oBamFile, param=param)
  # get the coverage
  cov = coverage(oGA)[[cSeqname]]
  # create a binned vector to create views on this coverage
  bins = f_bin_vector(start(which), end(which), bins=iBinSize)
  # create views on these bins
  ivCoverage = fViewSummary(Views(cov, bins$start, bins$end)) 
  ## there is an option to store only a sample of the read width to save on memory 
  ivReadWidth = qwidth(oGA)
  lWhat = vector('list', length=length(what))
  names(lWhat) = what
  # remove objects 
  rm(gr, oGA, cov)
  # clear memory
  gc()
  return(new('CBamScaffold', ivCoverage=ivCoverage, ivReadWidth=ivReadWidth, cSeqname=cSeqname, cBamFile=cBamFile, 
         lWhat=lWhat))
}

### slot accessors
CBamScaffold.getCoverage = function(obj) obj@ivCoverage
CBamScaffold.getReadWidthSample = function(obj, size=1000) sample(obj@ivReadWidth, size = size)
CBamScaffold.getReadWidth = function(obj) obj@ivReadWidth
CBamScaffold.getSeqname = function(obj) obj@cSeqname
CBamScaffold.getBamPath = function(obj) obj@cBamFile

### generic and plotting functions

setGeneric('getCoverageGammaParam', function(obj, ...)standardGeneric('getCoverageGammaParam'))
setMethod('getCoverageGammaParam', signature = 'CBamScaffold', definition = function(obj, ...){
  # get the vector of binned coverage
  t = CBamScaffold.getCoverage(obj)
  # trim the top 95%
  t = t[t < quantile(t, prob=c(0.95))]
  param = getalphabeta.poisson(t)
  names(param) = c('shape', 'rate')
  return(param)
})

## plot the distribution of binned coverage
setGeneric('plot.coverage.distribution', function(obj, title='', ...)standardGeneric('plot.coverage.distribution'))
setMethod('plot.coverage.distribution', signature = 'CBamScaffold', definition = function(obj, title='', ...){
  # get the vector of binned coverage
  t = CBamScaffold.getCoverage(obj)
  # trim the top 95%
  t = t[t < quantile(t, prob=c(0.95))]
  r = range(t)
  s = seq(max(0.5, floor(r[1]))-0.5, ceiling(r[2])+0.5, by=1)
  # calculate the mid points for histogram/discrete distribution
  h = hist(t, breaks=s, plot=F)
  param = getCoverageGammaParam(obj)
  dg = dgamma(h$mids, shape = param['shape'], rate = param['rate'])
  # which distribution can approximate the frequency
  hist(t, prob=T, sub='Distribution of Binned Coverage', breaks=s, main=paste(title, 'Sequence', obj@cSeqname),
       xlab='Coverage', ylab='', ylim=c(0, max(dg, h$density)))
  # parameterized on the means
  lines(h$mids, dg, col='black', type='b')
  points(qgamma(0.95, mean(t), 1), 0, pch=20, col='red')
  legend('topright', legend =c('Gamma'), fill = c('black'))
})

## plot the coverage vector
setGeneric('plot.coverage.vector', function(obj, title='', bLog=F, type='l', ...)standardGeneric('plot.coverage.vector'))
setMethod('plot.coverage.vector', signature = 'CBamScaffold', definition = function(obj, title='', bLog=F, type='l', ...){
  # get the vector of binned coverage
  t = CBamScaffold.getCoverage(obj)
  # threshold the top 95% if bLog is FALSE
  if (!bLog) {t[t > quantile(t, prob=c(0.95))] = quantile(t, prob=c(0.95))
  } else t = log(t+1)
  plot(t, type=type, sub='Binned Coverage', ylab='Coverage', xlab='Bins', main=paste(title, 'Sequence', obj@cSeqname))
})

## utility functions
# calculates the gamma prior parameters for a poisson sampling distribution
# see page 5 in notes here: https://www.evernote.com/shard/s288/res/659dc700-ccb9-457e-aea6-aa54bc3abbb9
# and for an example see page 154, chapter on Hierarchical Modeling Bayesian Computation with R.
## DESC: using the poisson sampling model, the data vector is used to count values of alpha (shape), beta (rate)
## parameters for the gamma prior
getalphabeta.poisson = function(lambda){
  m = mean(lambda)
  v = var(lambda)
  alpha = (m^2)/v
  beta = alpha/m
  return(c(alpha=alpha, beta=beta))
}
