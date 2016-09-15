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
CBamScaffold = function(oBamFile, cSeqname, iBinSize=1000, flag = scanBamFlag(), what = scanBamWhat(), fViewSummary=viewMaxs,...){
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
  ivReadWidth = width(oGA)
  lWhat = vector('list', length=length(what))
  names(lWhat) = what
  # remove objects 
  rm(gr, oGA, cov)
  # clear memory
  gc()
  return(new('CBamScaffold', ivCoverage=ivCoverage, ivReadWidth=ivReadWidth, cSeqname=cSeqname, cBamFile=cBamFile, 
         lWhat=lWhat))
}
