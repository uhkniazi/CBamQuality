# CBamQuality
Quality metrics and plots for a BAM file generated after alignment of FASTQ file.

# S4 Class: CBamScaffold
### Slots  
**ivCoverage** = binned vector with coverage over each bin  
**ivReadWidth** = integer vector of the width of each read  
**cSeqname** = character string with sequence/scaffold name  
**cBamFile** = file name and path  
**lWhat** = list, extendible for future use  

### Constructor  
```R
CBamScaffold = function(oBamFile, cSeqname, iBinSize=1000, flag = scanBamFlag(), what = scanBamWhat(), fViewSummary=viewMeans,...)
```  
### ARGS  
**oBamFile** = object of class BamFile  
**cSeqname** = character string with sequence/scaffold name  
**iBinSize** = the number of bins to bin the scaffold in  
**flag, what** = parameters for the ScanBamParam(flag=flag, what=what, which=which) function  
**fViewSummary** = the function to call to get the summary, other options are viewMeans, viewMins, viewSums, viewMaxs  

### Internal function  
```R
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
```
---
### Slot accessor functions  
```R
CBamScaffold.getCoverage = function(obj) 
CBamScaffold.getReadWidthSample = function(obj, size=1000) 
CBamScaffold.getReadWidth = function(obj) 
CBamScaffold.getSeqname = function(obj) 
CBamScaffold.getBamPath = function(obj) 
```
  
  
### Generic and plotting functions  
---
### getCoverageGammaParam  
```R
setGeneric('getCoverageGammaParam', function(obj, ...)standardGeneric('getCoverageGammaParam'))
```  

**RETS**  returns the shape parameter for the theoretical gamma distribution fitted to the distribution of the binned coverage.  

---
### plot.coverage.distribution  
```R
setGeneric('plot.coverage.distribution', function(obj, title='', ...)standardGeneric('plot.coverage.distribution'))
```  
**DESC** uses the vector of binned coverage (see *iBinSize=1000, fViewSummary=viewMeans*) to plot the histogram of discrete trimmed (at 0.95 quantile) distribution, and plots the gamma density on top to show the fit. 

![alt text](https://github.com/uhkniazi/CBamQuality/blob/master/SamplePlots/gamma.png "Binned empirical and Gamma distribution fit")  
---
### plot.coverage.vector  
```R
setGeneric('plot.coverage.vector', function(obj, title='', bLog=F, type='l', ...)standardGeneric('plot.coverage.vector'))
```
**DESC** Makes a plot of the binned coverage vector either on normal or log scale.  
One can also plot this as a lowess or any other curve fitting function.  
```R
# get matrix/vector of coverage
m = CBamScaffold.getCoverage(ob)
m2 = lowess(m, f=2/10)$y
plot(m2, type='l', main='lowess fit to average binned coverage')
```  
![alt text](https://github.com/uhkniazi/CBamQuality/blob/master/SamplePlots/binned_coverage.png "Binned Coverage") 
![alt text](https://github.com/uhkniazi/CBamQuality/blob/master/SamplePlots/binned_coverage_lowess.png "Binned Lowess Coverage") 

