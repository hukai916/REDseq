\name{REDseq-package}
\alias{REDseq-package}
\alias{REDseq}
\docType{package}
\title{
REDseq
}
\description{
REDSeq is a Bioconductor package for building genomic map of restriction enzyme sites REmap, assigning sequencing tags to RE sites using five different strategies, visualizing genome-wide distribution of differentially cut regions with the REmap as reference and the distance distribution of sequence tags to corresponding RE sites, generating count table for identifying statistically significant RE sites using edgeR or DEseq.
}
\details{
\tabular{ll}{
Package: \tab REDseq\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-05-10\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Lihua Julie Zhu

Maintainer:  
Lihua Julie Zhu <julie.zhu@umassmed.edu>
}
\references{
1. Roberts, R.J., Restriction endonucleases. CRC Crit Rev Biochem, 1976. 4(2): p. 123-64.
7. Pages, H., BSgenome package. http://bioconductor.org/packages/2.8/bioc/\cr
vignettes/BSgenome/inst/doc/GenomeSearching.pdf \cr
8. Zhu, L.J., et al., REDseq: A Bioconductor package for Analyzing High Throughput Sequencing Data from Restriction Enzyme Digestion. (In preparation) \cr

\keyword{ package }
\seealso{
buildREmap, assignSeq2REsit, plotCutDistribution, distanceHistSeq2RE, summarizeByRE, summarizeBySeq, compareREseq, binom.test.REDseq
}
\examples{
	if(interactive()){
	library(ChIPpeakAnno)
	REpatternFilePath = system.file("extdata", "examplePattern.fa", package="REDseq")
	library(BSgenome.Celegans.UCSC.ce2)
	buildREmap( REpatternFilePath, BSgenomeName=Celegans, outfile=tempfile())
	library(REDseq)
	data(example.REDseq)
	data(example.map)
	r.unique = assignSeq2REsite(example.REDseq, example.map, cut.offset = 1, 
seq.length = 36, allowed.offset = 5, min.FragmentLength = 60, 
max.FragmentLength = 300, partitionMultipleRE = "unique")
	r.average = assignSeq2REsite(example.REDseq, example.map, cut.offset = 1,
 seq.length = 36, allowed.offset = 5, min.FragmentLength = 60, 
max.FragmentLength = 300, partitionMultipleRE = "average")
	r.random = assignSeq2REsite(example.REDseq, example.map, cut.offset = 1,
 seq.length = 36, allowed.offset = 5, min.FragmentLength = 60, 
max.FragmentLength = 300, partitionMultipleRE = "random")
	r.best = assignSeq2REsite(example.REDseq, example.map, cut.offset = 1,
 seq.length = 36, allowed.offset = 5, min.FragmentLength = 60, 
max.FragmentLength = 300, partitionMultipleRE = "best")
	r.estimate = assignSeq2REsite(example.REDseq, example.map, cut.offset = 1,
 seq.length = 36, allowed.offset = 5, min.FragmentLength = 60, 
max.FragmentLength = 300, partitionMultipleRE = "estimate")
	r.estimate$passed.filter
	r.estimate$notpassed.filter
	data(example.assignedREDseq)
	plotCutDistribution(example.assignedREDseq,example.map, 
chr="2", xlim =c(3012000, 3020000))
	distanceHistSeq2RE(example.assignedREDseq,ylim=c(0,20))
	summarizeByRE(example.assignedREDseq,by="Weight",sampleName="example")
	REsummary  =summarizeByRE(example.assignedREDseq,by="Weight")
	binom.test.REDseq(REsummary)
}
}