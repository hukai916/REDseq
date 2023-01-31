buildREmap <-
function(REpatternFilePath,format="fasta", BSgenomeName, outfile, chr = c("all"))
{
		if(missing(REpatternFilePath))
		{
			stop("missing required parameter REpatternFilePath!")
		}
		if (!file.exists(REpatternFilePath))
		{
			stop("REpatternFilePath specified as ", REpatternFilePath, " does not exsist!")
		}
		if (format != "fasta" && format != "fastq")
		{
			stop("format needs to be either fasta or fastq!")
		}
		if (missing(BSgenomeName) || class(BSgenomeName) != "BSgenome")
		{
			stop("BSgenomeName is required as BSgenome object!")
		}
		if (file.exists(outfile))
		{
			outfile = paste("REDseq_buildREmap", outfile,sep="_")
			#stop("outfile specified as ", outfile, " already exists! Please rename the outfile!")
		}
		dict = readDNAStringSet(REpatternFilePath, format, use.names=TRUE)
		searchPattern(dict, BSgenomeName=BSgenomeName,outfile=outfile, chr = chr)
		REmap <- toGRanges(outfile, format="BED", header=TRUE, sep ="\t")
                ### BED file is 0 based. However, the outfile is 1 based
                start(REmap) <- start(REmap) - 1
		file.remove(outfile)
		REmap
}
