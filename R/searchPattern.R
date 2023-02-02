searchPattern <-function (dict0,  BSgenomeName, outfile = "", chr = c("all"), maskN = TRUE) {
  if (missing(dict0) || class(dict0) != "DNAStringSet") {
      stop("dictO is required as a DNAStringSet object!")
  }

	if (missing(BSgenomeName) || class(BSgenomeName) != "BSgenome")
		{
			stop("BSgenomeName is required as a BSgenome object!")
		}
    if (file.exists(outfile)) {
        stop("outfile specified as ", outfile, " already exists! Please rename the outfile!")
    }

    seqnames <- seqnames(BSgenomeName)

    if (chr != c("all")) {
      if (all(chr %in% seqnames)) {
        seqnames <- chr
      } else {
        stop("Not all chr names are valid! Please ensure chr is a subset of seqnames(BSgenomeName)!")
      }
    }

    append <- FALSE
    for (seqname in seqnames) {
        subject <- BSgenomeName[[seqname]]
        if (maskN) {
          subject <- maskMotif(subject, "N")
        }
        cat(">>> Finding all hits in sequences", seqname, "...\n")
        for (i in seq_len(length(dict0))) {
            patternID <- names(dict0)[i]
            if (length(patternID) < 1) {
                patternID = paste("pattern", i, sep = "")
            }
            pattern <- dict0[[i]]
            plus_matches <- matchPattern(pattern, subject, fixed = FALSE)
	    if (length(plus_matches)>0)
	   {
		 names(plus_matches) <- rep.int(paste(patternID, "+",sep="."), length(plus_matches))
           	 writeHits(seqname, plus_matches, "+", file = outfile, append = append)
            	append <- TRUE
	  }
            rcpattern <- reverseComplement(pattern)
            if (rcpattern != pattern) {
                minus_matches <- matchPattern(rcpattern, subject,fixed = FALSE)
		 if (length(minus_matches)>0)
           	{
                	names(minus_matches) <- rep.int(paste(patternID,"-", sep="."), length(minus_matches))
                	writeHits(seqname, minus_matches, "-", file = outfile, 	append = append)
			append <- TRUE
		}
            }
        }
        cat(">>> DONE searching\n")
    }
}
