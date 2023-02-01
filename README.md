# REDseq (dev)

### Install
```
library(devtools)
install_github("hukai916/REDseq")
```

### Version
#### v1.45.2
  - Added `chr` argument for `buildREmap` and `searchPattern` to restrict the searching space to reduce memory cost.
  ```
  library(REDseq)
  library(BSgenome.Celegans.UCSC.ce2)
  seqnames <- seqname(Celegans)

  # Search against the first chromosome only:
  chr_1 <- seqnames[[1]]

  REpatternFilePath = system.file("extdata", "examplePattern.fa", package="REDseq")
  res <- buildREmap(REpatternFilePath, BSgenomeName=Celegans, outfile=tempfile(), chr = chr_1)

  ```
  - Optimized the conversion from BED to GRanges by switching to `rtracklayer::import`. 
