# REDseq (dev)

### Install
```
library(devtools)
install_github("hukai916/REDseq")
```

### Version
#### v1.45.3
  - Added `maskN` argument for `buildREmap` and `searchPattern` to allow restricting the searching space to non-N-blocks regions in genome. Otherwise, hits might be dominated by N-mers. `maskN` default to `TRUE` meaning that N-block regions will be excluded from searching.
  ```
  library(REDseq)
  library(BSgenome.Hsapiens.UCSC.hg19)
  seqnames <- seqname(Hsapiens)

  # For demo, search against the second chromosome only:
  chr_2 <- seqnames[[2]]
  REpatternFilePath = system.file("extdata", "examplePattern.fa", package="REDseq")

  # setting maskN = FALSE will give the same results as of before v1.45.3
  res <- buildREmap(REpatternFilePath, BSgenomeName=Hsapiens, outfile=tempfile(), chr = chr_2, maskN = FALSE)
  ```

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
