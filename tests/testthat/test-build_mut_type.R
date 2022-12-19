library(BSgenome.Hsapiens.UCSC.hg19)

test_that("Build the final mutation type", {

  mutgr <- GRanges(seqnames = rep('chr22', 5),
                   ranges = IRanges(start = c(50300077,50300100,50300525,50300526,50300561),
                                  end = c(50300081,50300104,50300529,50300530,50300565)),
                   REF = c('A','G','G','G','A'),
                   ALT = c('G','A','A','A','G')
                  )

  expected <- c("CT[C>C]TG", "GG[C>T]CG", "TG[C>T]CG", "CT[G>T]CC", "TG[T>C]TG")

  expect_equal(build_mut_type(mutgr, BSgenome.Hsapiens.UCSC.hg19, 5, rev=TRUE), expected)
})
