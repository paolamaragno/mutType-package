vcffile <- system.file("extdata", "chr22.vcf.gz", package = "VariantAnnotation")
library(BSgenome.Hsapiens.UCSC.hg19)

test_that("Determine the mutation type of a set of single nucleotide variants in a genome", {

  expected_mut_types <- c("GC[C>T]CC", "GC[C>T]GG", "CC[C>T]GA", "CA[T>C]TG", "CC[T>C]GG", "TC[C>T]CT", "GC[C>A]TG",
                "GA[C>T]CA", "TT[C>T]CA", "GG[C>T]GA")

  expect_equal(mutation_type(vcffile, BSgenome.Hsapiens.UCSC.hg19, 5)[[1]][1:10], expected_mut_types)
})


test_that("Verify that the order of mutation types in the vector returned by mutation_type function is the same
          of the corresponding SNVs in the reordered GRanges object returned by the same function", {

  expected_mutgr <- GRanges(seqnames = rep('chr22', 5),
                    ranges = IRanges(start = c(50300085, 50300112, 50300165, 50300267, 50300437),
                                     end = c(50300087, 50300114, 50300167, 50300269, 50300439)),
                    REF = c('C', 'C', 'C', 'T', 'T'),
                    ALT = c('T', 'T', 'T', 'C', 'C'))

  expect_equal(mutation_type(vcffile, BSgenome.Hsapiens.UCSC.hg19, 3)[[2]][1:5], expected_mutgr)
})


test_that("Return error if context_length is even", {

  context_length <- 2

  expect_error(mutation_type(vcffile, BSgenome.Hsapiens.UCSC.hg19, context_length))
})


