## ----style, echo = FALSE, results = 'asis'------------------------------------
library(BiocStyle)

## ---- echo = FALSE------------------------------------------------------------
library(knitr)

## ---- echo = FALSE------------------------------------------------------------
library(mutType)

## -----------------------------------------------------------------------------
library('BSgenome.Hsapiens.UCSC.hg19')

vcf_file <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
results <- mutation_type(vcf_file, BSgenome.Hsapiens.UCSC.hg19, 5)
mut_types <- results[[1]]
mutgr_reordered <- results[[2]]
mut_types[1:10]
mutgr_reordered[1:10]

## -----------------------------------------------------------------------------
mut_types <- c("GT[C>A]CA", "GG[C>G]TC", "CC[T>C]TC", "GT[C>T]GT", "TA[C>A]CG", "GT[C>A]CA", "CA[C>G]CT", "CC[T>C]CT", "CA[T>C]AT", "TT[C>G]TC", "CC[T>C]CT", "CT[T>C]TG", "CC[T>C]CT", "TA[C>A]CG")
count_table(mut_types)

## -----------------------------------------------------------------------------
c_table <- data.frame(mut_types = c("AA[C>A]AA", "AA[C>A]AC", "AA[C>A]AG", "AA[C>A]AT", "AA[C>A]CA", "AA[C>A]CC",
                                     "AA[C>A]CT", "AA[C>A]GA", "AA[C>A]GC", "AA[C>A]GG", "AA[C>A]GT", "AA[C>A]TA", 
                                     "AA[C>A]TC"), 
                      Freq = c(15,  13,  32,  44,  25,  5,  67,  21,  14,  42,  52,  21,  19))
graphical_summary(c_table, 30, "Mut_types_visualization_30.pdf")

