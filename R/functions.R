#' \strong{Build the final mutation type}
#'
#' \code{build_mut_type} function takes a GRanges object containing a set of single nucleotide variants retrieved from a VCF file,
#' the corresponding reference genome, the parameter \emph{context_length} specified by the user and the indication whether
#' the reverse complement of the sequence containing the SNV must be computed.\cr
#' This function returns a vector with the mutation types `UP[REF>ALT]DOWN` of all the variants in the GRanges
#' object, such that all mutation types have C or T as mutated reference base.\cr
#' The overall length of the mutation type is determined by the context_length parameter.\cr
#' This function is not exported to be used by the user but it is defined only to be exploited in the \link{mutation_type}
#' function.
#'
#' @usage build_mut_type(mutgr, reference_genome, context_length, rev)
#' @param mutgr GRanges object containing a set of single nucleotide variants
#' @param reference_genome Reference genome
#' @param context_length Parameter specified by the user to indicate the overall length of the mutation type - it is
#' used to compute how many nucleotides upstream and downstream the SNV base include in the mutation type
#' @param rev Parameter to indicate whether the reverse complement of the sequence containing the SNV must be computed
#' (default is FALSE)
#' @return Vector with the mutation types `UP[REF>ALT]DOWN` of all the variants in the GRanges object with C or T as
#' mutated reference base.
#' @author Paola Maragno\cr Politecnico di Milano\cr Maintainer: Paola Maragno\cr E-Mail: <paola.maragno@@mail.polimi.it>
#'
#' @import Biostrings
#' @import BSgenome

build_mut_type <- function(mutgr, reference_genome, context_length, rev = FALSE) {

  # if the mutated reference base is C or T
  if (rev == FALSE) {

    # for each SNV contained in the GRanges object get the corresponding sequence from the specified reference genome
    seq <- as.character(getSeq(reference_genome, mutgr))

    # mutation type UP[REF>ALT]DOWN is obtained by pasting subparts of the sequence with C or T as mutated reference base
    # and the corresponding alternative nucleotide
    paste0(substr(seq, 1, (context_length-1)/2), '[',
           substr(seq, ((context_length-1)/2)+1, ((context_length-1)/2)+1), '>',
           mutgr$ALT, ']', substr(seq, ((context_length+1)/2)+1, context_length))

  } # if the mutated reference base is A or G
  else {

    # for each SNV contained in the GRanges object get the corresponding reverse complement sequence from the
    # specified reference genome
    seq <- as.character(reverseComplement(getSeq(reference_genome, mutgr)))

    # mutation type UP[REF>ALT]DOWN is obtained by pasting subparts of the sequence with C or T as mutated reference base
    # and the reverse complement of the corresponding alternative nucleotide
    paste0(substr(seq, 1, (context_length-1)/2), '[',
           substr(seq, ((context_length-1)/2)+1, ((context_length-1)/2)+1), '>',
           as.character(reverseComplement(DNAStringSet(mutgr$ALT))),
           ']', substr(seq, ((context_length+1)/2)+1, context_length))
}}


#' \strong{Determine the mutation type of a set of single nucleotide variants in a genome}
#'
#' \code{mutation_type} function takes a set of single nucleotide variants in VCF format, the corresponding reference
#' genome and a parameter \emph{context_length} specified by the user - that must be odd - and determines for each mutation
#' the corresponding mutation type `UP[REF>ALT]DOWN` such that all mutation types have C or T as mutated reference base.\cr
#' The overall length of the mutation type is determined by the context_length parameter.
#'
#' @usage mutation_type(vcf_file, reference_genome, context_length)
#' @param vcf_file VCF file containing a set of single nucleotide variants (like `chr22.vcf`)
#' @param reference_genome Reference genome
#' @param context_length Parameter specified by the user to indicate the overall length of the mutation type - it is
#' used to compute how many nucleotides upstream and downstream the SNV base include in the mutation type
#' @return List of two elements: the mutation types `UP[REF>ALT]DOWN` vector with C or T as mutated reference base
#' and the GRanges object containing first the SNVs with C or T as mutated reference base and then the SNVs with G or A as
#' mutated reference base. In this way the position of each mutation type in the vector is the same of the corresponding SNV in
#' the GRanges object returned by the function.
#' @author Paola Maragno\cr Politecnico di Milano\cr Maintainer: Paola Maragno\cr E-Mail: <paola.maragno@@mail.polimi.it>
#' @examples
#'
#' if (!'BSgenome.Hsapiens.UCSC.hg19' %in% installed.packages()) {
#' if (!requireNamespace("BiocManager", quietly = TRUE)) {
#' install.packages("BiocManager")
#' }
#' BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
#' }
#'
#' library('BSgenome.Hsapiens.UCSC.hg19')
#' vcf_file <- system.file("extdata", "chr22.vcf.gz", package = "VariantAnnotation")
#' results <- mutation_type(vcf_file, BSgenome.Hsapiens.UCSC.hg19, 3)
#' mut_types <- results[[1]]
#' mutgr_reordered <- results[[2]]
#' mut_types[1:10]
#' mutgr_reordered[1:10]
#'
#' @import VariantAnnotation
#' @import GenomicRanges
#' @import IRanges
#' @import Biostrings
#' @export

mutation_type <- function(vcf_file, reference_genome, context_length) {

  # read the VCF file with the reference genome specified in reference_genome parameter
  vcf <- readVcf(vcf_file, unique(genome(reference_genome)))

  # raise an error in case the context_length parameter is even
  if (context_length %% 2 == 0) {
    stop('Error: The context_length parameter must be odd!')
  }

  # create a GRanges object containing for each SNV the chromosome on which it is, the position on the chromosome
  # - including also a number of nucleotides before and after the SNV such that the overall sequence has a length
  # equal to the context_length value - the reference and alternative nucleotides
  mutgr <- GRanges(if (grepl('chr', levels(seqnames(vcf))[1])) {
                  seqnames = as.character(seqnames(vcf))} else {
                  seqnames = paste0('chr', as.character(seqnames(vcf)))},
                  ranges = IRanges(start = start(vcf) - (context_length-1)/2, end = end(vcf) + (context_length-1)/2),
                  REF = as.character(ref(vcf)),
                  ALT = as.character(unlist(alt(vcf))))

  # remove from the GRanges object all the variants for which the reference or alternative base is more than one
  # since we are only interested in single nucleotide variants
  mutgr <- mutgr[nchar(mutgr$REF) == 1 & nchar(mutgr$ALT) == 1]

  # split the GRanges object into two GRanges objects one containing the SNVs with C or T as mutated reference base and the
  # other containing the SNVs with G or A as mutated reference base
  mutgr_CT <- mutgr[mutgr$REF %in% c('C','T')]
  mutgr_AG <- mutgr[mutgr$REF %in% c('A','G')]

  # use the build_mut_type function of mutType package to compute the vector of mutation types UP[REF>ALT]DOWN
  # of all the variants in the two GRanges objects.
  # rev parameter is used to indicate whether the reverse complement of the sequence containing the SNV
  # must be computed so that by the end all mutation types have C or T as mutated reference base
  mut_type_CT <- build_mut_type(mutgr_CT, reference_genome, context_length)
  mut_type_AG <- build_mut_type(mutgr_AG, reference_genome, context_length, rev = TRUE)

  # return the concatenation of the two vectors of mutation types and the concatenation of the two GRanges objects
  mut_types <- c(mut_type_CT, mut_type_AG)
  mutgr_reordered <- append(mutgr_CT, mutgr_AG)
  results <- list(mut_types, mutgr_reordered)
  return(results)
}


#' \strong{Count table of mutation types}
#'
#' \code{count_table} function summarizes the mutation types for the set of mutations into a count table
#' reporting the number of mutations per mutation type.
#'
#' @usage count_table(mut_types)
#' @param mut_types Vector of mutation types `UP[REF>ALT]DOWN` returned by \link{mutation_type} function
#' @return Count table reporting the number of mutations per mutation type.
#' @author Paola Maragno\cr Politecnico di Milano\cr Maintainer: Paola Maragno\cr E-Mail: <paola.maragno@@mail.polimi.it>
#' @examples
#'
#' mut_types <- c("GT[C>A]CA", "GG[C>G]TC", "CC[T>C]TC", "GT[C>T]GT", "TA[C>A]CG", "GT[C>A]CA",
#' "CA[C>G]CT", "CC[T>C]CT", "CA[T>C]AT", "TT[C>G]TC", "CC[T>C]CT", "CT[T>C]TG", "CC[T>C]CT", "TA[C>A]CG")
#' count_table(mut_types)
#'
#' @export

count_table <- function(mut_types) {

  # convert into data frame the result of table function on the vector of mutation types
  # UP[REF>ALT]DOWN returned by mutation_type function
  c_table <- as.data.frame(table(mut_types))

  return(c_table)
}


#' \strong{Graphical visualization of the mutation types with a frequancy higher than a threshold}
#'
#' \code{graphical_summary} function generates a pdf file showing the barplot visualization of all the mutation
#' types with a frequency higher than a threshold specified by the user.
#'
#' @usage graphical_summary(c_table, freq, file_name)
#' @param c_table Count table returned by \link{count_table} function that summarizes the mutation types contained in a VCF file
#' @param freq Threshold of frequency that a mutation type must have at least in order to be visualized in the barplot
#' @param file_name Name of the pdf file in which the function will plot its graphical output
#' @return The function returns the name of the pdf file - that is named as chosen by the user and it is stored in the working
#' directory - showing the barplot visualization of the frequencies of the mutation types with a frequency at least equal to
#' \emph{freq} threshold.
#' @author Paola Maragno\cr Politecnico di Milano\cr Maintainer: Paola Maragno\cr E-Mail: <paola.maragno@@mail.polimi.it>
#' @examples
#'
#' c_table <- data.frame(mut_types = c("AA[C>A]AA", "AA[C>A]AC", "AA[C>A]AG", "AA[C>A]AT", "AA[C>A]CA", "AA[C>A]CC",
#'  "AA[C>A]CT", "AA[C>A]GA", "AA[C>A]GC", "AA[C>A]GG", "AA[C>A]GT", "AA[C>A]TA", "AA[C>A]TC"),
#'  Freq = c(15,  13,  32,  44,  25,  5,  67,  21,  14,  42,  52,  21,  19))
#' graphical_summary(c_table, 30, "Mut_types_visualization_30.pdf")
#'
#' @import ggplot2
#' @import grDevices
#' @export
#'

graphical_summary <- function(c_table, freq, file_name) {

  # extract from c_table the mutation types with frequency at least equal to the freq parameter
  c <- c_table[c_table$Freq >= freq,]

  # retrieve the mutation types column from the filtered count table
  x <- c$mut_types

  # generate a factor for the subsequent splitting of the initial data frame into different data frames each
  # containing at most 100 mutations
  f <- gl(ceiling(nrow(c)/100), 100, length=nrow(c))

  # open a pdf file with the name indicated by the user
  pdf(file_name)

  # for each iteration a new data frame is created with at most 100 mutations used to build the barplot, then
  # each barplot is plotted in a different page of the pdf file
  tapply(x, f, FUN = function(sublist){c2 <- c[which(c$mut_types %in% sublist),]
                                print(ggplot(c2, aes(x = mut_types, y = Freq)) +
                                        geom_bar(stat = "identity") +
                                        geom_col(aes(fill = Freq)) +
                                        scale_fill_gradient2(low = 'yellow',
                                                             mid ='red',
                                                             high = 'darkgreen',
                                                             midpoint = max(c2$Freq)/2) +
                                        theme_classic() +
                                        theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 5)) +
                                        labs(title = paste("Frequency of each mutation type with Freq >", as.character(freq)),
                                             subtitle = paste('Overall number of mutation types with Freq >', as.character(freq), ':',  as.character(nrow(c)),
                                                              '\nOverall number of mutation types with Freq >', as.character(freq), 'in this page:',  as.character(nrow(c2))),
                                             x = "Mutation type",
                                             y = "Frequency"))
  })

  # close the pdf file
  dev.off()
  # return the name of the pdf file
  return(file_name)
}
