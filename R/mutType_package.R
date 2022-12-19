#' mutType - Determine the mutation type for a set of single nucleotide variants in a genome
#'
#' The \code{mutType} package - starting from a set of single nucleotide variants in VCF format, the corresponding
#' reference genome and a parameter \emph{context_length} specified by the user - determines for each mutation the
#' corresponding mutation type `UP[REF>ALT]DOWN` such that all mutation types have C or T as mutated reference base.\cr
#' The overall length of the mutation type is determined by the context_length parameter specified by the user.\cr
#' The package also provides a function which summarizes the mutation types for the set of mutations into
#' a count table.\cr
#' Eventually the package generates a pdf file showing the barplot visualization of all the mutation types with
#' a frequency higher than a threshold specified by the user.
#'
#' \tabular{ll}{
#' Package: \tab mutType\cr
#' Type: \tab Package\cr
#' Version: \tab 0.99.0\cr
#' Date: \tab 2022-12-12\cr
#' License: \tab GPL (>=2)\cr
#' }
#'
#' @name mutType-package
#' @aliases mutType-package mutType
#' @docType package
#' @author Paola Maragno [aut, cre]\cr
#' Politecnico di Milano\cr
#' Maintainer: Paola Maragno \cr
#' E-Mail: <paola.maragno@@mail.polimi.it>
NULL
