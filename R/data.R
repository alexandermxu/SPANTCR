#' Modified BLOSUM62 matrix.
#'
#' A matrix containing thresholded BLOSUM values to calculate distance between k-mers.
#'
#' @format A matrix with 22 rows and 22 columns:
#'
"BLOSUMCappedDiff"

#' Amino acid parameter matrix.
#'
#' A matrix containing properties of amino acids.
#'
#' @format A matrix with 22 rows and 6 columns:
#'
"AminoAcidFilter"

#' Vector of tick centroids 1-100 to use with single CDR3 genes for distance comparisons.
#'
#' A vector from 0.005 to 0.995.
#'
#' @format A vector of length 100:
#'
"FineBinTicks"

#' Vector of tick centroids 1-200 to use with paired CDR3 for distance comparisons.
#'
#' A vector from 0.005 to 1.995.
#'
#' @format A vector of length 200:
#'
"FineBinTicksPaired"

#' Sample data table of paired VDJDB CMV pp65 epitope data
#'
#' A data.table properly formatted for SPAN-TCR
#'
#' @format A data.table with 444 rows and 6 columns:
#'
"VDJDBCMVPairedInput"
