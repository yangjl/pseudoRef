#' \code{pseudoRef}
#' Make a pseudo reference genome.
#'
#' @param fa Path for the fasta file. [string]
#' @param snpdt A data.table object with heterozygote IUPAC SNP codes. [data.table]
#' @param sampleidx A vector to indicate the sample columns. [vector].
#' @param sub_rules Nucleotide substitution rules defined by users. [data.frame]
#'                  For example, sub_rules <- data.frame(from=c("M", "Y", "R", "K"), to=c("C", "C", "G", "T")).
#' @param dir Output directory. Sample specific sub-folders will be created. [string]
#'
#' @return Fasta files in specified directory under sample-specific sub-folders. [fasta].
#'
#' @examples
#' sub_rules <- data.frame(from=c("M", "Y", "R", "K"), to=c("C", "C", "G", "T"))
#'
#' features <- c("C methylated in CHH context")
#' get_file2tab(files, features, replace=T )
#'
#' @export
pseudoRef <- function(fa, snpdt, sampleidx=5:24, sub_rules){
  chrid <-  gsub(" .*", "", names(fa))

  ### get one sample at a time
  for(s in sampleidx){

    sub <- snpdf[, c(1:4, s), with=FALSE]
    sid <- names(snpdf)[s]
    names(sub)[5] <- "SAMPLE"
    sub <- sub[SAMPLE != "./." & ref != SAMPLE]

    ### replace for each chrom
    for(i in 1:length(chrid)){
      subchr <- sub[chr == chrid[i]]
      ## remove missing sites and sites that are the same as ref
      idx <- which(chrid[i] == chrid)

      ## replaced base-pairs based on sub_rules
      for(k in 1:nrow(sub_rules)){
        subchr[SAMPLE == sub_rules$from[k], SAMPLE := sub_rules$to[k]]
      }

      fa[[idx]] <- replaceLetterAt(x=fa[[idx]], at=subchr[, pos], letter=subchr[, SAMPLE],
                                   if.not.extending="replace", verbose=FALSE)
    }

    ### write for each sample
  }
}
