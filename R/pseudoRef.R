#' \code{get_file2tab}
#' Scan txt files and extract information into a data.frame.
#'
#' @param files Full names of the files to scan [vector].
#' @param features Unique features for grep to search [vector].
#'
#' @return return Extracted values [data.frame].
#'
#' @examples
#' files <- list.files(path = dir, pattern = fileptn, full.names = TRUE)
#' features <- c("C methylated in CHH context")
#' get_file2tab(files, features, replace=T )
#'
#' @export
pseudo_ref <- function(snpdf, fa, sampleidx=5:24, sub_rules){
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
