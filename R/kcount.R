#' \code{kcount}
#' Counting kmer frequency for a given genome or pseudo genome.
#'
#' @param fa Path for the reference fasta file. [string or DNAStringSet/DNAString object]
#' @param feature A data.frame object (BED like, but all 1-based coordinates).
#'              [data.frame, 4 required columns: chr, start, end, geneid]
#' @param kmers A vector to indicate kmers to count. Note, k<=15. [vector, default=3:7].
#'
#' @return A data.frame of kmer as rows x geneid as columns. [data.frame].
#'
#' @examples
#'
#' feature <- data.frame(chr=1, start=seq(from=1, to=10000001, by=1000),
#'                       end=seq(from=1000, to=11000000, by=1000)[1:10001],
#'                       geneid=paste0("g", 1:10001))
#' kcount <- function(fa, feature, kmers=3:7)
#'
#' @export
kcount <- function(fa, feature, kmers=3:7){
  ### load reference genome fasta file into DNAStringSet
  library("Biostrings")
  #library("data.table")
  if(class(fa) == "character"){fa <- readDNAStringSet(filepath = fa, format="fasta")}
  if(class(fa) != "DNAStringSet" & class(fa) != "DNAString"){stop("fa should be a DNAStringSet or DNAString")}
  if(sum(names(feature) %in% c("chr", "start", "end", "geneid")) != 4)
  {stop("feature file should be a data.frame with three columns chr, start(1-based) and end(1-based)")}

  ### chr id, get only the first element of the blank seperated vector.
  chrid <-  gsub(" .*", "", names(fa))

  ### replace for each chrom
  allout <- data.frame()
  for(i in 1:10){
    subf <- subset(feature, chr==i)
    chridx <- which(chrid == i)

    ### for each kmer
    out <- data.frame()
    for(k in kmers){
      message(sprintf("###>>> counting [ %smer ] for [ chr%i ] ...", k, i))
      res <- sapply(1:nrow(feature), function(x){
        oligonucleotideFrequency(subseq(fa[[chridx]], start=subf$start[x], end=subf$end[x]),
                                 width=k, as.array=F)
      })
      out <- rbind(out, res)
      names(out) <- subf$geneid
    }
    allout <- cbind(out, allout)
  }
  return(allout)

}
