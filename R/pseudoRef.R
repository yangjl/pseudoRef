#' \code{pseudoRef}
#' Make a pseudo reference genome.
#'
#' @param fa Path for the reference fasta file. [string or DNAStringSet/DNAString object]
#' @param snpdt A data.table object with heterozygote SNPs coded with IUPAC ambiguity codes.
#'              [data.table, 4 required columns: chr, pos, ref, alt, (sample1, ..., sampleN)]
#' @param sidx A vector to indicate the sample columns. [vector, default=5:ncol(snpdt)].
#' @param arules Additional nucleotide substitution rules defined by users.
#'               [data.frame, 2 required columns: from, to, default=NULL]
#'               For example, arules <- data.frame(from=c("M", "Y", "R", "K"), to=c("C", "C", "G", "T")).
#' @param outdir Output directory. Sample specific sub-folders will be created. [string]
#'
#' @return A list of summary statistics of subsituted nucleotides. [list].
#'
#' @examples
#' # First of all, use BCFtools to convert VCF into IUPAC coded data.table:
#'
#' # bcftools view JRI20_filtered_snps_annot.bcf.gz -m2 -M2 -v snps -Oz -o JRI20_bi_snps_annot.vcf.gz
#' # bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%IUPACGT]\n' JRI20_bi_snps_annot.vcf.gz > JRI20_bi_snps_annot.txt
#' # bcftools query -f 'chr\tpos\tref\talt[\t%SAMPLE]\n' JRI20_bi_snps_annot.vcf.gz > JRI20_bi_snps_annot.header
#'
#' arules <- data.frame(from=c("M", "Y", "R", "K"), to=c("C", "C", "G", "T"))
#' res <- pseudoRef(fa, snpdt, sidx=5:24, arules, outdir)
#'
#' @export
pseudoRef <- function(fa, snpdt, sidx=5:ncol(snpdt), arules=NULL, outdir){

  ### load reference genome fasta file into DNAStringSet
  library("Biostrings")
  library("data.table")
  if(class(fa) == "character"){fa <- readDNAStringSet(filepath = fa, format="fasta")}
  if(class(fa) != "DNAStringSet" & class(fa) != "DNAString"){stop("fa should be a DNAStringSet or DNAString")}
  if(sum(names(arules) %in% c("from", "to")) !=2)
  {stop("arules should be a data.frame with at least two columns of from and to")}
  if(sum(names(snpdt) %in% c("chr", "pos", "ref", "alt")) !=4)
  {stop("snpdt should be a data.table with at least four columns of [chr, pos, ref, alt, (sample1,..., sampleN)]")}
  tab0 <- alphabetFrequency(fa)
  wd0 <- width(fa)
  nm0 <- names(fa)

  ### chr id, get only the first element of the blank seperated vector.
  chrid <-  gsub(" .*", "", names(fa))
  res <- list()
  ### substitute one sample at a time
  for(s in sidx){
    myfa <- fa
    #tab0 <- table(snpdt[, s, with=FALSE])
    sub <- snpdt[, c(1:4, s), with=FALSE]
    sid <- names(snpdt)[s]
    names(sub)[5] <- "SAMPLE"
    ## remove missing sites and sites that are the same as ref
    sub <- sub[SAMPLE != "./." & ref != SAMPLE]

    ### replace for each chrom
    for(i in 1:length(chrid)){
      subchr <- sub[chr == chrid[i]]
      if(nrow(subchr) > 0){

        idx <- which(chrid[i] == chrid)
        ## replaced base-pairs according to arules
        if(!is.null(arules)){
          for(k in 1:nrow(arules)){
            subchr[SAMPLE == arules$from[k], SAMPLE := arules$to[k]]
          }
        }
        myfa[[idx]] <- replaceLetterAt(x=myfa[[idx]], at=subchr[, pos], letter=subchr[, SAMPLE],
                                       if.not.extending="replace", verbose=FALSE)
      }
    }

    ### write for each sample
    tab1 <- alphabetFrequency(myfa)
    wd1 <- width(myfa)
    nm1 <- names(myfa)
    if(sum(wd0 - wd1) != 0 & sum(nm0 != nm1) != 0)
      {stop(sprintf("!!! changed chr length or chr id for sample [ %s ]", sid))}

    dir.create(paste0(outdir, "/", sid), showWarnings = FALSE)
    message(sprintf("###>>> printing fasta file for sample: [ %s ] ...", sid))
    writeXStringSet(myfa, paste0(outdir, "/", sid, "/", sid, ".fasta"), format="fasta", append=FALSE)

    ###output reprot
    res[[sid]] <- tab1 - tab0
  }
  return(res)
}
