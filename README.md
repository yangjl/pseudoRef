# pseudoRef
It is an `R` package to make a pseudo reference genome subsituted with the SNP variants for given samples. 


## Installation
Install [devtools](https://github.com/hadley/devtools) first, and then use devtools to install imputeR from github.
```
devtools::install_github("yangjl/pseudoRef")
```
It requires two other papckages:
```
library("data.table")
library("Biostrings")
```
## Usage

You can find help documentation by simply typing `?pseudoRef`, which is the major function have been implemented in this package.
```
Usage:

     pseudoRef(fa, snpdt, sidx = 5:ncol(snpdt), arules = NULL, outdir)

Arguments:

      fa: Path for the reference fasta file. [string or
          DNAStringSet/DNAString oject]

   snpdt: A data.table object with heterozygote SNPs coded with IUPAC
          ambiguity codes. [data.table, 4 required columns: chr, pos,
          ref, alt, (sample1, ..., sampleN)]

    sidx: A vector to indicate the sample columns. [vector,
          default=5:ncol(snpdt)].

  arules: Additional nucleotide substitution rules defined by users.
          [data.frame, 2 required columns: from, to, default=NULL] For
          example, arules <- data.frame(from=c("M", "Y", "R", "K"),
          to=c("C", "C", "G", "T")).

  outdir: Output directory. Sample specific sub-folders will be
          created. [string]
```

## Note

Before running the package, we should use `BCFtools` to convert VCF/BCF file into IUPAC coded tab seperated file:
```
# bcftools view JRI20_filtered_snps_annot.bcf.gz -m2 -M2 -v snps -Oz -o JRI20_bi_snps_annot.vcf.gz
# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%IUPACGT]\n' JRI20_bi_snps_annot.vcf.gz > JRI20_bi_snps_annot.txt
# bcftools query -f 'chr\tpos\tref\talt[\t%SAMPLE]\n' JRI20_bi_snps_annot.vcf.gz > JRI20_bi_snps_annot.header
```

## TODO
1. add `test` data.
2. write a help document and real usage examples.
3. report a more reasonable substitution statistics.
4. visulize the report.
5. handle indels.

## License

This package is free and open source software, licensed under GPL.
