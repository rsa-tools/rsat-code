---
title: "[RSAT](RSAT_home.cgi) - compare-reg-var manual"
output:
  html_document:
    toc: yes
    toc_depth: 3
  pdf_document:
    toc: yes
    toc_depth: 3
css: course.css
---

## NAME

compare-snps

## VERSION

$program\_version

## DESCRIPTION

Compare two or more files containing regulatory SNPs (predicted or
annotated).

Regulatory SNPs can be obtained from a variety of sources:

- _variation-scan_

    RSAT program that predicts rSNPs by scanning SNP sequences with
    position-specific scoring matrices.

    [http://www.rsat.eu/](http://www.rsat.eu/)

- _is-rsnp_

    A Web site allowing to predict rSNPs or to retrieve pre-computed
    predictions.

    [http://bioinformatics.research.nicta.com.au/software/is-rsnp/](http://bioinformatics.research.nicta.com.au/software/is-rsnp/)

- _haploreg_

    A web site combining annotations and predictions of regulatory SNPs,
    as well as integration of SNPs with regulatory information (ChIP-seq
    for transcription factors, histone modifications, ...).

    [http://www.broadinstitute.org/mammals/haploreg/haploreg.php](http://www.broadinstitute.org/mammals/haploreg/haploreg.php)

## AUTHORS

- Jacques.van-Helden\\@univ-amu.fr
- Yvon Mbouamboua

## CATEGORY

- variations

## USAGE

compare-snps \[-i inputfile\] \[-o outputfile\] \[-v #\] \[...\]

## INPUT FORMAT

## OUTPUT FORMAT

The output format is a tab-delimited file, with one row per rSNP.

A regulatory SNP is defined by 2 mandatory fields:

- _SNP ID_

    The identifier of the SNP shown or predicted to affect regulation.

- _TFBM_

    Identifier of the transcription factor binding motif (TFBM) affected
    by the SNP.

## SEE ALSO

- **variation-scan**

    The program _compare-rsnps_ takes as input the results of
    _variation-scan_.

## WISH LIST

- **wish 1**
- **wish 2**

## OPTIONS

- **-v #**

    Level of verbosity (detail in the warning messages during execution)

- **-h**

    Display full help message

- **-help**

    Same as -h

- **-i inputfile**

    If no input file is specified, the standard input is used.  This
    allows to use the command within a pipe.

- **-o outputfile**

    If no output file is specified, the standard output is used.  This
    allows to use the command within a pipe.
