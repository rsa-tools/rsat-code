---
title: "[RSAT](RSAT_home.cgi) - crer-scan manual"
output:
  html_document:
    toc: yes
    toc_depth: 3
  pdf_document:
    toc: yes
    toc_depth: 3
css: course.css
---

## Description

This tool takes as input a set of "sites" (genomic coordinates), and reports cis-regulatory enriched regions (CRER), i.e. genomic intervals containing a higher number of sites than expected by chance.

The tool can be used to predict cis-regulatory modules (CRM) from collections of transcription factor binding sites obatined by various methods.

1. Predictions obtained by scanning sequences with position-specific scoring matrices (for example the output of matrix-scan
2. Peaks from ChIP-seq experiments.
3. Annotated binding sites imported from a transcription factor database.
4. Any other data source that produces a set of genomic features.

The program proceeds by evaluating the statistical enrichment in sites for each possible interval (window) encompassing at least two sites. Window size can be restricted to reduce computation time. 

Enrichment is estimated by computing the binomial p-value.

## Authors

- Marie Artufel <sup>*</sup>
- Lucie Khamvongsa <sup>*</sup>
- Jacques van Helden

## Options

### Sites

Sites are described in tab-delimited format, with one row per site and one column per attribute.
The order of the columns depends on the selected format (see below). 

**crer-scan** uses the following information to detect enrichment.
- sequence ID, 
- start position, 
- end position, 
- site score
- site sequence

### Site formats

Sites can be described in *bed* or *ft* formats.

- **bed** see the official [definition at UCSC](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
- **ft** RSAT feature format, produced for example by matrix-scan.

See [convert-features](help.convert-features.html) for more information about feature formats.

### Thresholds

Thresholds can be specified on different statistics in order to select relevant CRERs.



| Field | Description |
|------------------------|--------------------------|
| Site score  | Score indicated in the "score" column of the input file. |



## References

The principle of CRER detection and interpretation of the results are described in the protocol about ***matrix-scan***. 

- Turatsinze, J.-V., Thomas-Chollier, M., Defrance, M. and van Helden, J. (2008) Using RSAT to scan genome sequences for transcription factor binding sites and cis-regulatory modules. Nature Protocols, 3, 1578–1588.

