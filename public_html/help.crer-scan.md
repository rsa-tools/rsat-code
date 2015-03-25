---
title: "[RSAT](RSAT_home.cgi) - crer-scan manual"
output:
  html_document:
    highlight: tango
    theme: united
    toc: yes
    toc_depth: 3
  pdf_document:
    highlight: zenburn
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

Field | Description<sup>1</sup>
-----------------|---------------
Score of input sites |   minimal and maximal site score to be considered
P-value of input sites |	maximal p-value of sites to be considered.Recommended to be the higher site p-value considered. Default : 1e-4. A lower threshold is not available because not considered sites with lower p-value has no bological sense.	
CRER size | minimal and maximal size of the enriched region (in bp). Default: minimal site size : 30bp and maximal site size : 500bp
Sites per CRER | minimal and maximal number of sites covered by the enriched region. Default: minimal number of sites : 2
Inter-site distance (bp) | distance between successive sites to be considered. A minimal inter-site distance can be used to prevent overlap between redundant matrices. Default : minimal distance : 1 bp. And a maximal inter-site distance can be used to prevent merging distinct modules into a single one. Note: the maximal inter-site distance is one of the most influential parameters in cluster-buster. Default : maximal distance : 35 bp
CRER significance	| minimal binomial significance. Default: 2. An upper threshold is not available because not considered CRER with higher significance has no bological sense.  
Overlap between sites |  maximal number of sites covered by the enriched region. A lower threshold is not available because two distinct sites is caracterised by no overlap. 

<sup>1</sup>Note: None value is used to skip the threshold

### Output options

It possible to report sequence limits, only if provided in site file. Different mode are available such as:

Mode | Description
-----------------|---------------
None | No return any limits. Recommended when the input contains a lot of sequences. 
Filtered | return the limits filtered of the sequence. Only the sequence limits of CRERs.
All | return all limits of sequences scanned. 

## References

The principle of CRER detection and interpretation of the results are described in the protocol about ***matrix-scan***. 

- Turatsinze, J.-V., Thomas-Chollier, M., Defrance, M. and van Helden, J. (2008) Using RSAT to scan genome sequences for transcription factor binding sites and cis-regulatory modules. Nature Protocols, 3, 1578–1588.

