"""
NAME
        matrix-scan-quick

VERSION
        200909

AUTHOR
        Matthieu Defrance <defrance@bigre.ulb.ac.be>

DESCRIPTION
        Faster and limited version of matrix-scan.

CATEGORY
        sequences
        pattern matching
        PSSM
        
USAGE        
        matrix-scan-quick -m matrix [-i sequences] [-bgfile background] [-h | --help]

INPUT FORMATS
  Sequence file
    Only sequences in FASTA format are supported.

  Matrix file
    Only the tab format is supported.
    see convert-matrix for details.

  Background file
    Only the INCLUSive format is supported.
    see convert-background-model for details.

OUTPUT FORMAT
    The output is a tab-delimited file, with one row per match.

SCORING SCHEME
    See matrix-scan -h for details.

ARGUMENTS
    -h, --help            show this help message and exit.

    -i #                  read sequence from filename # (FASTA format).
                          if not specified, the standard input is used.

    -o #                  print the output to filename #.
                          if not specified, the standard output is used.

    -m #                  read the matrix # (must be in tab format).
 
    -bgfile #             use # as background model (must be in INCLUSive format).
                          by default an equiprobable model is used.

    -2str                 scan both DNA strands.

    -1str                 scan only one DNA strand.

    -t #                  only capture site with a score >= #.

    -return distrib       output the weight score distribution.

    -return site          output the list of sites (default).

    -decimals #           precision parameter for the -return distrib option

    -pseudo #             pseudo-count for the matrix (1.0 by default)
"""
#    -V, --version         print version number and exit.

doc = globals()['__doc__']
for line in doc.splitlines():
    print '"%s\\n"' % line