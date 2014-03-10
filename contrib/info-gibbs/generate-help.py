"""
NAME
        info-gibbs

VERSION
        2008

AUTHOR
        Matthieu Defrance <defrance@scmbb.ulb.ac.be>

DESCRIPTION
        Gibbs sampling algorithm for motifs discovery.
        Searches for highly conserved motifs in a set of DNA sequences.
        Convervation is based on the motif information content (Relative Entropy).

CATEGORY
        sequences
        motif discovery
        
USAGE        
        info-gibbs -l motiflength [-i inputfile] [-h | --help]

ARGUMENTS
  GENERAL OPTIONS
    -h, --help            show this help message and exit

    -v #, --verbosity=#   set verbosity to level #
                              0 no verbosity
                              1 low verbosity
                              2 high verbosity
                              3 maximal verbosity + debug + trace

    -i #, --input=#       read sequence from # (must be in FASTA format)
                          if not specified, the standard input is used

    -l #, --length=#      set oligomer length to #
                          when the option dyad is used # represents the length monad1 + monad2
                          EXAMPLE: --length=7

    --maxspacing=#        set maximal spacing between motif monad to # (only for dyadic motif).
                          in this case the parameter -l corresponds to length of monad1 + monad2 (without spacing)

    --minspacing=#        set minimal spacing between motif monad to # (only for dyadic motif).
                          in this case the parameter -l corresponds to length of monad1 + monad2 (without spacing)

    -s #, --strand=#      search in foward strand + or in both strands +-
                          DEFAULT: +-
                          EXAMPLE: --strand=+

    -n #, --iter=#        maximum number of Gibbs sampling iterations
                          DEFAULT: 1000

    -w #, --words=#       number of motif occurrences that are expected to
                          be found (incompatible with -e)
                          DEFAULT: 10

    -e #, --expected=#    expected number of motif occurrences per sequence
                          that are expected to be found (incompatible with -w)
                          DEFAULT: 1

    -m #, --motifs=#      number of motifs to extract (one by default)
                          DEFAULT: 1

    -b #, --bgfile=#      use # predefined INCLUSive background model
                          [http://homes.esat.kuleuven.be/~thijs/help/help_motifsampler.html#background]    
                          EXAMPLE --bgfile=mybgfile

    -d #, --dmin=#        set minimal distance between 2 motif occurrences to #

    -t #                  set the temperature (should be in range [0.6 1.4])
                          DEFAULT: 1.0

    -r #  --nrun=#        try to run the Gibbs sampling seach # times
                          DEFAULT: 10

    --finalcycle          try to collect N best sites based on their weight score

    --rseed=#             set random seed to #

    --title=#             Add title # to output

    -V, --version         print version

"""
# -o #, --output=#      output results to #
#                       if not specified, the standard output is used

doc = globals()['__doc__']
for line in doc.splitlines():
    print '"%s\\n"' % line