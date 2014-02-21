################################################################
## This file is a template : please adapt the path [RSAT_PARENT_PATH]
## to your local configuration (should be the directory in which the
## rsat folder is installed.

################################################################ 
## Configuration for Regulatory Sequence Analysis Tools (RSAT)
export RSAT=[RSAT_PARENT_PATH]/rsat
export PATH=${RSAT}/bin:${PATH}
export PATH=${RSAT}/perl-scripts:${PATH}
export PATH=${RSAT}/perl-scripts/parsers:${PATH}
export PATH=${RSAT}/python-scripts:${PATH}

################################################################
## Class path for metabolic pathway analysis tools
if  [ ${CLASSPATH} ]; then
       export CLASSPATH=${CLASSPATH}:${RSAT}/java/lib/NeAT_javatools.jar
else
       export CLASSPATH=.:${RSAT}/java/lib/NeAT_javatools.jar
fi

################################################################
## Default path for the Ensembl Perl modules
export PERL5LIB=${RSAT}/lib/ensemblgenomes-21-74/ensembl/modules::${PERL5LIB}
export PERL5LIB=${RSAT}/lib/ensemblgenomes-21-74/ensembl-compara/modules::${PERL5LIB}
export PERL5LIB=${RSAT}/lib/ensemblgenomes-21-74/ensembl-external/modules::${PERL5LIB}
export PERL5LIB=${RSAT}/lib/ensemblgenomes-21-74/ensembl-functgenomics/modules::${PERL5LIB}
export PERL5LIB=${RSAT}/lib/ensemblgenomes-21-74/ensembl-tools/modules::${PERL5LIB}
export PERL5LIB=${RSAT}/lib/ensemblgenomes-21-74/ensembl-variation/modules::${PERL5LIB}
export PERL5LIB=${RSAT}/lib/bioperl-release-1-2-3/bioperl-live::${PERL5LIB}
