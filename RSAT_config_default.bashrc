################################################################
## This file is a template : please adapt the path [RSAT_PARENT_PATH]
## to your local configuration (should be the directory in which the
## rsat folder is installed.

################################################################ 
## Configuration for Regulatory Sequence Analysis Tools (RSAT)
export RSAT=[RSAT_PARENT_PATH]/rsat
export RSAT_SITE=localhost
export RSAT_WWW=http://localhost/rsat
export RSAT_WS=http://localhost/rsat
export PATH=${RSAT}/bin:${PATH}
export PATH=${RSAT}/perl-scripts:${PATH}
export PATH=${RSAT}/perl-scripts/parsers:${PATH}
export PATH=${RSAT}/python-scripts:${PATH}

################################################################
## Make sure R libs are accessible to the apache user.
##
## Note by JvH: installing R libraries in the RSAT package creates problems 
## because some libraries are installed both at the system-level and in the RSAT folder. 
## I thus comment this line, I will have to find a better solution. 
#
# if  [ ${R_LIBS_USER} ]; then
#   export R_LIBS_USER=${R_LIBS_USER}:${RSAT}/R-scripts/Rpackages/
#e lse
#  export R_LIBS_USER=${RSAT}/R-scripts/Rpackages/
# fi

################################################################
## Add RSAT python libraries to PYTHONPATH
if  [ ${PYTHONPATH} ]; then
    export PYTHONPATH=${PYTHONPATH}:${RSAT}/lib/python2.7/site-packages/
else
    export PYTHONPATH=${RSAT}/lib/python2.7/site-packages/
fi

################################################################
## Class path for metabolic pathway analysis tools
if  [ ${CLASSPATH} ]; then
       export CLASSPATH=${CLASSPATH}:${RSAT}/java/lib/NeAT_javatools.jar
else
       export CLASSPATH=.:${RSAT}/java/lib/NeAT_javatools.jar
fi

################################################################
## Use ssh as remote shell for CVS (required to install Ensembl API)
export CVS_RSH=ssh

################################################################
## Default path for the Ensembl Perl modules and sofwtare tools
export ENSEMBL_RELEASE=87
export ENSEMBLGENOMES_RELEASE=34
export PATH=${RSAT}/lib/ensemblgenomes-${ENSEMBLGENOMES_RELEASE}-${ENSEMBL_RELEASE}/ensembl-git-tools/bin:${PATH}
export PERL5LIB=${RSAT}/lib/bioperl-release-${BIOPERL_VERSION}/bioperl-live::${PERL5LIB}
export PERL5LIB=${RSAT}/lib/ensemblgenomes-${ENSEMBLGENOMES_RELEASE}-${ENSEMBL_RELEASE}/ensembl/modules::${PERL5LIB}
export PERL5LIB=${RSAT}/lib/ensemblgenomes-${ENSEMBLGENOMES_RELEASE}-${ENSEMBL_RELEASE}/ensembl-compara/modules::${PERL5LIB}
export PERL5LIB=${RSAT}/lib/ensemblgenomes-${ENSEMBLGENOMES_RELEASE}-${ENSEMBL_RELEASE}/ensembl-external/modules::${PERL5LIB}
export PERL5LIB=${RSAT}/lib/ensemblgenomes-${ENSEMBLGENOMES_RELEASE}-${ENSEMBL_RELEASE}/ensembl-functgenomics/modules::${PERL5LIB}
export PERL5LIB=${RSAT}/lib/ensemblgenomes-${ENSEMBLGENOMES_RELEASE}-${ENSEMBL_RELEASE}/ensembl-tools/modules::${PERL5LIB}
export PERL5LIB=${RSAT}/lib/ensemblgenomes-${ENSEMBLGENOMES_RELEASE}-${ENSEMBL_RELEASE}/ensembl-variation/modules::${PERL5LIB}
export PERL5LIB=${RSAT}/lib/ensemblgenomes-${ENSEMBLGENOMES_RELEASE}-${ENSEMBL_RELEASE}/ensemblgenomes-api/modules::${PERL5LIB}
export PERL5LIB=${RSAT}/lib/biomart-perl/lib/::${PERL5LIB}


################################################################
## Python libraries
## (required for external tool CEAS)
export PYTHONPATH=$PYTHONPATH:${RSAT}/lib/python2.7/site-packages
