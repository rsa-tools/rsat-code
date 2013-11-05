################################################################
#
# Generate documentation for the Regulatory Sequence Analysis Tools

MAKEFILE=${RSAT}/makefiles/doc.mk
include ${RSAT}/makefiles/util.mk


PERL_DOC_DIR=doc/perl-scripts_doc
SCRIPTS=perl-scripts
INFILE=${FILE}
OUTFILE=${FILE}.html
FILE=${SCRIPTS}/lib/RSAT/matrix.pm
DIRNAME=`dirname ${FILE}` 
OUTDIR=${PERL_DOC_DIR}/${DIRNAME}
#INFILE=${SCRIPTS}
#OUTFILE=${PERL_DOC_DIR}
pod_html:
	find ${SCRIPTS} -name '*.pm' -exec ${MAKE} one_pod_html FILE={} \;
	find ${SCRIPTS} -name '*.pl' -exec ${MAKE} one_pod_html FILE={} \;
	find ${SCRIPTS} -name '*.lib' -exec ${MAKE} one_pod_html FILE={} \;

one_pod_html:
	@mkdir -p ${PERL_DOC_DIR}
	@mkdir -p ${OUTDIR}
	pod2html --netscape --index --recurse --verbose	\
		--htmlroot=${PERL_DOC_DIR}		\
		--podroot ${SCRIPTS}/lib/RSAT		\
		--title '${FILE}'			\
		--outfile ${PERL_DOC_DIR}/${OUTFILE}	\
		--infile ${INFILE}
