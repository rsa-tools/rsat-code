
.SUFFIXES: .tex .dvi .log .aux .toc .pdf .ps .view .xdvi

############################################################
#
# $Id: makefile,v 1.3 2004/03/11 09:13:32 jvanheld Exp $
#
# Time-stamp: <2003-10-02 15:13:34 jvanheld>
#
############################################################

MAKE=make -s
RM=rm -f

TEX_FILES = \
	test.tex \
	install_guide.tex \
	tutorial_shell_rsat.tex

INCLUDE_TEX_FILES = \
	rsat_latex_commands.tex

AUX_FILES=${TEX_FILES:.tex=.aux} ${INCLUDE_TEX_FILES:.tex=.aux}
BBL_FILES=${TEX_FILES:.tex=.bbl}
BLG_FILES=${TEX_FILES:.tex=.blg}
DVI_FILES=${TEX_FILES:.tex=.dvi}
IDX_FILES=${TEX_FILES:.tex=.idx}
ILG_FILES=${TEX_FILES:.tex=.ilg}
IND_FILES=${TEX_FILES:.tex=.ind}
LOG_FILES=${TEX_FILES:.tex=.log}
TOC_FILES=${TEX_FILES:.tex=.toc}

PDF_FILES=${TEX_FILES:.tex=.pdf}
PS_FILES=${TEX_FILES:.tex=.ps}
VIEW_FILES=${TEX_FILES:.tex=.view}

usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' makefile

DOC=tutorial_shell_rsat
it:
	${MAKE} fullclean
	${MAKE} ${DOC}.xdvi
	${MAKE} tex_clean

ug:
	${MAKE} tutorial_shell_rsat.pdf
	${MAKE} tex_clean

_publish:
	@echo "${PUBLISH_TARGET}"
	@echo "Publishing to ${PUBLISH_SOURCE} => ${PUBLISH_TARGET}"
	@if [ -d $(dir ${PUBLISH_TARGET}) ] ; then \
		cp --remote-destination ${PUBLISH_SOURCE} ${PUBLISH_TARGET} ; \
	else \
		echo "	Error: missing directory. Please create it manually"; \
	fi

all_pdf: ${PDF_FILES}

force_pdf: 
	${MAKE} pdf_clean
	${MAKE} all_pdf

all_ps: ${PS_FILES}

force_ps: 
	${MAKE} ps_clean
	${MAKE} all_ps

all_view: ${VIEW_FILES}


############################################################
#
# Cleaning targets
#

tex_clean:
	${RM} ${AUX_FILES} ${BBL_FILES} ${BLG_FILES} ${DVI_FILES} ${IDX_FILES} ${ILG_FILES} ${IND_FILES} ${LOG_FILES} ${TOC_FILES}

pdf_clean:
	${RM} ${PDF_FILES}

ps_clean:
	${RM} ${PS_FILES}

fullclean:
	${MAKE} tex_clean
	${MAKE} pdf_clean
	${MAKE} ps_clean

############################################################
#
# Rules
#

.tex.dvi:
	latex $< 
	makeindex $< 
	latex $< 

.dvi.pdf:
	dvipdf $<

.dvi.ps:
	dvips -o $@ $<

.pdf.view:
	acroread $<

.dvi.xdvi:
	xdvi $<
