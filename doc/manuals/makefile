.SUFFIXES: .tex .dvi .log .aux .toc .pdf .ps .view .xdvi .html

############################################################
#
# $Id: makefile,v 1.19 2012/04/16 15:53:44 jvanheld Exp $
#
# Time-stamp: <2003-10-02 15:13:34 jvanheld>
#
############################################################

MAKE=make -s
RM=rm -f

#LATEX=pplatex
LATEX=latex

TEX_FILES = \
	install_guide.tex \
	tutorial_shell_rsat.tex \
	rsat_web_server.tex \
	neat_tutorial.tex \
	neat_web_server.tex 

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
OUT_FILES=${TEX_FILES:.tex=.out}
BRF_FILES=${TEX_FILES:.tex=.brf}

PDF_FILES=${TEX_FILES:.tex=.pdf}
PS_FILES=${TEX_FILES:.tex=.ps}
VIEW_FILES=${TEX_FILES:.tex=.view}

usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' makefile

DOC=tutorial_shell_rsat
it:
	@${MAKE} fullclean
	@${MAKE} ${DOC}.xdvi
	@${MAKE} tex_clean

rsat_tutorial:
	@${MAKE} tutorial_shell_rsat.pdf
	@${MAKE} tex_clean
	@echo tutorial_shell_rsat.pdf

install_guide:
	@${MAKE} install_guide.pdf
#	@${MAKE} tex_clean
	@echo install_guide.pdf

rsat_web_server:
	@${MAKE} rsat_web_server.pdf
	@${MAKE} tex_clean
	@echo rsat_web_server.pdf

neat_tutorial:
	@${MAKE} neat_tutorial.pdf
	@${MAKE} tex_clean
	@echo neat_tutorial.pdf

neat_web_server:
	@${MAKE} neat_web_server.pdf
	@${MAKE} tex_clean
	@echo neat_web_server.pdf

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
	@${MAKE} pdf_clean
	@${MAKE} all_pdf

all_ps: ${PS_FILES}

force_ps: 
	@${MAKE} ps_clean
	@${MAKE} all_ps

all_view: ${VIEW_FILES}


############################################################
#
# Cleaning targets
#

tex_clean:
	${RM} ${AUX_FILES} ${BBL_FILES} ${BLG_FILES} ${DVI_FILES} ${IDX_FILES} ${ILG_FILES} ${IND_FILES} ${LOG_FILES} ${TOC_FILES} ${OUT_FILES} ${BRF_FILES}

pdf_clean:
	${RM} ${PDF_FILES}

ps_clean:
	${RM} ${PS_FILES}

fullclean: tex_clean pdf_clean ps_clean emacs_clean

################################################################
## Clean the emacs temporary files
emacs_clean:
	find . -name '*~' -exec rm {} \;

############################################################
#
# Rules
#

SPLIT=4
.tex.html:
	latex2html -split ${SPLIT} $<
	@echo "$<.html"

.tex.dvi:
	${LATEX} $<; bibtex $<; ${LATEX} $<; ${LATEX} $<
	@echo "$<.dvi"
#	${LATEX} $< 
#	makeindex $< 
#	${LATEX} $< 

.dvi.pdf:
	dvipdf $<
	@echo "$<.pdf"

.dvi.ps:
	dvips -o $@ $<
	@echo "$<.ps"

.pdf.view:
	acroread $<

.dvi.xdvi:
	xdvi $<


################################################################
## Publish the pdf documents on the web site
DOC_SITE=rsat@rsat.ulb.ac.be:rsat/public_html/distrib/
#DOC_SITE2=jvanheld@itzamna.ccg.unam.mx:rsat/public_html/distrib/
publish:
	rsync -ruptvl -e ssh *.pdf ${DOC_SITE}
#	rsync -ruptvl -e ssh *.pdf ${DOC_SITE2}
