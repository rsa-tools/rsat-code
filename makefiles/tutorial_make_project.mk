################################################################
## Tutorial showing how to use "make" to manage a small project.
##
## Author: Jacques van Helden
## Date: 2014-03-13

## Contents
##
## 1. Hello world
## 2. Simple example: generating a random sequence
## 3. Conditional execution (if)
## 4. Iterating a target ("for" loop)


################################################################
## Specify the path of this makefile in a variable, in order to be
## able to call the file from itself.
MAKEFILE=${RSAT}/makefiles/tutorial_make_project.mk
MAKE=make -f ${MAKEFILE}

################################################################
## List of targets can be obtained by running this command:
##
##    make -f ${RSAT}/makefiles/tutorial_make_project.mk list_targets
##
## Note: it is convenient to place this target in first position in
## the makefile: when make is called without arguments, the target
## list is displayed.
##
##    make -f ${RSAT}/makefiles/tutorial_make_project.mk
##
list_targets:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}


################################################################
##                                                            ##
##                 Section 1: Hello world                     ##
##                                                            ##
################################################################

## The classical starting point: ask make to say "Hello World".
##
## To run this target, type the following command:
##    make -f ${RSAT}/makefiles/tutorial_make_project.mk hello_hello
##
## As you can see, the result is somewhat verbosy:"
## make first prints the command 'echo Hello World'"
## before running it."
##
## You can use the option "-s" to run make silently
##    make -s -f ${RSAT}/makefiles/tutorial_make_project.mk hello_hello
hello_hello:
	echo "Hello World"

## A more refined way to manage the verbosity is to indicate that a
## specific command sould run silently. Commands preceded by @ are
## run silently" (the command is not displayed before running)."
##
##    make -f ${RSAT}/makefiles/tutorial_make_project.mk hello
##
## This allows to display some commands and not other, in the same
## target.
hello:
	@echo "Hello World"

## Run another target
##    make -f ${RSAT}/makefiles/tutorial_make_project.mk bye
bye:
	@echo "Bye World"

## Run several targets subsequently. 
##
##    make -f ${RSAT}/makefiles/tutorial_make_project.mk hello bye
##
## Equivalent: declare a new target, followed by dependencies to
## other targets.
hello_bye: hello bye

##    make -f ${RSAT}/makefiles/tutorial_make_project.mk bye hello
bye_hello: bye hello

## If the same target is calle twice, it only runs once, because make
## manages dependencies between tasks. 
##
##    make -f ${RSAT}/makefiles/tutorial_make_project.mk bye bye
bye_bye: bye bye



################################################################
##                                                            ##
##  Section 2. Simple example: generating a random sequence   ##
##                                                            ##
################################################################

## Generate a random sequence. We first define some variables to
## specify some arguments (number of sequences and sequence lengths)
##
##    make -f ${RSAT}/makefiles/tutorial_make_project.mk randseq
SEQ_NB=2
SEQ_LEN=50
randseq:
	random-seq -n ${SEQ_NB} -l ${SEQ_LEN}

################################################################
## Variables can be overwritten on the make command line
##    make -f ${RSAT}/makefiles/tutorial_make_project.mk randseq SEQ_NB=9
##
## The target below overwrites a variable in the same way , by calling
## make from within the makefile, with the adapted command-line.
##    make -f ${RSAT}/makefiles/tutorial_make_project.mk randseq_n9
randseq_n9:
	@${MAKE} randseq SEQ_NB=9



################################################################
##                                                            ##
##         Section 3: Conditional execution (if)              ##
##                                                            ##
################################################################

## Create a directory for storing random sequences.  First test if the
## directory exists. If not, create it and display a message.
RANDSEQ_DIR=results/randseq
randseq_dir:
	@if [ -d "${RANDSEQ_DIR}" ] ; then \
		echo ""; \
	else \
		echo "Creating directory to store random sequences"; \
		mkdir -p ${RANDSEQ_DIR} ; \
		echo "	RANDSEQ_DIR	${RANDSEQ_DIR}" ; \
	fi

################################################################
## Generate a random sequence and store it to a file. By default, the
## repetition is set to 1.
##
## In the header of this target, we indicate that it depdends from the
## target randseq_dir.
##
## After command execution, we report the result file name.
REP=01
RAND_SEQ_FILE=${RANDSEQ_DIR}/randseq_L${SEQ_LEN}_N${SEQ_NB}_rep${REP}.fasta
RAND_SEQ_CMD=random-seq -n ${SEQ_NB} -l ${SEQ_LEN} -o ${RAND_SEQ_FILE}
randseq_to_file: randseq_dir
	${RAND_SEQ_CMD}
	@echo "	${RAND_SEQ_FILE}"

################################################################
##                                                            ##
##       Section 4: Task iteration ("for" loop)               ##
##                                                            ##
################################################################

## Random repetitions We use Perl to define repeat numbers having the
## same number of digits. This is convenient, since alphabetical and
## numerical orders are then identical.
RAND_REP_FROM=0
RAND_REP_TO=20
RAND_REP_DIGITS=02
RAND_REPEATS=`perl -le 'for $$i (${RAND_REP_FROM}..${RAND_REP_TO}) {printf "%${RAND_REP_DIGITS}d\n", $$i}' | xargs`
list_rep_numbers:
	@echo "Random repetitions from ${RAND_REP_FROM} to ${RAND_REP_TO}"
	@echo "	${RAND_REPEATS}"


## Run a task iteratively using a "for" loop
ITERATING_TASK=randseq_to_file
iterate_randseq:
	@echo
	@echo "Iterating task ${ITERATING_TAKS}	from ${RAND_REP_FROM} to ${RAND_REP_TO}"
	@for i in ${RAND_REPEATS}; do \
		${MAKE} -s ${ITERATING_TASK} REP=$$i; \
	done



