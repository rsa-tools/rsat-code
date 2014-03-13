################################################################
## Tutorial showing how to use "make" to manage a small project.
##
## Author: Jacques van Helden
## Date: 2014-03-13

################################################################
## The classical starting point: ask make to say "Hello World".
##
## To run this target, type the following command:
##    make -f makefiles/tutorial_make_project.mk hello_hello
##
## Since this is the first target of this file, the same result is
## obtained if no target if specified
##    make -f makefiles/tutorial_make_project.mk
##
## As you can see, the result is somewhat verbosy:"
## make first prints the command 'echo Hello World'"
## before running it."
hello_hello:
	echo "Hello World"

## Run a comment silently
##    make -f makefiles/tutorial_make_project.mk hello
##
## Commands preceded by @ are run silently"
## (the command is not displayed before running)."
hello:
	@echo "Hello World"

## Run another target
##    make -f makefiles/tutorial_make_project.mk bye
bye:
	@echo "Bye World"

## Run several targets subsequently. 
##
##    make -f makefiles/tutorial_make_project.mk hello bye
##
## Equivalent: declare a new target, followed by dependencies to
## other targets.
hello_bye: hello bye

##    make -f makefiles/tutorial_make_project.mk bye hello
bye_hello: bye hello

##    make -f makefiles/tutorial_make_project.mk bye bye
bye_bye: bye bye

################################################################
## Generate a random sequence.
##
## Use variables to specify some arguments (number of sequences and
## sequence lengths)
SEQ_NB=2
SEQ_LEN=50
randseq:
	random-seq -n ${SEQ_NB} -l ${SEQ_LEN}

## Specify the path of this makefile in a variable, in order to be
## able to call it from itself
MAKEFILE=${RSAT}/makefiles/tutorial_make_project.mk
MAKE=make -f ${MAKEFILE}

## Variables can be overwritten on the make command line
##    make -f makefiles/tutorial_make_project.mk randseq SEQ_NB=9
##
## The target below overwrites a variable in the same way , by calling
## make from within the makefile, with the adapted command-line.
##    make -f makefiles/tutorial_make_project.mk randseq_n9
randseq_n9:
	@${MAKE} randseq SEQ_NB=9




