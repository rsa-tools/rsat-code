################################################################
## Some utilities

## level of verbosity
V=1

################################################################
### commands
MAKE=make -s -f ${MAKEFILE}

SSH_OPT = -e ssh 
RSYNC_OPT= -ruptvlz  ${SSH_OPT} 
RSYNC = rsync  ${RSYNC_OPT}

WGET=wget --passive-ftp -np -rNL

################################################################
#### list of targets
usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}

################################################################
## Send a task, either direcoy (WHEN=now) or in a queue (WHEN=queue)
## for a cluster.
WHEN=now
JOB_DIR=jobs
my_command:
	@echo ${WHEN} command ${MY_COMMAND}
	${MAKE} command_${WHEN}

command_queue:
	@mkdir -p ${JOB_DIR}
	@for job in ${JOB} ; do						\
		echo "Job $${job}" ;					\
		echo "${MY_COMMAND}" > $${job} ;		\
		qsub -m e -q rsa@merlin.ulb.ac.be -N $${job} -j oe	\
			-o $${job}.log $${job} ;	\
	done

command_now:
	${MY_COMMAND}


################################################################
## Iterate over all organisms
iterate_organisms:
	@echo "Iterating task ${ORG_TASK} over organisms"
	@echo ${ORGANISMS}
	@for org in ${ORGANISMS} ; do			\
		${MAKE} ${ORG_TASK} ORG=$${org} ;	\
	done

