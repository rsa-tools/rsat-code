################################################################
## Some utilities

## level of verbosity
V=1

################################################################
### commands
MAKEFILE=makefile
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
my_command:
	@echo "${WHEN} command ${MY_COMMAND}"
	${MAKE} command_${WHEN}

JOB_DIR=jobs
JOB=${JOB_PREFIX}`mktemp ${JOB_DIR}/job.XXXXXX`
command_queue:
	@mkdir -p ${JOB_DIR}
	@for job in ${JOB} ; do						\
		echo "Job $${job}" ;					\
		echo "echo running on node "'$$HOST' > $${job}; \
		echo "${MY_COMMAND}" >> $${job} ;		\
		qsub -m e -q rsa@merlin.scmbb.ulb.ac.be -N $${job} -j oe	\
			-o $${job}.log $${job} ;	\
	done

command_now:
	${MY_COMMAND}


################################################################
## Iterate over all organisms
iterate_organisms:
	@echo "Iterating task ${ORG_TASK} over organisms"
	@echo "	${ORGANISMS}"
	@echo
	@for org in ${ORGANISMS} ; do			\
		${MAKE} ${ORG_TASK} ORG=$${org} ;	\
	done


################################################################
## Iterate over all organisms
iterate_oligo_lengths:
	@echo "Iterating task ${OLIGO_TASK} over oligonucleotide lengths ${OLIGO_LENGTHS}"
	@echo
	@for ol in ${OLIGO_LENGTHS} ; do			\
		${MAKE} ${OLIGO_TASK} OL=$${ol} ;	\
	done
