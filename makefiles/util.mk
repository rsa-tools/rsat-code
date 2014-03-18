################################################################
## Some utilities

## Load site-specific options for the cluster + other parameters
include ${RSAT}/RSAT_config.mk

################################################################
## Variables
V=1
MAKEFILE=makefile
MAKE=make -s -f ${MAKEFILE}
DATE=`date +%Y-%M-%d_%H:%M:%S`
SSH_OPT = -e ssh 
RSYNC_OPT= -ruptvlz  ${SSH_OPT} 
RSYNC = rsync  ${RSYNC_OPT}
WGET=wget --passive-ftp -np -rNL

################################################################
## List of targets
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

HOUR=`date +%Y%m%d_%Hh`
JOB_DIR=`pwd`/jobs/${HOUR}
JOB_PREFIX=job
#JOB=`mktemp ${JOB_DIR}/${JOB_PREFIX}.XXXXXX`
JOB=`mktemp ${JOB_PREFIX}.XXXXXX`

command_queue:
	${MAKE} command_queue_${QUEUE_MANAGER}

## Send a jobs to a cluster using the torque quee management system
##
## A unique name is obtained for the script file with the command
## mktemp.  I use an awful trick to ensure that this name is generated
## only once for this target, by running a "for" loop with a single
## element (JOB). Without this script, mktemp would be called once for
## creating the script, and another time when sending it to the queue
## (it would thus have a different name, and the task would fail).
command_queue_torque:
	@mkdir -p ${JOB_DIR}
	@for job in ${JOB} ; do	\
		echo "Job ${JOB_DIR}/$${job}" ;	\
		echo "echo running on node "'$$HOST' > ${JOB_DIR}/$${job}; \
		echo "hostname" >> ${JOB_DIR}/$${job}; \
		echo "${MY_COMMAND}" >> ${JOB_DIR}/$${job} ;	\
		chmod u+x ${JOB_DIR}/$${job} ;	\
		qsub -m a -q ${QUEUE} -N $${job} -d ${PWD} -o ${JOB_DIR}/$${job}.log -e ${JOB_DIR}/$${job}.err ${QSUB_OPTIONS} ${JOB_DIR}/$${job} ;	\
		rm $${job} ;\
	done

## Send a jobs to a cluster using the SGE queue management system
##
## A unique name is obtained for the script file with the command
## mktemp.  I use an awful trick to ensure that this name is generated
## only once for this target, by running a "for" loop with a single
## element (JOB). Without this script, mktemp would be called once for
## creating the script, and another time when sending it to the queue
## (it would thus have a different name, and the task would fail).
command_queue_sge:
	@mkdir -p ${JOB_DIR}
	@echo "job dir	${JOB_DIR}"
	@for job in ${JOB} ; do	\
		echo "job	$${job}" ; \
		echo "Job ${JOB_DIR}/$${job}" ;	\
		echo "echo running on node "'$$HOST' > ${JOB_DIR}/$${job}; \
		echo "hostname" >> ${JOB_DIR}/$${job}; \
		echo "echo Job started" >> ${JOB_DIR}/$${job}; \
		echo "date" >> ${JOB_DIR}/$${job}; \
		echo "${MY_COMMAND}" >> ${JOB_DIR}/$${job} ;	\
		echo "echo Job done" >> ${JOB_DIR}/$${job}; \
		echo "date" >> ${JOB_DIR}/$${job}; \
		chmod u+x ${JOB_DIR}/$${job} ;	\
		qsub -m a -q ${QUEUE} -N $${job} -cwd -o ${JOB_DIR}/$${job}.log -e ${JOB_DIR}/$${job}.err ${QSUB_OPTIONS} ${JOB_DIR}/$${job} ; \
		rm $${job} ;\
	done

command_now:
	${MY_COMMAND}


################################################################
## Watch the number of jobs in the cluster queue
watch_jobs: watch_jobs_${QUEUE_MANAGER}

watch_jobs_torque:
	@hostname
	@date
	@echo "`qstat | grep -v '^---'| grep -v '^Job id' | wc -l`	Jobs"
	@echo "`qstat | grep ' R ' | wc -l`	Running"
	@echo "`qstat | grep ' Q ' | wc -l`	Queued" 
	@echo "`qstat | grep ' C ' | wc -l`	Completed"

watch_jobs_sge:
	@hostname
	@date
	@echo "`qstat | grep -v '^---'| grep -v '^job-ID' | wc -l`	Jobs"
	@echo "`qstat | grep ' r ' | wc -l`	Running"
	@echo "`qstat | grep ' Eqw ' | wc -l`	Errors" 


################################################################
## Iterate over all organisms
iterate_organisms:
	@echo "Iterating task ${ORG_TASK} over organisms"
	@echo "	${ORGANISMS}"
	@for org in ${ORGANISMS} ; do	\
		echo "" ; \
		echo "Organism $${org}" ; \
		${MAKE} ${ORG_TASK} ORG=$${org} ;	\
	done


################################################################
## Iterate over all oligo lengths
OLIGO_LENGTHS=1 2 3 4 5 6 7 8
iterate_oligo_lengths:
	@echo "Iterating task ${OLIGO_TASK} over oligonucleotide lengths ${OLIGO_LENGTHS}"
	@echo
	@for ol in ${OLIGO_LENGTHS} ; do	\
		${MAKE} ${OLIGO_TASK} OL=$${ol} ;	\
	done
