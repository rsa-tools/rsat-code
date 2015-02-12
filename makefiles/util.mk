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
DAY=`date +%Y%m%d`
TIME=`date +%Y%m%d_%H%M%S`
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

################################################################
## Send a jobs to a cluster using either SGE or torque as queue
## management system.  The command is written in a shell script
## (stored in the JOB dir), which is submitted to qsub.
##
## Note: qsub is used by several queue management systems (SGE,
## torque, ...), but the options are slightly different. The targets
## should be adapted to use it with another queue manager. 
## 
## A unique name is obtained for the script file with the command
## mktemp.  I use an awful trick to ensure that this name is generated
## only once for this target, by running a "for" loop with a single
## element (JOB). Without this script, mktemp would be called once for
## creating the script, and another time when sending it to the queue
## (it would thus have a different name, and the task would fail).
JOB_DIR=`pwd`/jobs/${DAY}
JOB_PREFIX=job
JOB=`mktemp -u ${JOB_PREFIX}.XXXXXX`
command_queue:
	${MAKE} command_queue_${QUEUE_MANAGER}

## Print the job in a script file
##
## A unique name is obtained for the script file with the command
## mktemp.  I use an awful trick to ensure that this name is generated
## only once for this target, by running a "for" loop with a single
## element (JOB). Without this script, mktemp would be called once for
## creating the script, and another time when sending it to the queue
## (it would thus have a different name, and the task would fail).
command_queue_print_job:
	echo "job	${JOB}"
	echo "Job file ${JOB_DIR}/${JOB}.sh"
	echo "echo running on node "'$$HOST' > ${JOB_DIR}/${JOB}.sh
	echo "hostname" >> ${JOB_DIR}/${JOB}.sh
	echo "echo Job started" >> ${JOB_DIR}/${JOB}.sh
	echo "date" >> ${JOB_DIR}/${JOB}.sh
	echo "${MY_COMMAND}" >> ${JOB_DIR}/${JOB}.sh
	echo "echo Job done" >> ${JOB_DIR}/${JOB}.sh
	echo "date" >> ${JOB_DIR}/${JOB}.sh
	chmod u+x ${JOB_DIR}/${JOB}.sh

## Send a jobs to a cluster using the torque quee management system
QSUB_CMD_TORQUE=qsub -m a -q ${CLUSTER_QUEUE} -N $${job} -d ${PWD} -o ${JOB_DIR}/$${job}.log -e ${JOB_DIR}/$${job}.err ${QSUB_OPTIONS} ${JOB_DIR}/$${job}.sh
command_queue_torque:
	@mkdir -p ${JOB_DIR}
	@for job in ${JOB} ; do	\
		${MAKE} command_queue_print_job JOB=$${job}; \
		${QSUB_CMD_TORQUE}; \
		echo "	Qsub command:	${QSUB_CMD_TORQUE}"; \
		echo "	Job file: 	${JOB_DIR}/$${job}.sh"; \
		echo "	Log file: 	${JOB_DIR}/$${job}.log"; \
		echo "	Error log:	${JOB_DIR}/$${job}.err"; \
	done
#		echo; \
#		echo "Enqueued command:	${MY_COMMAND}" ;	\
#		echo "echo running on node "'$$HOST' > ${JOB_DIR}/$${job}.sh; \
#		echo "hostname" >> ${JOB_DIR}/$${job}.sh; \
#		echo "${MY_COMMAND}" >> ${JOB_DIR}/$${job}.sh ;	\
#		chmod u+x ${JOB_DIR}/$${job}.sh ;	\
#		echo "	Qsub command:	${QSUB_CMD_TORQUE}" ;	\
#		echo "	Job file: 	${JOB_DIR}/$${job}.sh" ;	\
#		echo "	Log file: 	${JOB_DIR}/$${job}.log" ;	\
#		echo "	Error log:	${JOB_DIR}/$${job}.err" ;	\


## Send a jobs to a cluster using the SGE queue management system
command_queue_sge:
	@mkdir -p ${JOB_DIR}
	@echo "job dir	${JOB_DIR}"
	@for job in ${JOB} ; do	\
		${MAKE} command_queue_print_job JOB=$${job}; \
		qsub -m a -q ${CLUSTER_QUEUE} -N $${job} -cwd -o ${JOB_DIR}/$${job}.log -e ${JOB_DIR}/$${job}.err ${QSUB_OPTIONS} ${JOB_DIR}/$${job}.sh ; \
	done

## Send a jobs to batch
command_queue_batch:
	@mkdir -p ${JOB_DIR}
	@echo "job dir	${JOB_DIR}"
	@for job in ${JOB} ; do	\
		${MAKE} command_queue_print_job JOB=$${job}; \
		batch -f ${JOB_DIR}/$${job}.sh ; \
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
