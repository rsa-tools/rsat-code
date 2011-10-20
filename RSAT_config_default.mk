################################################################
## Site-specific configuration for the RSAT makefiles
## 
## This file is automatically loade by the script
## ${RSAT}/makefiles/util.mk, which is itself loaded by all the RSAT
## makefile scripts.


## This parameter is essential, since the same program (qsub) has
## different parameters depending on the queue manager. Thus, a qsub
## command for torque will not be understood by sge and reciprocally;
##
## supported: torque | sge
##
QUEUE_MANAGER=sge

## Name of the queue where the jobs have to be sent
QUEUE=


