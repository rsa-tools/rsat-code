################################################################
## Manage RSAT tasks on a PC cluster

#include ${RSAT}/makefiles/util.mk
include ${RSAT}/makefiles/init_RSAT.mk
MAKEFILE=${RSAT}/makefiles/cluster.mk

## IMPORTANT: FOR THE CLUSTER WE INSTALL THE TOOLS IN /usr/local/bin BECAUSE THE ARCHITECTURE IS NODE-DEPENDENT

################################################################
## Node list on the PC cluster of the BiGRe laboratory
AMD_1500=`seq -w 25 56 | awk '{print "n"$$1}' | grep -v n41 | grep -v n48 | xargs`
NODES_XEON1=n61 n62 n64
NODES_XEON2=`seq -w 67 75 | awk '{print "n"$$1}' | xargs`
AMD_2006=`seq -w 100 130 | awk '{print "n"$$1}' | grep -v n102 | grep -v n105 | grep -v n106 | grep -v n112 | grep -v n114 | grep -v n116 | grep -v n121 | grep -v 121 | grep -v n123 | grep -v 125 | xargs`
NODES = ${AMD_1500} ${NODES_XEON1} ${NODES_XEON2} ${AMD_2006}
#NODES = ${NODES_XEON1} ${NODES_XEON2}
list_nodes:
	@echo "NODES"
	@echo ${NODES}

################################################################
## Specific setting for compiling RSAT on the cluster@bigre
COMPILE_NODES=n29 n103 n70 n25
NODE_TASK=compile_all
compile_nodes:
	${MAKE} iterate_nodes NODES='${COMPILE_NODES}' NODE_CMD='cd ${RSAT}; make -f makefiles/init_RSAT.mk compile_all BIN=/usr/local/bin' 

## Run a command on a single node
NODE_CMD=hostname
one_node:
	@echo "NODE=${NODE}	NODE_CMD=${NODE_CMD}"
	@ssh ${NODE} '${NODE_CMD}'


## Run a task on a list of nodes
iterate_nodes:
	@for node in ${NODES} ; do \
		${MAKE} NODE=$${node} one_node ; \
	done


