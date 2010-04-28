################################################################
## Manage RSAT tasks on a PC cluster

#include ${RSAT}/makefiles/util.mk
include ${RSAT}/makefiles/init_RSAT.mk
MAKEFILE=makefiles/cluster.mk


################################################################
## Node list on the PC cluster of the BiGRe laboratory
NODES_XEON1=n61 n62 n64
NODES_XEON2=`seq -w 67 75 | awk '{print "n"$$1}' | xargs`
AMD_2006=`seq -w 100 130 | awk '{print "n"$$1}' | xargs`
AMD_1500=`seq -w 25 57 | awk '{print "n"$$1}' | xargs`
NODES = ${AMD_1500} ${NODES_XEON1} ${NODES_XEON2} ${AMD_2006}
list_nodes:
	@echo "NODES"
	@echo ${NODES}

################################################################
## Specific setting for compiling RSAT on the cluster@bigre
COMPILE_NODES=n103 n71 n70 n25
NODE_TASK=compile_all
compile_nodes:
	${MAKE} iterate_nodes NODES='${COMPILE_NODES}' NODE_CMD='cd ${RSAT}; make -f makefiles/init_RSAT.mk compile_all'


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


