################################################################
#
# prepare a directory for conversion from RCS to CVS : 
# - check-in all files
# - move the ,v files one dir up from the RCS directory

## questions
##- how to specify the omitted files ?

MAKEFILE=makefiles/prepare_for_CVS.mk
MAKE=make -s -f ${MAKEFILE}
RSYNC=rsync -ruptvl
RSAT=${HOME}/rsa-tools



## default directory
DIR=perl-scripts

## all directories to submit
DIRS= perl-scripts				\
	public_html				\
	makefiles				\
	doc

SUBDIRS=`find ${TARGET_DIR}/${DIR} -type d | grep -v RCS | perl -pe 's|${TARGET_DIR}/||g'`

## excluded directories and files
EXCLUDED=					\
	--exclude perllib			\
	--exclude obsolete			\
	--exclude oldies			\
	--exclude old				\
	--exclude data				\
	--exclude tmp				\
	--exclude logs				\
	--exclude auto				\
	--exclude papers			\
	--exclude '*.pdf'			\
	--exclude '*.ps'			\
	--exclude '*~'


### tags
usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^(\S+):/){ print "\t$$1\n" }' ${MAKEFILE}

################################################################
#
#
TASK=list_subdirs
iterate_dirs:
	@for dir in ${DIRS}; do			\
		${MAKE} ${TASK} DIR=$${dir};		\
	done

list_subdirs:
	@echo ${SUBDIRS}

iterate_subdirs:
	@for dir in ${SUBDIRS}; do		\
		${MAKE} ${TASK} DIR=$${dir};	\
	done

all:
	${MAKE} iterate_dirs TASK=all_tasks_one_dir

all_tasks_one_dir: 
	${MAKE} rsync_one_dir 
	${MAKE} iterate_subdirs TASK=treat_one_dir 

TARGET_DIR=for_CVS
rsync_one_dir:
	mkdir -p ${TARGET_DIR}
	${RSYNC} ${EXCLUDED} ${RSAT}/${DIR} ${TARGET_DIR}

treat_one_dir: find_obsolete_RCS remove_obsolete_RCS check_in rcs_files_one_dir_up clean_tmp 

find_obsolete_RCS:
	@echo 
	@echo "Finding obsolete RCS fiels in dir ${TARGET_DIR}/${DIR}"
	mkdir -p ${TEMPORARY}
	find ${TARGET_DIR}/${DIR}/RCS -maxdepth 1 -type f  | perl -pe 's/,v//' | perl -pe 's|RCS/||' | sort -u > ${TEMPORARY}/in_RCS.txt
	find ${TARGET_DIR}/${DIR} -maxdepth 1 -type f  | sort -u  > ${TEMPORARY}/in_dir.txt
	wc ${TEMPORARY}/in_RCS.txt ${TEMPORARY}/in_dir.txt
	diff ${TEMPORARY}/in_RCS.txt ${TEMPORARY}/in_dir.txt | grep '^<' | wc
	diff ${TEMPORARY}/in_RCS.txt ${TEMPORARY}/in_dir.txt | grep '^>' | wc

remove_obsolete_RCS: find_obsolete_RCS
	@echo 
	@echo "Removing obsolete RCS fiels from dir ${TARGET_DIR}/${DIR}"
	diff ${TEMPORARY}/in_RCS.txt ${TEMPORARY}/in_dir.txt			\
		| grep '^<'				\
		| perl -pe 's|${TARGET_DIR}/${DIR}/|${TARGET_DIR}/${DIR}/RCS/|'	\
		| perl -pe 's/< //'			\
		| perl -pe 's/(\S+)/$$1,v/'		\
		| xargs rm -f
TEMPORARY=${TARGET_DIR}/tmp
clean_tmp:
	\rm -rf ${TEMPORARY}

check_in:
	(cd ${TARGET_DIR}/${DIR}; find . -maxdepth 1 -type f -exec ci -m'Before converting RCS to CVS' {} \;)

rcs_files_one_dir_up:
	mv ${TARGET_DIR}/${DIR}/RCS/*,v ${TARGET_DIR}/${DIR}/
	rmdir ${TARGET_DIR}/${DIR}/RCS/
