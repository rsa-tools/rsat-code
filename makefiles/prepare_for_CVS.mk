################################################################
#
# prepare a directory for conversion from RCS to CVS : 
# - check-in all files
# - move the ,v files one dir up from the RCS directory

## questions
##- how to specify the omitted files ?

MAKE=make -s -f ${MAKEFILE}
RSYNC=rsync -ruptvl
RSAT=${HOME}/rsa-tools

MAKEFILE=prepare_for_CVS.mk


## default directory
DIR=perl-scripts

## all directories to submit
DIRS= perl-scripts				\
	public_html				\
	makefiles				\
	doc

SUBDIRS=`find ${DIR} -type d | grep -v RCS | grep -v SWISS | grep -v perllib`

## excluded directories and files
EXCLUDED=					\
	--exclude obsolete			\
	--exclude data				\
	--exclude tmp				\
	--exclude oldies			\
	--exclude logs				\
	--exclude auto				\
	--exclude papers			\
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


all_tasks: 
	${MAKE} rsync_one_dir 
	${MAKE} treat_one_dir
	${MAKE} iterate_subdirs TASK=treat_one_dir 

rsync_one_dir:
	${RSYNC} ${EXCLUDED} ${RSAT}/${DIR} .

treat_one_dir: find_obsolete_RCS remove_obsolete_RCS check_in clean_tmp


find_obsolete_RCS:
	@echo 
	@echo "Finding obsolete RCS fiels in dir ${DIR}"
	mkdir -p tmp
	find ${DIR}/RCS -maxdepth 1 -type f  | perl -pe 's/,v//' | perl -pe 's|RCS/||' | sort -u > tmp/in_RCS.txt
	find ${DIR} -maxdepth 1 -type f  | sort -u  > tmp/in_dir.txt
	wc tmp/in_RCS.txt tmp/in_dir.txt
	diff tmp/in_RCS.txt tmp/in_dir.txt | grep '^<' | wc
	diff tmp/in_RCS.txt tmp/in_dir.txt | grep '^>' | wc

remove_obsolete_RCS: find_obsolete_RCS
	@echo 
	@echo "Removing obsolete RCS fiels from dir ${DIR}"
	diff tmp/in_RCS.txt tmp/in_dir.txt			\
		| grep '^<'				\
		| perl -pe 's|${DIR}/|${DIR}/RCS/|'	\
		| perl -pe 's/< //'			\
		| perl -pe 's/(\S+)/$$1,v/'		\
		| xargs rm -f
clean_tmp:
	\rm -rf tmp


check_in:
	(cd ${DIR}; find . -maxdepth 1 -type f -exec ci -m'Before converting RCS to CVS' {} \;)


