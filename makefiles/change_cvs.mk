
include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/change_cvs.mk

################################################################
## Change the CVS configuration for all directories

OLD_LOGIN=`logname`
OLD_SERVER=www.bigre.ulb.ac.be
OLD_ROOT=/cvs/rsat
OLD_CONNECT=${OLD_LOGIN}\@${OLD_SERVER}:${OLD_ROOT}

NEW_LOGIN=`logname`
NEW_SERVER=cvs.bigre.ulb.ac.be
NEW_ROOT=/cvs/rsat
NEW_CONNECT=${NEW_LOGIN}\@${NEW_SERVER}:${NEW_ROOT}

change_roots:
	find . -name Root -exec echo Modifying {} \;  -exec perl -pe "s|${OLD_CONNECT}|${NEW_CONNECT}|" -i {} \;

change_repositories:
	find . -name Repository -exec echo Modifying {} \;  -exec perl -pe "s|${OLD_ROOT}|${NEW_ROOT}|" -i {} \;

change_cvs: change_roots change_repositories