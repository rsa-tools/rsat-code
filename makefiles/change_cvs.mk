################################################################
## Change the CVS configuration for all directories

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/change_cvs.mk
OLD_LOGIN=jvanheld
OLD_SERVER=cvs.bigre.ulb.ac.be
OLD_ROOT=/cvs/rsat
OLD_CONNECT=${OLD_LOGIN}\@${OLD_SERVER}:${OLD_ROOT}

NEW_LOGIN=rsat
NEW_SERVER=cvs.bigre.ulb.ac.be
NEW_ROOT=/cvs/rsat
NEW_CONNECT=${NEW_LOGIN}\@${NEW_SERVER}:${NEW_ROOT}

change_roots:
	find . -name Root -exec echo Modifying {} \;  -exec perl -pe "s|${OLD_CONNECT}|${NEW_CONNECT}|" -i {} \;

change_repositories:
	find . -name Repository -exec echo Modifying {} \;  -exec perl -pe "s|${OLD_ROOT}|${NEW_ROOT}|" -i {} \;

change_cvs: change_roots change_repositories