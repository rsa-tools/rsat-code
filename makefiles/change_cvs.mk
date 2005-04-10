
MAKEFILE=${RSAT}/makefiles/change_cvs.mk

include ${RSAT}/makefiles/util.mk

OLD_LOGIN=jvanheld
OLD_SERVER=cvs.scmbb.ulb.ac.be
#OLD_SERVER=rubens.ulb.ac.be
#OLD_ROOT=/rubens/dsk2/cvs
OLD_ROOT=/cvs
OLD_CONNECT=${OLD_LOGIN}\@${OLD_SERVER}:${OLD_ROOT}

NEW_ROOT=/cvs
NEW_LOGIN=rsat
NEW_SERVER=cvs.scmbb.ulb.ac.be
NEW_CONNECT=${NEW_LOGIN}\@${NEW_SERVER}:${NEW_ROOT}

change_roots:
	find . -name Root -exec echo Modifying {} \;  -exec perl -pe 's|${OLD_CONNECT}|${NEW_CONNECT}|' -i {} \;

change_repositories:
	find . -name Repository -exec echo Modifying {} \;  -exec perl -pe 's|${OLD_ROOT}|${NEW_ROOT}|' -i {} \;
