################################################################
## template.mk makefile
## Usage: make -f template.mk

MAKEFILE=${RSAT}/makefiles/template.mk
MAKE = make -sk -f ${MAKEFILE}

### tags
usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}
