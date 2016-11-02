source 00_config.bash

cd ${RSAT}; source RSAT_config.bashrc ## Reload the (updated) RSAT environment variables


################################################################
## Installation of Perl modules required for RSAT
################################################################
## Notes
##
## 1) limxml2-dev is required to compile the Perl module XML::LibXML
# apt-get install limxml2-dev 
## 2) For some modules, installation failed until I used "force"
##	 force install SOAP::WSDL
##	 force install SOAP::Lite
##
## 3) Problem of dependency when installint REST::Client, even with "force".
##	 MCRAWFOR/REST-Client-291.tar.gz              : make_test FAILED but failure ignored because 'force' in effect
##	 NANIS/Crypt-SSLeay-0.64.tar.gz               : make NO
## To solve it I found this
## 	apt-get install libnet-ssleay-perl
##	apt-get install libcrypt-ssleay-perl
##	apt-get install libssl-dev


## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## NOTE: at this stage I had to inactivate the Kaspersky antivirus because it
## prevents CPAN from downloading and installing Perl modules !
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

cpan YAML ## Answer "yes" !!!
cpan CPAN ## Update CPAN

## Set working directory to RSAT
cd ${RSAT}

## Get the list of Perl modules to be installed
make -f makefiles/install_rsat.mk  perl_modules_list

## Check which Perl modules are already installed
make -f makefiles/install_rsat.mk perl_modules_check
## The result file name will be displayed at the end of the tests

## Check the result of perl modules test
# cat check_perl_modules_eval.txt
## On Ubuntu 14.04, Object::InsideOut has status "Fail" but there is
## apparently no problem

## Check missing modules
# grep Fail  check_perl_modules_eval.txt

## Identify Perl modules that were not OK after the ubuntu package installation
#grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;'
MISSING_PERL_MODULES=`grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;' | cut -f 2 | xargs`
echo "Missing Perl modules:     ${MISSING_PERL_MODULES}"

## Beware: the _noprompt suffix is optional. It has the advantage to
## avoid for the admin to confirm each installation step, but the
## counterpart is that errors may be overlooked.
make SUDO='' -f makefiles/install_rsat.mk perl_modules_install PERL_MODULES="${MISSING_PERL_MODULES}"

## Note: I had to force installation for the some modules, because
## there seem to be some circular dependencies.
make -f makefiles/install_rsat.mk perl_modules_check
MISSING_PERL_MODULES=`grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;' | grep -v "Object::InsideOut"`
echo "Missing Perl modules:     ${MISSING_PERL_MODULES}"
make SUDO='' -f makefiles/install_rsat.mk perl_modules_install_by_force PERL_MODULES_TO_FORCE="`grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;' | grep -v Object::InsideOut| cut -f 2 | xargs`"

## Last check for Perl modules. 
## If some of them still fail (except Object::InsideOut), manual intervention will be required.
make -f makefiles/install_rsat.mk perl_modules_check
more check_perl_modules_eval.txt

## Measure remaining disk space
df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_perl_modules_installed.txt
grep ${DEVICE} ${RSAT_PARENT_PATH}/install_logs/df_*.txt

