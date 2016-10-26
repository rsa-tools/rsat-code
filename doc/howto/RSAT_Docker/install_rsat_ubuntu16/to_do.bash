

# ################################################################
# ## TO BE CHECKED: TO WE STILL NEED TO DO ALL THE TRICKY STUFF BELOW ?
# ## NOT SURE: on 2014/08/15 Jv installed a VM on IFB cloud, and the
# ## SOAP/WSDL Web services seem to work without this.
#
# ## The installation of SOAP:WSDL under cpan is particularly tricky. 
# ## In Ubuntu, there is a way to install it with ${INSTALLER}. 
# ## http://www.installion.co.uk/ubuntu/trusty/universe/l/libsoap-wsdl-perl/fr/install.html
# emacs -nw /etc/apt/sources.list
#
# ## Ensure that the following line is set to "universe"
# deb http://us.archive.ubuntu.com/ubuntu trusty main universe
# ## You can now quit emacs

apt-get update
apt-get --quiet --assume-yes install libmodule-build-perl
apt-get --quiet --assume-yes install libsoap-wsdl-perl

################################################################



################################################################
## Next steps require to be done as rsat administrator user

## compile RSAT programs written in C
cd ${RSAT}
make -f makefiles/init_rsat.mk compile_all
df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_app_compiled.txt

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!  BUG    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!! I HAVE A PROBLEM TO COMPILE KWALKS. SHOULD BE CHECKED !!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


## Install two model organisms, required for some of the Web tools.
download-organism -v 1 -org Saccharomyces_cerevisiae \
 -org Escherichia_coli_K_12_substr__MG1655_uid57779

## Optionally, install some pluricellular model organisms
# download-organism -v 1 -org Drosophila_melanogaster
# download-organism -v 1 -org Caenorhabditis_elegans
# download-organism -v 1 -org Arabidopsis_thaliana

## Get the list of organisms supported on your computer.
supported-organisms

df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_organism_installed.txt

################################################################
## Install selected R librairies, required for some RSAT scripts
################################################################

## Installation of R packages
cd $RSAT

make -f makefiles/install_rsat.mk install_r_packages


# ## The command R CMD INSTALL apparently does not work at this stage.
# ##	root@rsat-tagc:/workspace/rsat# R CMD INSTALL reshape
# ##	Warning: invalid package 'reshape'
# ##	Error: ERROR: no packages specified
# cd $RSAT; make -f makefiles/install_rsat.mk  r_modules_list 
# ### I install them from the R interface. This should be revised to
# ### make it from the bash, but I need to see how to specify the CRAN
# ### server from the command line (for the time being, I run R and the
# ### programm asks me to specify my preferred CRAN repository the first
# ### time I install packages).
# R
# ## At the R prompt, type the following R commands.
# ##
# ## Beware, the first installation of bioconductor may take a while,
# ## because there are many packages to install
# ##
# ## Note: since this is the first time you install R packages on this
# ## VM, you need to choose a RCRAN server nearby to your site.
# install.packages(c("reshape", "RJSONIO", "plyr", "dendroextras", "dendextend"))
# source('http://bioconductor.org/biocLite.R'); biocLite("ctc")
# quit()
# ## At prompt "Save workspace image? [y/n/c]:", answer "n"

## Check remaining disk space
df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_R_packages_installed.txt

################################################################
## Install some third-party programs required by some RSAT scripts.
cd ${RSAT}
make -f makefiles/install_software.mk
make -f makefiles/install_software.mk install_ext_apps
df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_extapp_installed.txt



################################################################
## tests on the Web site

## Run the demo of the following tools
##
## - retrieve-seq to check the access to local genomes (at least
##   Saccharomyces cerevisiae)
##
## - feature-map to check the GD library
##
## - retrieve-ensembl-seq to check the interface to Ensembl
##
## - fetch-sequences to ceck the interface to UCSC
##
## - some NeAT tools (they rely on web services)
##
## - peak-motifs because it mobilises half of the RSAT tools -> a good
##   control for the overall installation.
##
## - footprint-discovery to check the tools depending on homology
##   tables (blast tables).




################################################################
########### Remove unnecessary files to clean space  ###########
###########   BEFORE DELIVERY for VirtualBox         ###########
################################################################

cd ${RSAT_PARENT_PATH}; rm -f  ${RSAT_ARCHIVE} ## Free space
cd ${RSAT}; make -f makefiles/server.mk clean_tmp ## Clean temporary directory
rm -rf ${RSAT}/app_sources
rm -rf ${RSAT}/public_html/tmp/serialized_genomes
find   ${RSAT_HOME}/public_html/tmp/  -maxdepth 1 -mindepth 1 -type d -exec rm -rf {} \;
rm -rf public_html/tmp/www-data ## Clean Apache user temporary directory
rm -rf public_html/tmp/serialized_genomes ## Clean serialized organisms
rm -rf  ~/.ssh/ ## Make sure that the installer ssh config is not included in the distribution


## THE INSTALLATION OF THE RSAT SERVER IS HOW DONE. THE REST IS OPTIONNAL
################################################################

