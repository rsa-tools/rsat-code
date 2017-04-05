source installer/00_config.bash

cd ${RSAT}; source RSAT_config.bashrc ## Reload the (updated) RSAT environment variables

################################################################
########### Remove unnecessary files to clean space  ###########
###########   BEFORE DELIVERY for VirtualBox         ###########
################################################################

cd ${RSAT}; make -f makefiles/server.mk clean_tmp ## Clean temporary directory
find   ${RSAT}/public_html/tmp/  -maxdepth 1 -mindepth 1 -type d -exec rm -rf {} \;
rm -rf ${RSAT}/public_html/tmp/www-data ## Clean Apache user temporary directory
rm -rf ${RSAT}/public_html/tmp/serialized_genomes  ## Clean serialized organisms
rm -rf  ~/.ssh/ ## Make sure that the installer ssh config is not included in the distribution

## Remove the archives and source code of external program
cd ${RSAT_PARENT_PATH}; rm -f  ${RSAT_ARCHIVE} ## Free space
rm -rf ${RSAT}/app_sources ## Remove source code of external programs

df -m . > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_cleaned_unnecessary_files.txt

## THE INSTALLATION OF THE RSAT SERVER IS HOW DONE. THE REST IS OPTIONNAL
################################################################
