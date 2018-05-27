source installer/00_config.bash

cd ${RSAT}; source RSAT_config.bashrc ## Reload the (updated) RSAT environment variables

################################################################
########### Remove unnecessary files to clean space  ###########
###########   BEFORE DELIVERY for VirtualBox         ###########
################################################################

df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_before_cleaned_unnecessary_files.txt

echo "Running clean_tmp task"
cd ${RSAT}; make -f makefiles/server.mk clean_tmp ## Clean temporary directory

echo "Removing www-data in tmp dir"
rm -rf ${RSAT}/public_html/tmp/www-data ## Clean Apache user temporary directory

echo "Removing serialized genomes from tmp dir"
rm -rf ${RSAT}/public_html/tmp/serialized_genomes  ## Clean serialized organisms

echo "Removing all other sub-directories from tmp dir"
find   ${RSAT}/public_html/tmp/ -maxdepth 1 -mindepth 1 -type d -exec rm -rf {} \;

## Remove user-specific temoprary directory
echo "Removing ~/.rsat_tmp_dir"
rm -rf ~/.rsat_tmp_dir/

## Remove the archives and source code of external program
# cd ${RSAT_PARENT_PATH}; rm -f  ${RSAT_ARCHIVE} ## Free space
echo "Removing ${RSAT}/app_sources"
rm -rf ${RSAT}/app_sources ## Remove source code of external programs

## Removing all the install_tests directory
echo "Removing ${RSAT}/install_tests directory"
rm -rf ${RSAT}/install_tests

df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_cleaned_unnecessary_files.txt

## THE INSTALLATION OF THE RSAT SERVER IS HOW DONE. THE REST IS OPTIONNAL
################################################################
