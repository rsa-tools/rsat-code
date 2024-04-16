#!/usr/bin/env bash

source $(dirname $0)/00_config.bash

echo
echo "================================================================"
echo "====    Installing some organisms on this RSAT instance     ===="
echo "================================================================"
echo

cd ${RSAT}; source RSAT_config.bashrc ## Reload the (updated) RSAT environment variables
## Install two model organisms, required for some of the Web tools.
download-organism -v 1 -org Saccharomyces_cerevisiae \
    -org Escherichia_coli_K_12_substr__MG1655_uid57779 \
    -org Escherichia_coli_GCF_000005845.2_ASM584v2

## Optionally, install some pluricellular model organisms
# download-organism -v 1 -org Drosophila_melanogaster
# download-organism -v 1 -org Caenorhabditis_elegans
# download-organism -v 1 -org Arabidopsis_thaliana

## Get the list of organisms supported on your computer.
supported-organisms

df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_organism_installed.txt
