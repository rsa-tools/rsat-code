

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




