## BEWARE: this protocol is reserved to the members of the RSAT
## development team who want to pass from the CVS to the GIT
## distribution. The transfer should be done only once, in a careful
## step-by-step way, because missing one step may cause problems for
## the subsequent ones. 

## BEFORE DOING THE FOLLOWING, you must configure your new git RSAT
## distribution by following the instructions in
## rsat/doc/manuals/install_guide.pdf

### After all the configuration steps of the manual have been
### performed, you can recuperate the data in this way:

## You must start this from the "rsat" folder of the git distribution
## cd rsat

## Load the new environment variables
source RSAT_config.bashrc

## Check that your RSAT environment variables points to the new distribution folder
echo $RSAT

## Move the supported genomes + other data types from the old to the new location
## and create links in oder for the old distribution to be still functional. 
TO_MOVE=genomes supported_organisms.tab taxon_frequencies taxon_frequencies published_data neat_tuto_data Stored_networks phylogeny
for f in ${TO_MOVE}  ; do \
	echo "Moving $${f} from ../rsa-tools/public_html/data to public_html/data" ; \
	mv ../rsa-tools/public_html/data/$${f} public_html/data/ ; \
	cd ../rsa-tools/public_html/data/;  ln -s ../../../rsat/public_html/data/$${f} .; cd -; \
done 

## Check that the supported organisms are taken into account by your new config  
supported-organisms

## (optional) 
## Transfer the log files
## We synchronize them, rather than moving them, so that from now on the logs can be pursued separately on the old and new instances.
rsync --exclude CVS -ruptvl ../rsa-tools/public_html/logs public_html/

