
# Regulatory Sequence Analysis Tools Docker container

If you want to fetch a ready to run Docker container with RSAT please choose
one from https://hub.docker.com/r/biocontainers/rsat/tags .

Please follow these steps to install it:

 1. Download a container

        docker pull biocontainers/rsat:20230602_cv1

 2. Create local folders for input data and results, outside the container, as these might be large

        mkdir -p ~/rsat_data/genomes ~/rsat_results
        chmod -R a+w ~/rsat_data/genomes ~/rsat_results

 3. Launch Docker RSAT container:

        docker run --rm -v ~/rsat_data:/packages/rsat/public_html/data/ -v ~/rsat_results:/home/rsat_user/rsat_results -it biocontainers/rsat:20230602_cv1


 4. (From container terminal) Download organism from public RSAT server, such as the Plants server. Other available servers are http://fungi.rsat.eu, http://metazoa.rsat.eu, http://protists.rsat.eu and http://teaching.rsat.eu

        download-organism -v 2 -org Prunus_persica.Prunus_persica_NCBIv2.38 -server https://rsat.eead.csic.es/plants

 5. Test container

        cd rsat_results 
        make -f ../test_data/peak-motifs.mk RNDSAMPLES=2 all

 6. Install any organism, please follow instructions at [managing-RSAT](https://rsa-tools.github.io/managing-RSAT/genome_installation/install_organisms_FASTA_GTF.html)

 7. To connect to RSAT Web server running from Docker. If you run Docker as a normal user Apache will not start properly and you will see these messages:

        * Starting Apache httpd web server apache2
        (13)Permission denied: AH00091: apache2: could not open error log file /var/log/apache2/error.log.
        AH00015: Unable to open logs
        Action 'start' failed.
        The Apache error log may have more information

 If you really want lo launch the Docker Web server launch tge container and do (see RSATPASS above):
 
        sudo service apache2 restart
        hostname -I

 Finally open the following URL in your browser, using the obtained the IP address: http://172.17.0.2/rsat

## Dockerfiles

The Dockerfiles for the containers at [biocontainers/rsat](https://hub.docker.com/r/biocontainers/rsat) 
can be found at https://github.com/BioContainers/containers/tree/master/rsat

Those are based on this [Dockerfile](./Dockerfile), which we use for testing and builds a standalone Docker container of RSAT.
