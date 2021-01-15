# Installing genomes from Ensembl Genomes

**Author: Bruno Contreras Moreira**

## Introduction

This document explains how to install genomes and annotations from Ensembl Genomes using mostly the [FTP](ftp://ftp.ensemblgenomes.org/pub) site.
Gene Ontology terms are optional and are obtained from [BioMart](http://plants.ensembl.org/biomart/martview) instead.

Note that while Ensembl covers Vertebrates, Ensembl Non Vertebrates (NV, Ensembl Genomes) includes the other divisions (Protists, Fungi, Plants, Bacteria, Metazoa). These instructions have mostly been tested with **Ensembl Plants**.

This installation procedure can actually be used to **install genomes from other sources** as well, see below.

This protocol uses a makefile found in the rsat directory [makefiles/ensemblgenomes_FTP_client.mk](makefiles/ensemblgenomes_FTP_client.mk) for these tasks.

This document **does not** use the following scripts, which are documented elsewhere: 

- `install-ensembl-genome`, 
- `download-ensembl-genome`,
- `download-ensembl-features`, 
- `download-ensembl-variations`.


## Downloading a genome from another RSAT web server

Usually the fastest way of getting a genome installed is to fetch it from other RSAT server. This way oligo frequencies are not computed; instead thay are copied over:

```{r, engine='bash', eval=FALSE}
download-organism \\
   -server http://rsat.eead.csic.es/plants \\
   -org Arabidopsis_thaliana.TAIR10.29
```

## Installing genome sequences and annotations from Ensembl Genomes

### Check FTP site URL

The current FTP site is at ftp://ftp.ensemblgenomes.org ; should it change in the future, or its folder structure, [../../makefiles/ensemblgenomes_FTP_client.mk](makefiles/ensemblgenomes_FTP_client.mk) will have to be updated.


### Find out current Ensembl Genomes release 

Visit the Web site of your division of interest, such as http://plants.ensembl.org, and check the release statement at the bottom.
Alternatively you can use the REST [info/eg_version endpoint](http://rest.ensembl.org/documentation/info/eg_version)

Ensembl Genome releases are integers and have an offset of 53 with Ensembl releases. 
For instance, Ensembl release 99 corresponds to Ensembl Genomes 46.

```{r, engine='bash', eval=FALSE}
export EGRELEASE=46
```

### Get supported species 

We shall download the current list of supported genomes in the respective division (GROUP).
This can be done by running the following command:

```{r, engine='bash', eval=FALSE}
cd $RSAT
export EGDIVISION=Plants

make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=$EGDIVISION RELEASE=$EGRELEASE organisms

# in the example, this will create $RSAT/data/ensemblgenomes/plants/release-46/species_EnsemblPlants.txt
```

Once this is done then you can install genomes from that release.

## Install all species 

You can install all genomes with these commands:

```{r, engine='bash', eval=FALSE}
cd $RSAT

# download FASTA and GTF files
nohup make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=$EGDIVISION RELEASE=$EGRELEASE download_all_species

# parse input files, extract genomic features and compute oligo frequencies
nohup make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=$EGDIVISION RELEASE=$EGRELEASE install_all_species

# check upstream sequences can be retrieved
nohup make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=$EGDIVISION RELEASE=$EGRELEASE check_all_species
```

The newly installed species will be added to $RSAT/data/supported_organisms.tab and should be listed with the following 
command-line:

```{r, engine='bash', eval=FALSE}
supported-organisms
```


**Note:** this will take a very long time, weeks for Ensembl Plant releases, so it might not be a good idea.

## Install selected species

In case you want to install a single genome you can do that by:

```{r, engine='bash', eval=FALSE}
cd $RSAT
make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=$EGDIVISION RELEASE=$EGRELEASE SPECIES=oryza_longistaminata download_one_species
make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=$EGDIVISION RELEASE=$EGRELEASE SPECIES=oryza_longistaminata install_one_species
```

You can also install selected genomes from older releases:

```{r, engine='bash', eval=FALSE}
cd $RSAT
make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=$EGDIVISION RELEASE=42 organisms
make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=$EGDIVISION RELEASE=42 SPECIES=oryza_sativa download_one_species
make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=$EGDIVISION RELEASE=42 SPECIES=oryza_sativa install_one_species
```

As indicated earlier, the newly installed species are added to $RSAT/data/supported_organisms.tab and should appear in the list produced by:

```{r, engine='bash', eval=FALSE}
supported-organisms
```

## Compute genome stats report

Run this to generate a report of descriptive stats of genomes currently in your system, such as http://rsat.eead.csic.es/plants/data/stats :
```{r, engine='bash', eval=FALSE}
make -f makefiles/ensemblgenomes_FTP_client.mk calc_stats
```

### Install variation data

If variation data is available for your species of interest you can download it with:

```{r, engine='bash', eval=FALSE}
cd $RSAT
make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=$EGDIVISION RELEASE=$EGRELEASE SPECIES=oryza_sativa variations_one_species
```

**Note:** this will update file $RSAT/data/supported_organisms.tab


### Download and Install Compara homologies

Script `get-orthologs-compara` can be used to retrieve homologues (orthologues by default) precomputed at Ensembl Compara.
In order to use it you must first install Compara in your system, which you can do with:

```{r, engine='bash', eval=FALSE}
cd $RSAT
export ENSRELEASE=99

make -f makefiles/ensemblgenomes_FTP_client.mk RELEASE=$EGRELEASE ENSEMBL_RELEASE=$ENSRELEASE GROUP=$EGDIVISION download_compara
make -f makefiles/ensemblgenomes_FTP_client.mk RELEASE=$EGRELEASE ENSEMBL_RELEASE=$ENSRELEASE GROUP=$EGDIVISION parse_compara_match
make -f makefiles/ensemblgenomes_FTP_client.mk RELEASE=$EGRELEASE ENSEMBL_RELEASE=$ENSRELEASE GROUP=$EGDIVISION install_compara

```

All going well you can check the species with installed homologies with:

```{r, engine='bash', eval=FALSE}
get-orthologs-compara -supported_organisms
```

### Install Gene Ontology terms from BioMart

First we shall get the current registry for the division of interest. Only Plants, Metazoa and Fungi are supported.
For instance, for Plants the registry is at http://plants.ensembl.org/biomart/martservice?type=registry

The XML content of the registry must be copied and pasted into file *$RSAT/ext_lib/biomart-perl/conf/martURLLocation.xml* .
For instance, for Plants release 46 it will look like:
```
<?xml version="1.0" encoding="UTF-8"?>

<MartRegistry>
  <MartURLLocation database="plants_mart_46" default="" displayName="Ensembl Plants Genes 46" host="plants.ensembl.org" includeDatasets="" martUser="" name="plants_mart" path="/biomart/martservice" port="80" serverVirtualSchema="plants_mart" visible="1" />
  <MartURLLocation database="plants_snp_mart_46" default="" displayName="Ensembl Plants Variations 46" host="plants.ensembl.org" includeDatasets="" martUser="" name="plants_variations" path="/biomart/martservice" port="80" serverVirtualSchema="plants_mart" visible="1" />
  <MartURLLocation database="plants_sequence_mart_46" default="" displayName="Ensembl Plants Sequences 46" host="plants.ensembl.org" includeDatasets="" martUser="" name="plants_sequences" path="/biomart/martservice" port="80" serverVirtualSchema="plants_mart" visible="0" />
  <MartURLLocation database="plants_genomic_features_mart_46" default="" displayName="Ensembl Plants Genomic Features 46" host="plants.ensembl.org" includeDatasets="" martUser="" name="plants_genomic_features" path="/biomart/martservice" port="80" serverVirtualSchema="plants_mart" visible="0" />
  <MartURLLocation database="ontology_mart_99" default="" displayName="Ontology Mart 99" host="plants.ensembl.org" includeDatasets="" martUser="" name="ontology" path="/biomart/martservice" port="80" serverVirtualSchema="plants_mart" visible="0" />
</MartRegistry>
```

Once this is done we shall update the registry cache with:
```{r, engine='bash', eval=FALSE}
download-ensembl-go-annotations-biomart -reg
```
This will take a few minutes, but then we can do:
```{r, engine='bash', eval=FALSE}
make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=$EGDIVISION RELEASE=$EGRELEASE download_go
```

And for each species for which we want GO terms we can now do:
```{r, engine='bash', eval=FALSE}
make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=$EGDIVISION RELEASE=$EGRELEASE SPECIES=oryza_sativa download_go_annotations
```

### Installing from other sources

The procedures described can also be used to install arbitrary genomes from other sources,
provided that 4 input files are available with the following extensions: 

* SPECIES_RSAT_ID.dna.toplevel.fa : raw genomic sequence
* SPECIES_RSAT_ID.dna_rm.genome.fa : repeat-hard-masked genomic sequence 
* SPECIES_RSAT_ID.gtf : annotation file 
* SPECIES_RSAT_ID.pep.all.fa : peptide sequences of CDS features

where SPECIES_RSAT_ID is a string identifying this organism and its annotation. 
For instance, for assembly Wm82.a2.v1 of *Glycine max* from [JGI Phytozome](https://phytozome.jgi.doe.gov), we could install with:
 
```{r, engine='bash', eval=FALSE}
cd $RSAT

mkdir -p $RSAT/data/genomes/Glycine_max.Wm82.a2.v1.JGI/genome
# put there those 4 files (dna.toplevel.fa,dna_rm.genome.fa,.gtf,.pep.all.fa)

make -f makefiles/ensemblgenomes_FTP_client.mk SPECIES=Glycine_max \
    SPECIES_DIR=/var/www/html/rsat/data/genomes/Glycine_max.Wm82.a2.v1.JGI \
    SPECIES_RSAT_ID=Glycine_max.Wm82.a2.v1.JGI TAXON_ID=3847 GTF_SOURCE=JGI \
    install_from_gtf
```

Note that TAXON_ID can be obtained at <https://www.ncbi.nlm.nih.gov/taxonomy>

*********
