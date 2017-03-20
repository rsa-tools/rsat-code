# Introduction

This document explains how to install and configure the Regulatory Sequence Analysis Tools (RSAT). 

# Downloading RSAT

1. Go to the RSAT portal (**<http:rsat.eu>**)
2. Select any server. 
3. In the left-sided menu,  click on the **Download** link. 
4. Fill in your name and coordinates and accept the license.
You now have access to the download page containing the links to  
     * the latest stable release named  
     `rsat_20XX-XX-XX.tar.gz`  
     (XX-XX-XX are replaced by the release date).
  
     * the previous release archives available in the folder  
       `previous_versions`.
	 
## Downloading the latest RSAT release via the Web installer  

5. Download the tar archive named  
`rsat_20XX-XX-XX.tar.gz`  
where `XX-XX-XX` is the latest release date and put it in your chosen directory.  

6. Uncompress the archive. This will create a directory named `rsat` where we will continue the install procedure.

```
tar -xpvzf rsat_20XX-XX-XX.tar.gz
cd rsat
```

7. In the `rsat` directory, type:

```

```






## Files

rsat_YYYY-MM-DD.tar.gz

   Stand-alone versions + web servers for RSAT and NeAT (except the
   metabolic tools).

metabolic-tools_YYYYMMDD.tar.gz

   Metabolic pathway analysis tools (supported on some NeAT servers).


================================================================

RSAT/NeAT installation
======================

After having uncompressed the archive, you will find the installation
and user guides in the directory

      rsa-tools/doc/manuals/*.pdf

Regulatory Sequence Analysis Tools (RSAT)
=========================================
RSAT installation guide:   RSAT_install_guide.pdf
Web configuration guide:   rsat_web_server.pdf
Command-linde user guide:  tutorial_shell_rsat.pdf

Network Analysis Tools (NeAT)
=============================
Web server configuration:  neat_web_server.pdf
Command-line user guide:   neat_tutorial.pdf
