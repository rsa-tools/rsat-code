#!/usr/bin/perl
############################################################
#
# $Id: parse_genbank.pl,v 1.3 2002/05/10 20:55:20 jvanheld Exp $
#
# Time-stamp: <2002-05-10 22:54:33 jvanheld>
#
############################################################

### parse_genbank.pl
### type parse_genbank.pl -h for info

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`"); ### add the program's directory to the lib path
}
require "PFBP_classes.pl";
require "PFBP_config.pl";
require "PFBP_util.pl";
require "PFBP_loading_util.pl"; ### for converting polypeptide IDs into ACs
require "PFBP_parsing_util.pl";



package main;
{
    
    #### initialization
    $verbose = 0;
    $out_format = "obj";
    
    
    #### files to parse
    @selected_organisms = ();
    push @selected_organisms, qw ( Mycoplasma_genitalium ); 
    push @selected_organisms, qw ( Escherichia_coli_K12 );

    $dir{output} = $parsed_data."/genbank_parsed/".$delivery_date;
    $dir{genbank} = "/lin/genomics/genbank/ftp.ncbi.nlm.nih.gov/genbank/genomes";


    #### classes and classholders
    @classes = qw( PFBP::Gene );
    $genes = PFBP::ClassFactory->new_class(object_type=>"PFBP::Gene",
					   prefix=>"gene_");

    $contigs = PFBP::ClassFactory->new_class(object_type=>"PFBP::Contig",
					  prefix=>"contig_");

    
    
    &ReadArguments;

    ### default output fields for each class
    if ($rsa) {
	warn "; Exporting in special format for rsa-tools.\n" if ($verbose >=1);
	#### specific export format for RSA-tools
	$single_name = 1;
	$genes->set_out_fields(qw( id type name chromosome start end strand description position names ));
    } else {
	$genes->set_out_fields(qw( id source organism type position chromosome strand start end description names exons introns db_xref ));
	#@{$out_fields{'PFBP::Gene'}} = qw( id source organism raw_position chromosome strand start_base end_base description names exons );
    }
    
    ### outfile names
    unless (defined($suffix)) {
	$suffix .= "_test" if ($test);
    }
    
    unless (-d $dir{output}) {
	warn "Creating output dir $dir{output}\n";
	mkdir $dir{output}, 0775 
	    || die "Error: cannot create directory $dir\n";
    }
    chdir $dir{output};
    $out_file{genes} = "$dir{output}/gene".$suffix.".obj" if ($export{obj});
    $out_file{error} = "$dir{output}/gene".$suffix.".errors.txt";
    $out_file{stats} = "$dir{output}/gene".$suffix.".stats.txt";
    
    ### open error report file
    open ERR, ">$out_file{error}" || die "Error: cannot write error file $out_file{error}\n";
    
    
    ### select all organisms if none was selected (-org)
    if (($export{all}) || 
	($#selected_organisms < 0)) {
	@selected_organisms = &SelectAllOrganisms();
    }
    
    foreach my $org (@selected_organisms) {
	my $org_dir = "Bacteria/${org}/";
	my $data_dir = $dir{genbank}."/${org_dir}";
	die ("Error: directory ", $data_dir, " does not exist.\n") 
	    unless (-d $data_dir);

	@genbank_files = glob($data_dir."/*.gbk");
	push @genbank_files, glob($data_dir."/*.gbk.gz"); #### compressed files are supported
	
#die join "\n", $#genbank_files, @genbank_files;
	if ($#genbank_files <= -1) {
	    die ("Error: there is not genbank file in the directory ", 
		 $data_dir,
		 "\n");
	} elsif ($#genbank_files > 0) {
	    warn ("Error: there are several genbank files in the directory ", 
		 $data_dir, "\n\t",
		 join ("\n\t", @genbank_files),
		 "\n");
	    next;
	} else {
	    $in_file{$org} = $genbank_files[0];
	    push @parsed_organisms, $org;
	}

	$short_file{$org} = `basename $in_file{$org} .gbk`;
	chomp $short_file{$org};
	
	if (-e $in_file{$org}) {
	    if ($in_file{$org} =~ /\.gz$/) {
		$in_stream{$org} = "gunzip -c $in_file{$org} | ";
	    } else {
		$in_stream{$org} = "cat $in_file{$org} | ";
	    }
	} elsif (-e "$in_file{$org}.gz") {
	    $in_stream{$org} = "gunzip -c $in_file{$org}.gz | ";
	} else {
	    die ("Error: cannot find data file for organism ", $org, "\n",
		 "\t", $in_file{$org}, "\n");
	}
       
	### test conditions
	if ($test) {
	    $in_stream{$org} .= " head -1000 |";
	}
    }
    
    if ($verbose >=1) {
	warn ";TEST\n" if ($test); 
	
	&DefaultVerbose;

	warn "; Selected organisms\n;\t", join("\n;\t", @selected_organisms), "\n";
	warn "; Parsed organisms\n";
	foreach my $org (@parsed_organisms) {
	    warn ";\t$org\t",$short_file{$org},"\n";
	}
    }
    
    ### parse data from original files
    foreach $org (@parsed_organisms) {
	&ParseGenbankFile($in_stream{$org}, 
			  $genes, 
			  $contigs, 
			  source=>"genbank:".$short_file{$org},
			  no_seq=>1);
	
#  	### check organism attribute
#  	foreach $gene ($genes->get_objects()) {
#  	    unless ($gene->get_attribute("organism")) {
#  		&ErrorMessage("Warning: gene ", $gene->get_attribute("id"), " has no organism attribute\n");
#  		$gene->set_attribute("organism",$org);
#  		next;
#  	    }
#  	}
    }
    
    &ParsePositions($genes);
    
    foreach $gene ($genes->get_objects()) {
	foreach my $name ($gene->get_attribute("gene")) {
	    $gene->push_attribute("names",$name);
	};
	#$gene->push_attribute("names", $gene->get_attribute("gene"));
	
	### define a single name  (take the first value in the name list)
	if ($single_name) {
	    if ($name = $gene->get_name()) {
		$gene->set_attribute("name",$name);
	    } else {
		$gene->set_attribute("name",$gene->get_id());
	    }
	}
	
	#### check for genes without description
	if (($gene->get_attribute("description") eq $null) 
	    || ($gene->get_attribute("description") eq "")) {
	    $gene->set_attribute("description",$gene->get_attribute("product"));
	}

	#### use GI as gene identifier
	my @xrefs = $gene->get_attribute("db_xref");
	my $gi = "";
	foreach my $xref (@xrefs) {
	    if ($xref =~ /GI:/) {
		$gi = $';
		last;
	    } 
	}

	if ($gi) {
	    $gene->force_attribute("id",$gi);
	} else {
	    &ErrorMessage("; Error\tgene ".$gene->get_attribute("id")." has no GI.\n"); 
	}

	#### use genbank name as chromosome name
	my $source = $gene->get_attribute("source");
	if ($source =~ /genbank:/) {
	    my $chromosome = $';
	    $chromosome =~ s/\.gz$//;
	    $chromosome =~ s/\.gbk$//;
	    $gene->force_attribute("chromosome",$chromosome);
	}



    }
    
    
    ### print result
    chdir $dir{output};
    &PrintStats($out_file{stats}, @classes);
    $genes->dump_tables($suffix);
    $contigs->dump_tables($suffix);
    &ExportClasses($out_file{genes}, $out_format, PFBP::Gene) if $export{obj};
    
    
    ### report execution time
    if ($verbose >= 1) {
	$done_time = &AlphaDate;
	warn ";\n";
	warn "; job started $start_time";
	warn "; job done    $done_time\n";
    }
    
    close ERR;
    
    
    warn "; compressing the files\n" if ($verbose >= 1);
    system "gzip -f $dir{output}/*.tab $dir{output}/*.txt";
    system "gzip -f $dir{output}/*.obj" if ($export{obj});
    
    
    exit(0);
    
}

### subroutines for the main package

### print the help message 
### when the program is called with the -h or -help option
sub PrintHelp {
  open HELP, "| more";
  print <<EndHelp;
NAME
	parse_genbank.pl

DESCRIPTION
	Parse genes from a Genbank file 
	
AUTHOR
	Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)  

VERSION
	0.01
	Created		1999/12/16
	Last modified	2000/01/08
	
OPTIONS	
	-h	detailed help
	-help	short list of options
	-test	fast parsing of partial data, for debugging
	-outdir output directory
	-v #	warn level
		Warn level 1 corresponds to a restricted verbose
		Warn level 2 reports all polypeptide instantiations
		Warn level 3 reports failing get_attribute()
	-obj    Export results in "obj" format (human readable)
	-all	export all polypeptides for the selected organisms
		(default)
	-org	select an organism for exportation
		can be used reiteratively in the command line 
		to select several organisms
		by default, all organisms found in the input directory
		are selected
	-suffix suffix
		add a suffix to output file names
		    -suffix '' prevents from adding any suffix
	-name
		exports a name as single value attribute in
		the main table (this is redundant but can be useful)

EndHelp
  close HELP;
}

  



### read arguments from the command line
sub ReadArguments {
  for my $a (0..$#ARGV) {
    
    ### warn level
    if (($ARGV[$a] eq "-v" ) && 
	($ARGV[$a+1] =~ /^\d+$/)){
      $main::verbose = $ARGV[$a+1];
      
      ### test run
    } elsif ($ARGV[$a] eq "-test") {
      $main::test = 1;
      
      ### export single name in main table
    } elsif ($ARGV[$a] eq "-name") {
      $main::single_name = 1;
      
      ### specific export format  for RSA-tools
    } elsif ($ARGV[$a] eq "-rsa") {
      $main::rsa = 1;
      
      ### suffix
    } elsif ($ARGV[$a] eq "-suffix") {
      $a++;
      $main::suffix = $ARGV[$a];
      
      ### output dir
    } elsif ($ARGV[$a] eq "-outdir") {
      $a++;
      $main::dir{output} = $ARGV[$a];
      
      ### help
    } elsif (($ARGV[$a] eq "-h") ||
	     ($ARGV[$a] eq "-help")) {
      &PrintHelp;
      exit(0);

      ### select enzymes for exportation
    } elsif ($ARGV[$a] =~ /^-org/) {
      push @selected_organisms, $ARGV[$a+1];

      #### export all organisms
    } elsif ($ARGV[$a] =~ /^-all/) {
      $main::export{all} = 1;

      #### export object file
    } elsif ($ARGV[$a] =~ /^-obj/) {
      $main::export{obj} = 1;
    }
    
  }
}





#  sub ParsePositions {
#  ### clean up positions
#      warn ("; ",
#  	  &AlphaDate(),
#  	  "\tparsing gene positions\n")
#  	if ($verbose >= 1);

#      foreach my $gene ($genes->get_objects()) {
#  	my $position = $gene->get_attribute("position");

#  	if ($position eq $null) {
#  	    $gene->set_attribute("position",$null);
#  	    $gene->set_attribute("chromosome", $null);
#  	    $gene->set_attribute("strand",$null);
#  	    $gene->set_attribute("start",$null);
#  	    $gene->set_attribute("end",$null);
#  	    &ErrorMessage("Warning: gene ", $gene->get_attribute("id"), " has no position attribute\n");
#  	    next;
#  	}
#  	my $coord = $null;
#  	my $chomosome = $null;
#  	my $chrom_pos = $null;
#  	my $strand = $null;
#  	my $start = $null;
#  	my $end = $null;

#  	if ($position =~ /^(\S+)\:(.*)/) {
#  	    $chromosome = $1;
#  	    $chrom_pos = $2;
#  	} elsif ($position =~ /^([^\:]*)$/) {
#  	    $chromosome = "genome";
#  	    $chrom_pos = $1;
#  	} else {
#  	    &ErrorMessage("Warning: invalid position",
#  			  "\t", $gene->get_attribute("id"),
#  			  "\t", $gene->get_attribute("organism"),
#  			  "\t", $position, "\n");
#  	    $gene->set_attribute("chromosome", $null);
#  	    $gene->set_attribute("strand",$null);
#  	    $gene->set_attribute("start",$null);
#  	    $gene->set_attribute("end",$null);
#  	    next;
#  	}

#  	if ($chrom_pos =~ /complement\((.*)\)/) {
#  	    $strand = "R";
#  	    $coord = $1;
#  	} else {
#  	    $strand = "D";
#  	    $coord = $chrom_pos;
#  	}

#  	if ($coord =~ /^(\d+)\.\.(\d+)$/) {
#  	    $start = $1;
#  	    $end = $2;
#  	} elsif ($coord =~ /^join\((.*)\)$/) { ### exons
#  	    my @exons = split ",", $1;
#  	    my @exon_starts = ();
#  	    my @exon_ends = ();

#  	    #### exon limits
#  	    foreach my $exon (@exons) {
#  		$exon =~ s/\s+//g;
#  		$gene->push_attribute("exons", $exon);
#  		if ($exon =~ /(\d+)\.\.(\d+)/) {
#  		    push @exon_starts, $1;
#  		    push @exon_ends, $2;
#  		} else {
#  		    &ErrorMessage("Error gene\t",$gene->get_id(),"\tinvalid exon\t$exon\n");
#  		}
#  	    }
	    
#  	    #### gene start and end 
#  	    $start = $exon_starts[0] || $null;
#  	    $end = $exon_ends[$#exons] || $null;

#  	    #### introns
#  	    my @introns = ();
#  	    for my $e (0..$#exon_starts - 1) {
#  		my $intron = $exon_ends[$e] + 1;
#  		$intron .= "..";
#  		$intron .= $exon_starts[$e+1] -1;
#  		$gene->push_attribute("introns", $intron);
#  	    }
#  	} else {
#  	    &ErrorMessage("Warning : gene ",$gene->get_attribute("id"),"\tinvalid gene position $position\n");
#  	    $strand = $null;
#  	    $start = $null;
#  	    $end = $null;
#  	}
#  	$gene->set_attribute("chromosome", $chromosome);
#  	$gene->set_attribute("strand",$strand);
#  	$gene->set_attribute("start",$start);
#  	$gene->set_attribute("end",$end);
#      }
#  }



sub SelectAllOrganisms {
    my @organisms = ();
    my $dir = $dir{genbank}."/Bacteria/";
    die "Error: directory $dir does not exist.\n"
	unless (-d $dir);
    chdir $dir;
    @organisms = glob "*";
    die "Error: there are no organisms in the directory $dir.\n"
	unless $#organisms >=0;
    return @organisms;
}

