#!/usr/bin/perl -w
############################################################
#
# $Id: process-pathwayinference-output.pl,v 1.2 2012/03/09 09:29:06 rsat Exp $
#
############################################################

## use strict;

=pod

=head1 Extract Pathway from gene list

=head1 VERSION 1.0

=head1 DESCRIPTION

This code processes the pathway inference ouput subgraphe file:
1. Annotation of the inferred pathway: identify the EC numbers,
enzymes and genes associated to each reaction (seed + inferred
reactions) of the inferred pathway. This documentation relies on the
same GPR file as the gene to reaction mapping.

2. Create an dot file and an image files

This implementation requires no database or Internet connection and
works with local files only.  The PathwayInference tool wraps a number
of algorithms to infer pathways: k shortest paths (REA), kWalks and a
hybrid approach combining both (Faust, et al., 2010). In addition, two
Steiner tree algorithms are available (Takahashi-Matsuyama and
Klein-Ravi), each of them alone or in combination with kWalks.

=head1 AUTHORS

The java PathwInference tool was developed by Karoline Faust. This
Perl wrapper was developed by Didier Croes. The doc was written by
Didier Croes and Jacques van Helden.

=head1 REFERENCES

Faust, K., Croes, D. and van Helden, J. (2011). Prediction of
metabolic pathways from genome-scale metabolic networks. Biosystems.

Faust, K., Dupont, P., Callut, J. and van Helden, J. (2010). Pathway
discovery in metabolic networks by subgraph extraction. Bioinformatics
26, 1211-8.

Faust, K. and van Helden, J. (2011). Predicting metabolic pathways by
sub-network extraction.  Methods in Molecular Biology in press, 15.

=head1 CATEGORY

Graph tool

=head1 USAGE

process_pathwayinference_output -h  [-i inputfile] [-o output_directory] [-v verbosity] \
    -a gene2ec -b ec/rxn/cpd2rxnid/cpdid [-d unique_descriptor] [-t temp_directory] [-show]

=head1 INPUT FORMAT



----------------------------------------------------------------
=head2 Seed mapping file

The seed mapping file makes the link between different types of seeds
(genes, EC numbers, proteins, compound names) and nodes of the network
(reactions or compounds depending on the seed type).

=head3 EC-REACTION-COMPOUND (ECR) file (option a<-ecr>)

Mandatory.

The ECR files makes the link between EC numbers/rxn name/cpd name and reaction id/compound id .

These files are used for the reaction ids/compound ids to gene annotation (backward).


=head3 Example of ECR file

 #query  id      qualifier       name
1.-.-.- RXN1G-1486      EC      3-oxo-C78-α-mycolate-reductase
1.-.-.- RXN1G-1527      EC      3-oxo-C85-cis-methoxy-mycolate reductase
1.-.-.- RXN1G-1528      EC      3-oxo-C86-trans-methoxy-mycolate-reductase
10-deoxysarpagine       10-DEOXYSARPAGINE       compounds       10-deoxysarpagine
10-DEOXYSARPAGINE       10-DEOXYSARPAGINE       compounds       10-deoxysarpagine
 ...

=head3 Gene-EC (GE) file (option b<-ge>)

Not Mandatory. the graph will be annotated with EC only

The GE files makes the link between GENE and EC .

These files are used for the reaction ids/compound ids to gene annotation (backward).


=head3 Example of GE file

#query  id      qualifier       name
bisc	1.1.1.-	GENE	bisc
ahpc	1.11.1.15	GENE	bisc
kdu	1.1.1.125	GENE	kdu

=head1 EXAMPLES

=head2 With an input file

=head3 Motivation

Get methionine-related genes in Escherichia coli genome. This
generates a file containing one line per gene and one column per
attribute (ID, start, end, name, ...).


=head3 Commands

Extract all E.coli genes whose name starts with met
 gene-info -org Escherichia_coli_K_12_substr__MG1655_uid57779 -feattype CDS -full -q '^met.*' -o met_genes.tab

Select the first column, containing gene Ids.
 grep -v "^;" met_genes.tab | cut -f 1 > met_genes_IDs.txt

Extract a pathway connecting at best the reactions catalyzed by these gene products
  pathway-extractor -i met_genes_IDs.txt \
     -g data/networks/MetaCyc/MetaCyc_directed_141.txt \
     -ger ${RSAT}/data/metabolic_networks/GER_files/GPR_Uniprot_112011_Escherichia_coli_K12.tab \
     -o result_dir \
     -t temp_dir

----------------------------------------------------------------

=head2 Using standard input

=head3 The script pathway-extractor can also use as input the STDIN. This allows to use it in aconcatenation of commands. For example, all the commands above could be combined in a single pipeline as follows.

 gene-info -org Escherichia_coli_K_12_substr__MG1655_uid57779 -feattype CDS -q 'met.*' \
   | grep -v "^;" met_genes.tab | cut -f 1 \
   | pathway-extractor -g data/networks/MetaCyc/MetaCyc_directed_141.txt \
       -ger data/networks/MetaCyc/METACYC_GPR_EC_20110620.txt \
       -o result_dir -t temp_dir

----------------------------------------------------------------

=head1 OUTPUT FILES:

*_converted_seeds.txt : pathway inference input file.

*_pred_pathways.txt : result graph file of the pathway inference

*_annot_pred_pathways.txt : : result file of the pathway inference with gene and EC annotation

*_annot_pred_pathways.dot : result file in dot format, which can be converted to a figure using the
automatic layout program dot (included in the software suite graphviz).

*._annot_pred_pathways.png : png image file of the inferred pathway
----------------------------------------------------------------

=cut


BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
	push (@INC,"$ENV{RSAT}/perl-scripts/lib/");
    }
}
require "RSA.lib";


################################################################
## Main package
package main;
{

  ################################################################
  ## Initialise parameters
  #
  local $start_time = &RSAT::util::StartScript();
  $program_version = do { my @r = (q$Revision: 1.2 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r };
  #    $program_version = "0.00";

  ## Input/output files
  our %infile = ();	     # File name containing a list of genes ID
  our %outfile = ();

  ## Directories
  our %dir = ();
  $dir{output} = "."; # output directory
  $dir{temp}= "";     # temporary directory

  our $verbose = "";
  our $in = STDIN;
  our $out = STDOUT;


  $infile{gec} =""; # GPR Gene -> EC -> REACTION annotation file path. Default (METACYC_GPR_EC.tab)
  $infile{ecr}=""; 
  $show = 0;		# Open png image in gwenview
#   $group_descriptor= ""; # Unique name to differenciate output files

  my $organism = "Unknown";
  my $organism_id;
  # my $working_dir = "";
  my $query_ids;
  my @query_id_list;

  ################################################################
  ## Read argument values
  &ReadArguments();

  ################################################################
  #generate ouput files name
  $outfile{prefix} = $dir{output}."/";
  my $outputfile =  $infile{input};
  $outputfile =~ s{.*/}{};# remove path
  $outputfile =~ s{\.[^.]+$}{};# remove file extension
  $outfile{prefix} =~ s|//|/|g; ## Suppress double slashes
   print  $outfile{prefix}."\n";
  $outfile{prefix} .= $outputfile;
  $outfile{graph_png} = $outfile{prefix}."_annot.png";
  $outfile{graph_dot} = $outfile{prefix}."_annot.dot";
  $outfile{graph_annot} = $outfile{prefix}."_annot.txt";
  print  $outfile{graph_annot}."\n";
  
    ################################################################
    # Loading reactions from extracted  graph
    
#     &RSAT::message::TimeWarn("Loading reactions from extracted graph") if ($verbose >= 1);
    open (INFILE, '<'.$infile{input}) or die "couldn't open the file!";
    my $i = 0;
    my $stop = "";
    my $line;
    my $reactioncpdquery="";
    while ($line=<INFILE>) {
      #$line = $_;
      chomp  ($line );
      if (length($line)>0 && !($line=~m/^;/)) {
	my @tempdata = split(/\t/,$line);
	if ($tempdata[6] &&(($tempdata[6] eq "Reaction")|| ($tempdata[6] eq "Compound"))) {
	  #       print "|".$line."|"."\n";
	  $tempdata[0]=~s/<$|>$//;
	  $i++;
	  $reactioncpdquery = $reactioncpdquery."(\$2~\"^$tempdata[0]\")||";
	}
      } elsif ($i>0) {
	last;
      }
    }
    close (INFILE);
    # End of Loading results graph reactions
    ################################################################

    ################################################################
    # Searching all reactions information for the reaction that are in  the inferred pathway graph
    &RSAT::message::TimeWarn("Searching information about extracted reactions") if ($verbose >= 1);
    $reactioncpdquery =~s/\|+$//;
   

    my $command_ = "awk -F'\\t+' ' $reactioncpdquery {print \$2\"\\t\"\$1\"\\t\"\$3\"\\t\"\$4}' $infile{ecr}|sort|uniq";
#    print "$command_\n";
    &RSAT::message::TimeWarn($command_) if ($verbose >= 1);
    my @conversiontable = qx ($command_);
    chomp(@conversiontable);

    ## Storing reaction infos in a hash for faster search
    my %reactioninfos=();
    undef @previousarray;
    my @reacinfoarray=();
    foreach my $content (@conversiontable) {
      my @currentarray = split(/\t/,$content);
      if ( @previousarray && !($previousarray[0] eq $currentarray[0])) {
	# 	  print $previousarray[0]."\n";
	my @truc = @reacinfoarray;
	$reactioninfos{$previousarray[0]}=\@truc;
	undef @reacinfoarray;
      }
      push (@reacinfoarray,\@currentarray);
      @previousarray = @currentarray;
    }
    $reactioninfos{$previousarray[0]}=\@reacinfoarray;
    
    ## if gene file get ECs from conversion table
    my %gecinfos=();
    undef @previousarray;
    my @gecinfoarray=();
    if ($infile{gec}){
      my $ec2genequery="";
      while(my ($rxnid, $infoarray) = each(%reactioninfos)) {
	my @values = @{$infoarray};
	foreach my $info_ref (@values) {
	  my @info = @{$info_ref};  
	  if($info[2] eq "EC"){
	    &RSAT::message::Info("EC2GENE:", $info[0]) if ($verbose >= 3);
	    $ec2genequery = $ec2genequery."(\$2~\"^$info[1]\")||"
	  }
	}
      }
       $ec2genequery =~s/\|+$//;
      
      #buid an hash for ec 2 genes
      
      $command_ = "awk -F'\\t+' ' $ec2genequery {print \$2\"\\t\"\$1\"\\t\"\$3\"\\t\"\$4}' \"$infile{gec}\"|sort|uniq";
      &RSAT::message::TimeWarn($command_) if ($verbose >= 1);
      my @geconversiontable = qx ($command_);
      foreach my $content (@geconversiontable) {
# 	&RSAT::message::Info("EC2GENE:", $content) if ($verbose >= 3);
	my @currentarray = split(/\t/,$content);
	if ( @previousarray && !($previousarray[0] eq $currentarray[0])) {
	# 	  print $previousarray[0]."\n";
	  
	  my @truc = @gecinfoarray;
	  &RSAT::message::Info("ECS2GENE:", $previousarray[0] ."\t" . $truc[0][1]) if ($verbose >= 3);
	  $gecinfos{$previousarray[0]}=\@truc;
	  undef @gecinfoarray;
	}
# 	&RSAT::message::Info("GENE:", $previousarray[0]) if ($verbose >= 3);
	push (@gecinfoarray,\@currentarray);
	@previousarray = @currentarray;
      }
      $gecinfos{$previousarray[0]}=\@gecinfoarray;
      
       &RSAT::message::Info("GENE:", "|".$previousarray[0] ."|\t" . $gecinfos{"5.3.1.16"}[0][1]) if ($verbose >= 3);
    }
#     exit 1;
 

    # End of Searching all reactions information for the reaction that are in  the infered pathway graph
    ################################################################

    ################################################################
    # Adding description to the pathway graph
    &RSAT::message::TimeWarn("Adding descriptions to pathway graph") if ($verbose >= 1);
    open (INFILE, '<'.$infile{input}) or die "couldn't open the file!";
#    my $outfile{graph_annot} = $dir{output}.(join "_",$group_descriptor, $groupid, $graph, "annot_pred_pathways.txt");
    # my $outfilename = `mktemp $outfile{graph_annot}`;
    open (OUTFILE, '>'.$outfile{graph_annot}) or die "couldn't open the file!";
    #     print $group_descriptor;
    while (<INFILE>) {
      my($line) = $_;
      chomp($line);
      my @tempdatab = split(/\t/,$line);
      if (length($line)==0 || $line=~ m/^;/) {
	print OUTFILE $line. "\n";
      } else {

	my $tempstring = $tempdatab[0];
	$tempstring=~s/<$|>$//;
	$tempdatab[0] = $tempstring; #remove directionality from reaction node id to merge the nodes
	# 	      print "TEMPSTRING = $tempstring\n";
	my $values_ref = $reactioninfos{$tempstring};
	if (defined $values_ref) {
	  my @values = @{$values_ref};
	  if ($tempdatab[6] && ($tempdatab[6] eq "Reaction")) {

	    my $label="";
	    my $labelb;
	    my $ecs;
	    my $reactionid;
	    foreach my $info_ref (@values) {
	      my @info = @{$info_ref};

	      #  	      print "JOIN=".join("\t", $myarray[0][0])."\n";
	      my($reacid,$ec,$qualif) = @info;
	      # 		    print "ec: $ec\n";
	      if ($ec) {
		if ($qualif eq "EC"){
		  if (%gecinfos){
		    chomp($ec);
		    my $genes = $gecinfos{$ec};
		    &RSAT::message::Info("EC:", "|".$ec."|",  $gecinfos{"$ec"}[0][1]) if ($verbose >= 4);
		    if (defined $genes) {
		      my @genesarray = @{$genes};
		      foreach my $genearrayref (@genesarray) {
			my @genearray = @{$genearrayref};
			my($id,$genename,$qualif) = @genearray;		  
			chomp($genename);
			$label.= "$genename,";
		      }
		    }
		  }
		  $ecs .= $ec;
		}
		if (!defined $reactionid) {
		  $reactionid = $reacid;
		}
	      }
	    }
	    $labelb = "\\n$ecs\\n($reactionid)";
	    
	    $tempdatab[3] = $label.$labelb;
	  } elsif ($tempdatab[6] &&($tempdatab[6] eq "Compound")) {
	    $tempdatab[3] =  $values[0][1];
	  } else {
	    $tempdatab[1]=~s/<$|>$//;
	  }
	}
	print OUTFILE (join "\t",@tempdatab). "\n";
      }
    }
    close (OUTFILE);
    # End of Adding description to the pathway graph
    ################################################################

    ################################################################
    # Converting graph to dot format
    &RSAT::message::TimeWarn("Converting graph to dot format") if ($verbose >= 1);
#    my $outfile{graph_dot} = $dir{output}.(join "_",$group_descriptor, $groupid, $graph, "annot_pred_pathways.dot");
    my $convert_graph_cmd = "convert-graph -from path_extract -to dot -i $outfile{graph_annot} -o $outfile{graph_dot}";
    system  $convert_graph_cmd;
    # End of Converting graph to dot graph format
    ################################################################

    ################################################################
    # Converting dot graph to image with graphviz dot
    &RSAT::message::TimeWarn("Creating graph image with graphviz dot") if ($verbose >= 1);
#    my $outfile{graph_png} = $dir{output}.(join "_",$group_descriptor, $groupid, $graph, "annot_pred_pathways.png");
    my $graph_image_cmd = "dot -Tpng -Kdot -o $outfile{graph_png} $outfile{graph_dot}";
    system $graph_image_cmd;
    #exit 0;
    if ($show) {
      system "gwenview $outfile{graph_png} &";
    }
 
  # End of Converting dot graph to image with graphviz dot
  ################################################################

  ################################################################
  ## Insert here output printing
  ## Report execution time and close output stream
  &RSAT::message::TimeWarn("Ending") if ($verbose >= 1);
  my $exec_time = &RSAT::util::ReportExecutionTime($start_time); ## This has to be exectuted by all scripts
  print $out $exec_time if ($verbose); ## only report exec time if verbosity is specified
  close $out if ($outfile{output});

  exit(0);
}

################################################################
################### SUBROUTINE DEFINITION ######################
################################################################



################################################################
## Display full help message
sub PrintHelp {
  system "pod2text -c $0";
  exit()
}

################################################################
## Display short help message
sub PrintOptions {
  &PrintHelp();
}

################################################################
## Read arguments
sub ReadArguments {
  my $arg;
  my @arguments = @ARGV; ## create a copy to shift, because we need ARGV to report command line in &Verbose()
  while (scalar(@arguments) >= 1) {
    $arg = shift (@arguments);
    ## Verbosity

=pod

=head1 OPTIONS

=over 4

=item B<-v>

Verbose mode

=cut
    if ($arg eq "-v") {
      if (&IsNatural($arguments[0])) {
	$verbose = shift(@arguments);
      } else {
	$verbose = 1;
      }


=pod

=item B<-h>

Display full help message

=cut
    } elsif ($arg eq "-h") {
      &PrintHelp();


=pod

=item B<-help>

Same as -h

=cut
    } elsif ($arg eq "-help") {
      &PrintOptions();

=pod

=item B<-show>

execute gwenview to display the pathway results in png format

=cut
    } elsif ($arg eq "-show") {
     $show = 1;

=pod

=item B<-i inputfile>

If no input file is specified, the standard input is used.  This
allows to use the command within a pipe.

=cut
    } elsif ($arg eq "-i") {
      $infile{input} = shift(@arguments);

=pod

=item	B<-gec GE Genes file>

Gene -> EC (GE) annotation file.

=cut
    } elsif ($arg eq "-gec") {
      $infile{gec} = shift(@arguments);
=pod

=item	B<-ecr ECR file>

EC -> REACTION and COUMPOUNDS (ECR) annotation file.

=cut
    } elsif ($arg eq "-ecr") {
      $infile{ecr} = shift(@arguments);      

# =item	B<-b GR Gene -> REACTION annotation>
#
# An gene annotation file with diredt link gene to reaction. Does not rely on the EC number annotation
#
#
# =pod
#
# =cut
#     } elsif ($arg eq "-b") {
#       $outfile{gr} = shift(@arguments);

# =pod

# =item	B<-n Graph name>
#
# Name of the Graph (default: Name of the graph file).
#
# =cut
#     } elsif ($arg eq "-n") {
#       $graph = shift(@arguments);
#

=pod

=item	B<-o output Directory>

If no output file is specified, the current directory is used.

=cut
    } elsif ($arg eq "-o") {
      $dir{output} = shift(@arguments);

=pod

=item	B<-t temp Directory>

If no output file is specified, the current directory is used.

=cut
    } elsif ($arg eq "-t") {
      $dir{temp} = shift(@arguments);

    } else {
      &FatalError(join("\t", "Invalid pathway_extractor option", $arg));

    }
  }

=pod

=back

=cut

}
################################################################
## Verbose message
sub Verbose {
    print $out "; template ";
    &PrintArguments($out);
    printf $out "; %-22s\t%s\n", "Program version", $program_version;
    if (%infile) {
	print $out "; Input files\n";
	while (my ($key,$value) = each %infile) {
	  printf $out ";\t%-13s\t%s\n", $key, $value;
	}
    }
    if (%outfile) {
	print $out "; Output files\n";
	while (my ($key,$value) = each %outfile) {
	  printf $out ";\t%-13s\t%s\n", $key, $value;
	}
    }
}


__END__
