#!/usr/bin/perl -w
############################################################
#
# $Id: PathwayExtraction.pm,v 1.14 2012/10/04 07:57:26 rsat Exp $
#
############################################################

## use strict;

=pod

=head1 Extract Pathway from gene list

=head1 VERSION 1.0

=head1 DESCRIPTION

This tool is a Perl Wrapper around the stand-alone application
(PathwayInference) developed by Karoline Faust.

The PathwayInference tool can also be used via the Web interface
"Pathway Extraction" (http:// rsat.ulb.ac.be/neat).

The Perl wrapper performs several steps to enable the extraction of
pathways from sets of functionally related genes (e.g. co-expressed,
co-regulated, members of the same operon, …).

1. Gene to reaction mapping. Identify the set of reactions ("seed
reactions") catalyzed by the set of input genes ("seed genes"). The
mapping relies on a user-specified file describing the mapping of
genes to reactions (GNN and NNN, Gene-Node Name and Network Node Name<>NodeID file).

2. Pathway extraction (=Pathway inference). PathwayInference takes as
input a network (typically a metabolic network made of compounds +
reactions) and a set of "seed" nodes. The program attempts to return a
subnetwork that interconnects the seed nodes at a minimal “cost",pathway-extractor_seeds.php
where the cost is a weighted sum of the intermediate compounds and
reactions used to link the seed nodes (Faust, et al., 2011; Faust and
van Helden, 2011).

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

pathway_extractor -h -hp [-i inputfile] [-o output_directory] [-v verbosity] \
    -g infile{graph} -a gene2ec -b ec/rxn/cpd2rxnid/cpdid [-d unique_descriptor] [-t temp_directory] [-show]

=head1 INPUT FORMAT

Warning: the same gene identifiers should be used in all input files.

=head2 1) List of seed genes (gene identifiers):

(Warning at least 2 gene ids must be present in the graph file see
below) in this example we use gene IDs. Beware, the gene IDs must be
compatible with the genome version installed on RSAT. Depending on the
organism source, the IDs can correspond to RefSeq, EnsEMBL, or other
references.

Example of seed gene file:
NP_414739
NP_414740
NP_414741
NP_416617
NP_417417
NP_417481
NP_418232
NP_418272
NP_418273
NP_418373
NP_418374
NP_418375
NP_418376
NP_418437
NP_418443

----------------------------------------------------------------

=head2 2)Graph file format:

see Pathwayinference tools helppathway_extractor

java graphtools.algorithms.Pathwayinference –h

The same result can be obtained by typing

pathway-extractor -hp

----------------------------------------------------------------



=head2 Seed mapping file

The seed mapping file makes the link between different types of seeds
(genes, EC numbers, proteins, compound names) and nodes of the network
(reactions or compounds depending on the seed type).

=head3 Network Node Names (nnn) file (option I<-nnn>)

Mandatory.

The NNN files makes the link between EC numbers/rxn name/cpd name and node id in the network .

These files are used for the reaction ids/compound ids to gene annotation (backward).


=head3 Example of NNN file

 #query  id      qualifier       name
1.-.-.- RXN1G-1486      EC      3-oxo-C78-α-mycolate-reductase
1.-.-.- RXN1G-1527      EC      3-oxo-C85-cis-methoxy-mycolate reductase
1.-.-.- RXN1G-1528      EC      3-oxo-C86-trans-methoxy-mycolate-reductase
10-deoxysarpagine       10-DEOXYSARPAGINE       compounds       10-deoxysarpagine
10-DEOXYSARPAGINE       10-DEOXYSARPAGINE       compounds       10-deoxysarpagine
 ...

=head3 Gene to Node Names (gnn) file (option I<-gnn>)

not mandatory, in this cas the queyr file or stdin must contains only Node Names.

The GNN files makes the link between external identifier and node names (example: gene name, refseq, locuslink) . 


=head3 Example of NNN file

#query	id	qualifier	name	taxonomy_id	species
aas	2.3.1.40	GENE_NAME	aas	83333	Escherichia coli (strain K12)
aas	6.2.1.20	GENE_NAME	aas	83333	Escherichia coli (strain K12)
aat	2.3.2.6	GENE_NAME	aat	83333	Escherichia coli (strain K12)


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
     -gnn ${RSAT}/data/metabolic_networks/GER_files/GPR_Uniprot_112011_Escherichia_coli_K12.tab \
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
use RSAT::util;

## Guess RSAT path from module full name
unless ($ENV{RSAT}) {
    $ENV{RSAT} = $0;
    $ENV{RSAT} =~ s|/perl-scripts/.*||;
    $ENV{RSAT} =~ s|/public_html/.*||;
}

################################################################
## pathwayinference package
package RSAT::PathwayExtraction;

# {  
#   #other pathwayinference otptions : specific opions that will be directly pass to the java pathway inference app
#   our @pioptions=('-s','-i','-g','-f','-b','-n','-q','-o','-O','-E','-a','-y','-p','-F','-Z','-m','-w','-t',
# 		  '-l','-T','-W','-P','-C','-e','-d','-u','-x','-k','-U','-B','-r','-D','-M','-I','-N','-G',
# 		  '-H','-X','-A','-K','-S','-R','-j','-J','-Q','-L','-Y','-v','-h','-V'); 
#   our %otherPIoptions=();
#    ## Input/output files
#   our %infile = ();	     # input file names container
#   our %outfile = ();	     # output file names container
# 
#   ## Directories
#   our %dir = ();
#   $dir{output} = "."; # output directory
#   $dir{temp}= "";     # temporary directory
# 
#   our $verbose = "3";
#   our $in = STDIN;
#   our $out = STDOUT;
# 
#   $infile{gnn} =""; # GPR Gene -> EC -> REACTION annotation file path. Default (METACYC_GPR_EC.tab)
#   $infile{nnn}=""; 
#   $infile{graph}="";		# File containing the graph
# 
#   our $graph = "";		# Graph Name
#   our $group_descriptor= ""; # Unique name to differenciate output files
#   
#    ################################################################
#   ## Read argument values
#   &ReadArguments();
# 
#    ################################################################
#   ## Check argument values
#   
#   my $input = $infile{input};
#   my $isInputFile=0;
#   if ($input){
#     $isInputFile=1;
#   }else{
#     my @query_id_list = <$in>;
#     chomp(@query_id_list);
#     $input = join("\t", @query_id_list);
#   }
#   &RSAT::message::Info("--INPUT ", $input) if ($verbose >= 3);
#    &Inferpathway($input,$isInputFile, $dir{output},$infile{gnn},$infile{nnn},$infile{graph},$dir{temp},$group_descriptor,\%otherPIoptions);
# # &Inferpathway();
# }

sub Inferpathway{

    ################################################################
    ## Initialise parameters
    #
    local $start_time = &RSAT::util::StartScript();
    $program_version = do { my @r = (q$Revision: 1.14 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r };
    #    $program_version = "0.00";
    my ($input,
	$isinputfile,
	$outputdir,
	$gnnfile,
	$nnnfile,
	$graphfile,
	$directed,
	$tempdir,
	$localgroup_descriptor,
	$ECgenericitylevel,
	$verbose,
	$piparameters) = @_;

    my %localotherPIparameters = %{$piparameters} if ($piparameters);
#   if the parameters comes from the function use those one
#   $dir{output} = $outputdir if ($outputdir);    
#   $dir{temp} = $tempdir if ($tempdir);
#   $infile{gnn} = $gnnfile if ($gnnfile);
#   $infile{nnn} = $nnnfile if ($nnnfile); 
#   $infile{graph} = $graphfile if ($graphfile);
#  

    if ($isinputfile){
	($main::in) = &RSAT::util::OpenInputFile($input);
	$input="";
	while (<$main::in>) {
	    $input .= $_;
	}
	#   if no group descriptor use the input file name
	if (!$localgroup_descriptor) {
	    $localgroup_descriptor = $input;
	    $localgroup_descriptor =~ s{.*/}{};     # removes path
	    $localgroup_descriptor=~ s{\.[^.]+$}{}; # removes extension
	}
	close $main::in;
	
    }
#close $in;
    if (!$localgroup_descriptor) {
	$localgroup_descriptor ="stdin";
    }
    
    
    $localgroup_descriptor=~s/(\s|\(|\))+/_/g;
    
    
    
    # if no graph name take the graph file name
    if (!$graph) {
	$graph = $graphfile;
	$graph =~ s{.*/}{};	    # removes path
	$graph=~ s{\.[^.]+$}{};   # removes extension
    }
    $graph=~s/(\s|\(|\))+/_/g;
    
    my $organism = "Unknown";
    my $organism_id;
    # my $working_dir = "";

    ################################################################
    ## Print verbose
#   &Verbose() if ($verbose);

    ################################################################
    ## Execute the command

    ## Check the existence of the output directory and create it if
    ## required
    unless ($outputdir =~ /\/$/) {
	$outputdir =$outputdir."/";
    }
    $outputdir =~s |//|/|g; ## Suppress double slashes
    &RSAT::util::CheckOutDir($outputdir);

    ## Check the existence of the temp directory and create it if
    ## required
    $tempdir=$outputdir unless ($tempdir);
    if (!($tempdir=~m/\/$/)) {
	$tempdir = $tempdir."/";
    }
    &RSAT::util::CheckOutDir($tempdir);

    if (!defined $group_id) {
	$groupid =$localgroup_descriptor;
    }
    $groupid=~s/(\s|\(|\))+/_/g;
    
    ################################################################
    ## Define output file names
    $outfile{gr} = "";		# GR Gene -> REACTION annotation
    $outfile{prefix} =$outputdir."/";
    #removed  from name file $localgroup_descriptor,
    $outfile{prefix} .= join("_", $groupid, $graph, "pred_pathways");
    $outfile{prefix} =~ s|//|/|g; ## Suppress double slashes

    $outfile{predicted_pathway} = $outfile{prefix}.".txt";
    $outfile{seeds_converted} = $outfile{prefix}."_seeds_converted.txt";

    ################################################################
    ## ECR Mapping
    my $patrial;
    ## DIDIER: the ECR mapping should be redone with the new program match-names.
    my ($rxnoutput,$cpdoutput) = &RSAT::PathwayExtraction::MapSeeds($input,$gnnfile,$nnnfile,$patrial,$verbose);
    my %mappedseeds = %{$rxnoutput};
    my %compounds = %{$cpdoutput};
    
    my $seednum= 0; 
#processing reaction seeds
    &RSAT::message::TimeWarn("building seed file: ",$outfile{seeds_converted}) if ($verbose >= 1);
    open (MYFILE, '>'.$outfile{seeds_converted});
    print MYFILE "# Batch file created", &RSAT::util::AlphaDate()."\n";
    print MYFILE "# ECR groups parsed from file: "."\n";
    print MYFILE "# EC number grouping: true". "\n";

    while (my ($key, $val) = each(%mappedseeds)){
	my $dashcount = $key =~ tr/-//;
	if ($dashcount <= $ECgenericitylevel){
	    foreach my $reaction (@{$val}) {
		if ($directed){
		    print MYFILE $reaction .">\t".$reaction. "\n";
		    print MYFILE $reaction ."<\t".$reaction. "\n";
		}
		print MYFILE $reaction ."\t".$key. "\n";
		$seednum++;
	    }
	    print MYFILE "$key\t$groupid\n";
	}
    }

# processing compound seeds   
    while (my ($cpd, $val) = each(%compounds)){
#       print MYFILE $cpd ."\t".$cpd. "\n";
	print MYFILE "$cpd\t$cpd\n";
	print STDERR "ERROR:$cpd\t$cpd\n";
	$seednum++;
    }
    close MYFILE;
    
    
#  if (@previousarray && !($tempdata[0] eq $previousarray[0])) {
#       print MYFILE "$previousarray[0]\t$groupid\n";
#       $seednum++;
#     }
#     	print "$tempdata[1] eq $previousarray[1]\n";
#     if ($directed){
#       print MYFILE $tempdata[1] .">\t".$tempdata[1]. "\n";
#       print MYFILE $tempdata[1] ."<\t".$tempdata[1]. "\n";
#     }else{
#       print MYFILE $tempdata[1] ."\t".$tempdata[1]. "\n";
#     }
#     print MYFILE $tempdata[1] ."\t".$tempdata[0]. "\n";
# 
# 
#   if (@conversiontable) {
#     print MYFILE "$previousarray[0]\t$groupid\n";
#     $seednum++;
#   }
    # END OF  Creating reaction file fo pathway inference
    ################################################################
## TO DO WITH DIDIER: CHECK SEED NUMBER AND SEND WARNING IF TOO BIG.
    ##
    ## Define a parameter max{seed_numbers}. By default, program dies
    ## with error if seeds exceed max{seed_number}, but the max can be
    ## increased by the user with the option -max seeds #.

    if ($seednum > 1) {
	################################################################
	# Running PatywayInference
	&RSAT::message::TimeWarn("Running pathway inference with ", $seednum, "seeds") if ($verbose >= 1);
#    $outfile{predicted_pathway} =$outputdir.(join "_",$localgroup_descriptor, $groupid, $graph, "_pred_pathways.txt");
	our $maxinterseedslength = $localotherPIparameters{"-F"} || "5";
	delete($localotherPIparameters{"-F"});
	our $graphfileformat = $localotherPIparameters{"-f"} || "flat";
	delete($localotherPIparameters{"-f"});
	our $weightpolicy= $localotherPIparameters{"-y"} || "con";
	delete($localotherPIparameters{"-h"});
	delete($localotherPIparameters{"-g"});
	delete($localotherPIparameters{"-o"});
	delete($localotherPIparameters{"-p"});
	delete($localotherPIparameters{"-v"});
	#after the program has handled the mandatory parameters, it will add the remaining ones 
	my $piparameters =" ";

	while( my($key, $val) = each(%localotherPIparameters) ) {
	    $piparameters .= " $key $val ";
	}

	my $pathway_infer_cmd = "java ";
	$pathway_infer_cmd .= " -cp ".$ENV{RSAT}."/java/lib/NeAT_javatools.jar";
	$pathway_infer_cmd .= " -Xmx1000M graphtools.algorithms.Pathwayinference";
	$pathway_infer_cmd .= " -i ".$outfile{seeds_converted};
	$pathway_infer_cmd .= " -A ".$ENV{REA_ROOT};
	$pathway_infer_cmd .= " -K ".$ENV{KWALKS_ROOT};
	$pathway_infer_cmd .= " -F $maxinterseedslength -C -f $graphfileformat";
	$pathway_infer_cmd .= " -p $tempdir";
	$pathway_infer_cmd .= " -E $outputdir";
	$pathway_infer_cmd .= " -b";
	$pathway_infer_cmd .= " -d" if ($directed);
	$pathway_infer_cmd .= " -g $graphfile";
	$pathway_infer_cmd .= " -y $weightpolicy ";
	$pathway_infer_cmd .=  $piparameters;
	$pathway_infer_cmd .= " -o $outfile{predicted_pathway}";
#     -v $verbose

	&RSAT::message::TimeWarn("Pathway inference command", $pathway_infer_cmd) if ($verbose >= 2);
	&RSAT::message::Info("Predicted pathway file", $outfile{predicted_pathway}) if ($verbose >= 1);
	#exit 1;
	&RSAT::util::doit($pathway_infer_cmd, $dry, $die_on_error, $verbose, $batch, $job_prefix);

	## TO DO WITH DIDIER: redirect STDERR/STDOUD to log and err files in the output directory

	# END of Running patywayinference
	################################################################
    } else {
	
	print STDERR "NOT ENOUGH SEEDS. Min 2. I stop here!!\n";
	return "NOT ENOUGH SEEDS. Min 2. I stop here!!\n";
    }
    # End of Converting dot graph to image with graphviz dot
    ################################################################

    ################################################################
    ## Insert here output printing
    ## Report execution time and close output stream
    &RSAT::message::TimeWarn("Ending") if ($verbose >= 1); 
# return an error whenusing Soap:
#    my $exec_time = &RSAT::util::ReportExecutionTime($start_time); ## This has to be exectuted by all scripts
#    print $exec_time if ($verbose); ## only report exec time if verbosity is specified
#   close $out if ($outfile{output});

    return $outfile{predicted_pathway};
}


sub ProcessOutputFiles{
    my ($inputfile,
	$outputdir,
	$gnnfile,
	$nnnfile,
	$directed,
	$verbose) = @_;
    
    my $undirected = "-undirected";
    undef $undirected if ($directed); # if directed $undirected = null
    
    my %outfile = ();
    $outfile{prefix} = $outputdir."/";
    my $outputfile =  $inputfile;
    $outputfile =~ s{.*/}{};# remove path
    $outputfile =~ s{\.[^.]+$}{};# remove file extension
    $outfile{prefix} =~ s|//|/|g; ## Suppress double slashes
#   print STDERR $outfile{prefix}."\n";
    $outfile{prefix} .= $outputfile;
    $outfile{graph_png} = $outfile{prefix}."_annot.png";
    $outfile{graph_dot} = $outfile{prefix}."_annot.dot";
    $outfile{graph_annot} = $outfile{prefix}."_annot.txt";
#   print  STDERR $outfile{graph_annot}."\n";
    &RSAT::message::TimeWarn("Graph annotated file: ",$outfile{graph_annot}, "PWD:$ENV{PWD}" ) if ($verbose >= 1);
    ################################################################
    # Loading reactions from extracted  graph
#     &RSAT::message::Info("RSATPATH: ", $ENV{RSAT}) if ($verbose >= 1);
#     &RSAT::message::Info("REA_ROOT: ", $ENV{REA_ROOT}) if ($verbose >= 1);
#     &RSAT::message::Info("KWALKS_ROOT: ", $ENV{KWALKS_ROOT}) if ($verbose >= 1);
#     &RSAT::message::Info("CLASSPATH: ",$ENV{CLASSPATH}) if ($verbose >= 1); 
#     &RSAT::message::TimeWarn("Loading reactions from extracted graph") if ($verbose >= 1);
    open (INFILE, '<'.$inputfile) || &RSAT::error::FatalError("Could not open file", $inputfile);
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
		$reactioncpdquery = $reactioncpdquery."(\$2~\"^$tempdata[0]\$\")||";
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
    &RSAT::message::TimeWarn("Searching information about extracted reactions") if ($verbose >= 1);my $dashcount = $EC =~ tr/-//;
    $reactioncpdquery =~s/\|+$//;
    
    
    my $command_ = "awk -F'\\t+' ' $reactioncpdquery {print \$2\"\\t\"\$1\"\\t\"\$3\"\\t\"\$4}' $nnnfile|sort|uniq";
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
    my %gnninfos=();
    undef @previousarray;
    my @gnninfoarray=();
    if ($gnnfile){
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
	
	$command_ = "awk -F'\\t+' ' $ec2genequery {print \$2\"\\t\"\$1\"\\t\"\$3\"\\t\"\$4}' \"$gnnfile\"|sort|uniq";
	&RSAT::message::TimeWarn($command_) if ($verbose >= 1);
	my @gnnonversiontable = qx ($command_);
	foreach my $content (@gnnonversiontable) {
# 	&RSAT::message::Info("EC2GENE:", $content) if ($verbose >= 3);
	    my @currentarray = split(/\t/,$content);
	    if ( @previousarray && !($previousarray[0] eq $currentarray[0])) {
		# 	  print $previousarray[0]."\n";
		
		my @truc = @gnninfoarray;
		&RSAT::message::Info("ECS2GENE:", $previousarray[0] ."\t" . $truc[0][1]) if ($verbose >= 3);
		$gnninfos{$previousarray[0]}=\@truc;
		undef @gnninfoarray;
	    }
	    &RSAT::message::Info("GENE:", $previousarray[0]) if ($verbose >= 3);
	    push (@gnninfoarray,\@currentarray);
	    @previousarray = @currentarray;
	}
	$gnninfos{$previousarray[0]}=\@gnninfoarray;
	
    }

    # End of Searching all reactions information for the reaction that are in  the infered pathway graph
    ################################################################

    ################################################################
    # Adding description to the pathway graph
    &RSAT::message::TimeWarn("Adding descriptions to pathway graph") if ($verbose >= 1);
    open (INFILE, '<'.$inputfile) || &RSAT::error::FatalError("Could not open file", $inputfile);
    open (OUTFILE, '>'.$outfile{graph_annot}) || &RSAT::error::FatalError("Could not write file", $outfile{graph_annot});
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
				if (%gnninfos){
				    chomp($ec);
				    my $genes = $gnninfos{$ec};
				    &RSAT::message::Info("EC:", "|".$ec."|",  $gnninfos{"$ec"}[0][1]) if ($verbose >= 4);
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
		    if ($values[0][3]){
			$tempdatab[3] =  $values[0][3];
		    }else {
			$tempdatab[3] = $values[0][1];
		    }
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
    &RSAT::message::TimeWarn("Converting graph to dot format in $outfile{graph_dot} ") if ($verbose >= 1);     
    my $convert_graph_cmd = $ENV{RSAT}."/perl-scripts/convert-graph -from path_extract -to dot -i $outfile{graph_annot} -o $outfile{graph_dot} $undirected";
    $convert_graph_cmd .= " -v $verbose" if ($verbose >= 1); 
    &RSAT::message::TimeWarn($convert_graph_cmd) if ($verbose >= 1);
    &RSAT::util::doit($convert_graph_cmd, $dry, $die_on_error, $verbose, $batch, $job_prefix);
    # End of Converting graph to dot graph format
    ################################################################

    ################################################################
    # Converting dot graph to image with graphviz dot
    &RSAT::message::TimeWarn("Creating graph image with graphviz dot") if ($verbose >= 1);
    my $dot_cmd = &RSAT::server::GetProgramPath("dot", 1, $ENV{dot_bin_dir});

    my $graph_image_cmd = $dot_cmd." -Tpng -Kdot -o $outfile{graph_png} $outfile{graph_dot}";
    &RSAT::message::TimeWarn($graph_image_cmd) if ($verbose >= 1);
    &RSAT::util::doit($graph_image_cmd, $dry, $die_on_error, $verbose, $batch, $job_prefix);

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
#   print  $exec_time if ($verbose); ## only report exec time if verbosity is specified
    
}
sub MapSeeds{

    ################################################################
    ## Initialise parameters
    #
    local $start_time = &RSAT::util::StartScript();
    $program_version = do { my @r = (q$Revision: 1.14 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r };
    #    $program_version = "0.00";
    my $query_ids;
    my @query_id_list;
    my ($input,
	$gnnfile,
	$nnnfile,
	$partial,
	$verbose) = @_;

#   if the parameters comes from the function use those one
#   $dir{output} = $outputdir if ($outputdir);    
#   $dir{temp} = $tempdir if ($tempdir);
#   $infile{gnn} = $gnnfile if ($gnnfile);
#   $infile{nnn} = $nnnfile if ($nnnfile); 
#   $infile{graph} = $graphfile if ($graphfile);
#  
    @query_id_list = split(/[\t\n;\|\r|\n|\t]+/,$input);

    ################################################################
    ## Print verbose
#   &Verbose() if ($verbose);

    ################################################################
    ## ECR Mapping

    ## DIDIER: the ECR mapping should be redone with the new program match-names.

    &RSAT::message::TimeWarn("Mapping seeds to reactions") if ($verbose >= 1);
    chomp(@query_id_list);
    if ( @query_id_list){
	
# $query_ids = (join "\$|^",@query_id_list );	# build a query from the input file or stdin
	my $separator = "\$|^";
	if ($partial){
	    $separator = "|";
	    # build a query from the input file or stdin
	}
	
	$query_ids = (join $separator,@query_id_list );
	&RSAT::message::Info("Query IDs", join("; ", @query_id_list)) if ($verbose >= 3);

	my @ercconversiontable;
	my %querylist=();
	## search into the GEC or GR file to find EC or Reactionid from gene input
	
	my $seed_converter_cmd = "awk -F'\\t+' '\$1~\"".$query_ids."\" {print \$2\"\\t\"\$1\"\\t\"\$3\"\\t\"\$4}' \"$gnnfile\"";

	&RSAT::message::TimeWarn("Seed conversion:", $seed_converter_cmd) if ($verbose >= 2);

	if($gnnfile){
	    my @gnnconversiontable = qx ($seed_converter_cmd) ;
	    chomp(@gnnconversiontable);
	    foreach my $line (@gnnconversiontable){
		@tempdata = split(/\t/,$line);
		$query_ids .= $separator.$tempdata[0]; # complete the query with id found in gnn file
	    }
	}
	# search in the ECR file with the complete query this normaly map ec, compound and reaction_id already presenet in the query
	$seed_converter_cmd = "awk -F'\\t+' '\$1~\"^".$query_ids."\" {print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4}' \"$nnnfile\"";
	&RSAT::message::TimeWarn("Seed conversion:", $seed_converter_cmd) if ($verbose >= 2);
	my @conversiontable = qx ($seed_converter_cmd) ;
	chomp(@conversiontable);
	&RSAT::message::Info("Mapped seeds\n", join("\n", @conversiontable)) if ($verbose >= 1);
	
# getting organism information
# End of ECR mapping
################################################################
################################################################
# preparing output data
	&RSAT::message::TimeWarn("Creating reaction file for pathway inference") if ($verbose >= 1);
#  my $outfile{seeds_converted} =$outputdir.(join "_",$localgroup_descriptor,$groupid,$graph, "_converted_seeds.txt");
	my $seednum= 0;
	my @previousarray;
	my %conversiontablehash= ();
	my %compounds = ();
	foreach my $val (@conversiontable) {
	    @tempdata = split(/\t/,$val);
	    #separating reaction nodes from compound nodes
	    if ($tempdata[2] =~/^compound.*|cpd.*/i){
		$compounds{$tempdata[1]} = $tempdata[0];
	    }else{
		my $ECgroup = $conversiontablehash{$tempdata[1]};
		if($ECgroup){
		    $ECgroup.="_$tempdata[0]";
		}else{
		    $ECgroup=$tempdata[0];
		}
		$conversiontablehash{$tempdata[1]} = $ECgroup;
	    }
#     @previousarray = @tempdata;
	}
	my %invertedconversiontablehash = ();
#   foreach my ($key,$val) (%conversiontablehash) {
	while (my ($key, $val) = each(%conversiontablehash)){

	    my $ECgroup = $invertedconversiontablehash{$val};
	    if($ECgroup){
		push (@{$ECgroup},$key);
	    }else{
		$ECgroup=[$key];
	    }
	    $invertedconversiontablehash{$val} = $ECgroup;
#      &RSAT::message::Info("$val:\n", join("\t", @{$ECgroup})) if ($verbose >= 1);
	    $seednum++;
	}
	
	return \%invertedconversiontablehash,\%compounds; 
    } else {
	&RSAT::message::Info("Empty Input\n");
	return ;
    }

}


################################################################
## Map query names to Network node ID using nnn et gnn files
sub QueryExactMetabNames{
    
    ################################################################
    ## Initialise parameters
    #
    local $start_time = &RSAT::util::StartScript();
    $program_version = do { my @r = (q$Revision: 1.14 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r };
    #    $program_version = "0.00";
    my $query_ids;
    my @query_id_list;
    my ($input,
	$gnnfile,
	$nnnfile,
	$partial,
	$verbose) = @_;

#   if the parameters comes from the function use those one
#   $dir{output} = $outputdir if ($outputdir);    
#   $dir{temp} = $tempdir if ($tempdir);
#   $infile{gnn} = $gnnfile if ($gnnfile);
#   $infile{nnn} = $nnnfile if ($nnnfile); 
#   $infile{graph} = $graphfile if ($graphfile);
#  
    @query_id_list = split(/[\t\n;\|\r|\n|\t]+/,$input);
    my @filteredlist =();
    foreach my $val (@query_id_list){
	push @filteredlist,$val if (length($val) > 2); 
    }
    @query_id_list = @query_id_list;

    ################################################################
    ## ECR Mapping

    ## DIDIER: the ECR mapping should be redone with the new program match-names.

#    &RSAT::message::TimeWarn("Searching info") if ($verbose >= 2);
    chomp(@query_id_list);
    
# $query_ids = (join "\$|^",@query_id_list );	# build a query from the input file or stdin
    my $separator = "\$|^";
    my $start = "^";
    my $end = "\$";
    if ($partial){
	$separator = "|";
	$start="";
	$end="";
    	# build a query from the input file or stdin
    }
    
    $query_ids = (join $separator,@query_id_list );
    &RSAT::message::Info("Query IDs", join("; ", @query_id_list)) if ($verbose >= 3);
    if ($partial){
	$separator = "|";
    	# build a query from the input file or stdin
    }
    my @gnnconversiontable;
    my %querylist=();

    ################################################################
    ## Seed conversion step 1:
    ## search into the GEC or GR file to find EC or Reactionid from gene input
    my $seed_converter_cmd = "awk -F'\\t+' '\$1~\"$start".$query_ids."$end\" {print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4}' \"$gnnfile\"";

    &RSAT::message::TimeWarn("Getting EC for gene seeds:", $seed_converter_cmd) if ($verbose >= 2);

    if($gnnfile){
	@gnnconversiontable = qx ($seed_converter_cmd) ;
	chomp(@gnnconversiontable);
	my $addedquery;
	foreach my $line (@gnnconversiontable){
	    @tempdata = split(/\t/,$line);
	    $addedquery .= "|^".$tempdata[1]."\$"; # complete the query with id found in gnn file
	}
	
	$query_ids .= $addedquery;
    }

    ################################################################
    ## Seed conversion step 2:
    ## Search in the ECR file with the complete query this normaly map
    ## ec, compound and reaction_id already present in the query.
    $seed_converter_cmd = "awk -F'\\t+' '\$1~\"$start".$query_ids."$end\" {print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4}' \"$nnnfile\"";
    &RSAT::message::TimeWarn("Getting reactions for seed ECs:", $seed_converter_cmd) if ($verbose >= 2);
    my @nnnconversiontable = qx ($seed_converter_cmd) ;
    chomp(@nnnconversiontable);

    &RSAT::message::Info("Mapped seeds\n", join("\n", @nnnconversiontable)) if ($verbose >= 3);
    
    my @conversiontable = @nnnconversiontable;
    @conversiontable = (@gnnconversiontable, @nnnconversiontable) if (@gnnconversiontable);

    # End of ECR mapping
    ################################################################


    ################################################################
    ## Preparing output data
    &RSAT::message::TimeWarn("Generating reaction table for pathway inference") if ($verbose >= 2);
    my $seednum= 0;
    my @previousarray;
    our %nnnconversiontablehash= ();
    my %compounds = ();
    foreach my $val (@nnnconversiontable) {
	@tempdata = split(/\t/,$val);
	#separating reaction nodes from compound nodes

	
	my $ECgroupref= $nnnconversiontablehash{$tempdata[0]};
	if($ECgroupref){
#       &RSAT::message::Info("ECgroupREFF ",$ECgroupref) if ($verbose >= 2);
#       &RSAT::message::Info("ECgroupBB ", $tempdata[0], $ECgroupref {"reactions"}) if ($verbose >= 2);
	    $nnnconversiontablehash {$tempdata[0]} {"ids"}.=", $tempdata[1]" if (!($nnnconversiontablehash {$tempdata[0]} {"ids"}=~m/$tempdata[1]/));
	    if ($tempdata[3]){
		$nnnconversiontablehash {$tempdata[0]} {"name"}.= ", $tempdata[3]";
	    }
	}else{
	    my %ECgroup=();
	    $ECgroup{"ids"} = $tempdata[1];
	    $ECgroup{"type"} = $tempdata[2];
	    if($tempdata[3]){
		$ECgroup{"name"} = $tempdata[3];
	    }else{
		$ECgroup{"name"} = $tempdata[1];
	    }
	    $nnnconversiontablehash{$tempdata[0]} = \%ECgroup; 
	    
	}
	&RSAT::message::Info("ECgroupA ", $tempdata[0],$nnnconversiontablehash {$tempdata[0]} {"ids"}) if ($verbose >= 3);
    }
    
#    my $testret = "";
#    
#    while (my($key, $val) = each(%nnnconversiontablehash)) { 
#      $testret .= $key ." | ". join("; ", @{$val})."\n" if ($val);
#    
#    }
#    
#   &RSAT::message::Info("TESTRETURN", $testret) if ($verbose >= 2);
    my $result = "";
    my @todelete =();
    foreach my $value (@gnnconversiontable) {
	my @tempdata = split(/\t/,$value);
	my $genemapping = $nnnconversiontablehash{"$tempdata[1]"};
	if ($genemapping){
	    my %hashreac = %{$genemapping};
	    $result .=  join ("\t", $tempdata[0], $tempdata[1], $hashreac{"ids"}, $tempdata[2], $hashreac{"name"});
	    $result .= "\n";
	    push (@todelete,  $tempdata[1]); 
	    &RSAT::message::Debug("GNN", "value=".$value) if ($main::verbose >= 4);
	}else {
	    &RSAT::message::Info("EC? ", $tempdata[1]) if ($verbose >= 3);
	}
    }
    # delete gnn mapping
    foreach my $gnn (@todelete){
	delete $nnnconversiontablehash{$gnn};
    }

    ## Add nnn mapping
    while (my($key, $value) = each(%nnnconversiontablehash)) {
	my %hashreac = %{$value};
	my $EC = "NA";
	$EC = $key if ($hashreac{"type"} eq "EC");
	$result .= join ("\t", $key, $EC, $hashreac{"ids"}, $hashreac{"type"}, $hashreac{"name"});
	$result .= "\n";	
	&RSAT::message::Debug("NNN", "key=".$key, "value=".$value) if ($main::verbose >= 4);
    } 
    &RSAT::message::Info("ret? ", $result) if ($verbose >= 3);

    return $result."\n";
}


__END__
