#!/usr/bin/perl -w
require "RSA.lib";
BEGIN {
   if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
 	push (@INC,"$ENV{RSAT}/perl-scripts/lib/");
 	push (@INC,"$ENV{RSAT}/perl-scripts/");
  }
}
# require "RSA.lib";
# use lib '../../perl-scripts/lib/';


use RSAT::util;
use File::Temp qw/tempfile/;
use RSAT::PathwayExtraction;
use Cwd 'abs_path';

# use Pathwayinference;

  package  PathwayExtractor_WS;
  sub inferpathway {

#   return "RSATPATH: $ENV{RSAT}\nREA_ROOT: $ENV{REA_ROOT}\nKWALKS_ROOT: $ENV{KWALKS_ROOT}\nCLASSPATH:  $ENV{CLASSPATH}\n";

    my $basedir="/home/rsat/rsa-tools/public_html/data/metabolic_networks";
    my $outputdir = "$ENV{RSAT}/public_html/tmp/";
    my $outputurl = $ENV{rsat_ws_tmp};
    open(STDERR,">","$outputdir/stdERR.log");
    
    # workaround try to guess REA_ROOT and KWALKS_ROOT if unable to get it from prperties
     $ENV{"REA_ROOT"} = "$ENV{RSAT}/contrib/REA" if (not exists($ENV{"REA_ROOT"}));
     $ENV{"KWALKS_ROOT"} = "$ENV{RSAT}/contrib/kwalks/kwalks/bin" if (not exists($ENV{"KWALKS_ROOT"}));
     $ENV{"CLASSPATH"} .= ":$ENV{RSAT}/java/lib/NeAT_javatools.jar" if (not exists($ENV{"CLASSPATH"})|| !($ENV{"CLASSPATH"}=~ /NeAT_javatools.jar/));
#     
    print STDERR "RSATPATH: $ENV{RSAT}\n";
    print STDERR "REA_ROOT: $ENV{REA_ROOT}\n";
    print STDERR "KWALKS_ROOT: $ENV{KWALKS_ROOT}\n";
    print STDERR "CLASSPATH: $ENV{CLASSPATH}\n";
    
    my ($class, $seeds) = @_;
    
    &RSAT::message::TimeWarn("running inferpathway:START");# if ($verbose >= 2);
    SOAPAction: "running inferpathway:START";
    # build random id for file uniqueness
    my $filename;
    (undef, $file) =  File::Temp::tempfile('XXXXXX', OPEN=>0);
#      return "starting in $basedir \n";
#     my $outputfile = &RSAT::Pathwayinference::Inferpathway(
 

      $outputdir = Cwd::abs_path($outputdir);
      print STDERR "outputdir: $outputdir\n";
      #my $outputfile = PathwayExtraction::InferpathwaySub(
      my $outputfile = &RSAT::PathwayExtraction::Inferpathway(
#       "NP_416523.1\tNP_416524.1\tNP_416525.1\tNP_416526.4\tNP_416527.1\tNP_416528.2\tNP_416529.1\tNP_416530.1",
      $seeds,
      0,
      $outputdir,
      "$basedir/GER_files/Escherichia_coli_strain_K12_gene_refseq_ec.tab",
      "$basedir/GER_files/MetaCyc_EC_cpds_2rxns.tab",
      "$basedir/networks/MetaCyc_directed_141.txt",
      $outputdir, "WS$file",3);
#       
        
      &RSAT::message::TimeWarn("running inferpathway:END");# if ($verbose >= 2);
#      --return  $outputfile;
      &RSAT::message::TimeWarn("processing files:START in $outputdir");# if ($verbose >= 2);
      &RSAT::PathwayExtraction::ProcessOutputFiles(
      $outputfile,
      $outputdir,
      "$basedir/GER_files/Escherichia_coli_strain_K12_gene_refseq_ec.tab",
      "$basedir/GER_files/MetaCyc_EC_cpds_2rxns.tab",0);
      &RSAT::message::TimeWarn("processing files:END");# if ($verbose >= 2); 
      
      my $cmd = "zip -j $outputdir/$file"."_results $outputdir/WS$file*";
      print STDERR $cmd."\n";
      print STDERR qx($cmd)."\n";

    return $outputurl."/$file"."_results.zip";     
  }
 
  sub hi{
   return "hello world\n";
  }
  
  sub InferpathwaySub{

  ################################################################
  ## Initialise parameters
  #
   
  local $start_time = &RSAT::util::StartScript();
  $program_version = do { my @r = (q$Revision: 1.6 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r };
  #    $program_version = "0.00";
   my $query_ids;
   my @query_id_list;
   my ($input,
       $isinputfile,
       $outputdir,
       $gnnfile,
       $nnnfile,
       $graphfile,
       $tempdir,
       $localgroup_descriptor,
       $verbose,
       $piparameters) = @_;
  
  my %localotherPIparameters = %{$piparameters} if ($piparameters);
  
 
  if ($isinputfile){
    ($in) = &RSAT::util::OpenInputFile($input);
     
    @query_id_list = <$in>;
    close $in;
    #   if no group descriptor use the input file name
    if (!$localgroup_descriptor) {
      $localgroup_descriptor = $input;
      $localgroup_descriptor =~ s{.*/}{};     # removes path
      $localgroup_descriptor=~ s{\.[^.]+$}{}; # removes extension
    }
    
  }elsif ($input){
     @query_id_list = split(/\t|\s|;/,$input);
  }else{
    @query_id_list = <$in>;
  }

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



  ################################################################
  ## ECR Mapping

  ## DIDIER: the ECR mapping should be redone with the new program match-names.

  &RSAT::message::TimeWarn("Mapping seeds to reactions") if ($verbose >= 1);
  
  chomp(@query_id_list);
  $query_ids = (join "\$|^",@query_id_list );	# build a query from the input file or stdin
  &RSAT::message::Info("Query IDs", join("; ", @query_id_list)) if ($verbose >= 3);

  my @ercconversiontable;
  my %querylist=();
  ## search into the GEC or GR file to find EC or Reactionid from gene input
  my $seed_converter_cmd = "awk -F'\\t+' '\$1~\"^".$query_ids."\" {print \$2\"\\t\"\$1\"\\t\"\$3\"\\t\"\$4}' \"$gnnfile\"";
    
  &RSAT::message::TimeWarn("Seed conversion:", $seed_converter_cmd) if ($verbose >= 2);

  if($gnnfile){
    my @gnnconversiontable = qx ($seed_converter_cmd) ;
    chomp(@gnnconversiontable);
    foreach my $line (@gnnconversiontable){
      @tempdata = split(/\t/,$line);
      $query_ids .= "\$|^".$tempdata[0]; # complete the query with id found in gnn file
    }
  }
  # search in the ECR file with the complete query this normaly map ec, compound and reaction_id already presenet in the query
  $seed_converter_cmd = "awk -F'\\t+' '\$1~\"^".$query_ids."\" {print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4}' \"$nnnfile\"";
   &RSAT::message::TimeWarn("Seed conversion:", $seed_converter_cmd) if ($verbose >= 2);
  my @conversiontable = qx ($seed_converter_cmd) ;
  chomp(@conversiontable);
  &RSAT::message::Info("Predicted pathway file\n", join("\n", @conversiontable)) if ($verbose >= 1);
  
# getting organism information
# End of ECR mapping
################################################################
  if (! defined $group_id) {
  $groupid =$localgroup_descriptor;
  }
  $groupid=~s/(\s|\(|\))+/_/g;
  
  ################################################################
  ## Define output file names
  $outfile{gr} = "";		# GR Gene -> REACTION annotation
  $outfile{prefix} =$outputdir."/";
  $outfile{prefix} .= join("_", $localgroup_descriptor, $groupid, $graph, "pred_pathways");
  $outfile{prefix} =~ s|//|/|g; ## Suppress double slashes

  $outfile{predicted_pathway} = $outfile{prefix}.".txt";
  $outfile{seeds_converted} = $outfile{prefix}."_seeds_converted.txt";
  ################################################################
  # Creating reaction file fo pathway inference
  &RSAT::message::TimeWarn("Creating reaction file for pathway inference") if ($verbose >= 1);

  open (MYFILE, '>'.$outfile{seeds_converted});
  
  print MYFILE "# Batch file created", &RSAT::util::AlphaDate();
  print MYFILE "# ECR groups parsed from file:". "\n";
  print MYFILE "# EC number grouping: true". "\n";
  my $seednum= 0;
  my @previousarray;
  foreach my $val (@conversiontable) {

    @tempdata = split(/\t/,$val);
    if (@previousarray && !($tempdata[0] eq $previousarray[0])) {
      print MYFILE "$previousarray[0]\t$groupid\n";
      $seednum++;
    }
    # 	print "$tempdata[1] eq $previousarray[1]\n";
    print MYFILE $tempdata[1] .">\t".$tempdata[1]. "\n";
    print MYFILE $tempdata[1] ."<\t".$tempdata[1]. "\n";
    print MYFILE $tempdata[1] ."\t".$tempdata[0]. "\n";
    @previousarray = @tempdata;
  }

  if (@conversiontable) {
    print MYFILE "$previousarray[0]\t$groupid\n";
    $seednum++;
  }
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

    our $minpathlength = $otherPIoptions{"-m"} || "5";
    delete($otherPIoptions{"-m"});
    our $graphfileformat = $otherPIoptions{"-f"} || "flat";
    delete($otherPIoptions{"-f"});
    our $weightpolicy= $otherPIoptions{"-y"} || "con";
    delete($otherPIoptions{"-h"});
    delete($otherPIoptions{"-g"});
    delete($otherPIoptions{"-o"});
    delete($otherPIoptions{"-p"});
    delete($otherPIoptions{"-v"});
    #after the program has handled the mandatory parameters, it will add the remaining ones 
    my $piparameters =" ";
    
    while( my($key, $val) = each(%localotherPIparameters) ) {
	$piparameters .= "$key $val ";
    }
    &RSAT::message::Info("RSATPATH: ", $ENV{RSAT}) if ($verbose >= 1);
    &RSAT::message::Info("REA_ROOT: ", $ENV{REA_ROOT}) if ($verbose >= 1);
    &RSAT::message::Info("KWALKS_ROOT: ", $ENV{KWALKS_ROOT}) if ($verbose >= 1);
    &RSAT::message::Info("CLASSPATH: ",$ENV{CLASSPATH}) if ($verbose >= 1);
    
#     print STDERR "#########################################################\n";
#     print STDERR "RSATPATH: $ENV{RSAT}\n";
#     print STDERR "c $ENV{REA_ROOT}\n";
#     print STDERR "KWALKS_ROOT: $ENV{KWALKS_ROOT}\n";
#     print STDERR "#########################################################\n";
#     
    my $pathway_infer_cmd = "java -Xmx1000M graphtools.algorithms.Pathwayinference";
    $pathway_infer_cmd .= " -i ".$outfile{seeds_converted};
    $pathway_infer_cmd .= " -m $minpathlength -C -f $graphfileformat";
    $pathway_infer_cmd .= " -p $tempdir";
    $pathway_infer_cmd .= " -E $outputdir";
#     $pathway_infer_cmd .= " -E temp -p temp "; 
    $pathway_infer_cmd .= " -d -b";
    $pathway_infer_cmd .= " -g $graphfile";
    $pathway_infer_cmd .= " -y $weightpolicy -v $verbose ";
    $pathway_infer_cmd .=  $piparameters;
    $pathway_infer_cmd .= " -o $outfile{predicted_pathway}";
    $pathway_infer_cmd .= " -v  true 1> temp/test.log";
#     $pathway_infer_cmd .= " -A  $ENV{REA_ROOT}";
#     $pathway_infer_cmd .= " -K  $ENV{KWALKS_ROOT}";
    &RSAT::message::TimeWarn("Pathway inference command", $pathway_infer_cmd) if ($verbose >= 2);
    &RSAT::message::Info("Predicted pathway file", $outfile{predicted_pathway}) if ($verbose >= 1);
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
  
   my $exec_time = &RSAT::util::ReportExecutionTime($start_time); ## This has to be exectuted by all scripts

#   print STDERR $exec_time if ($verbose); ## only report exec time if verbosity is specified
#   close $out if ($outfile{output});

 return $outfile{predicted_pathway};
}

#   {
#      print $ENV{PERL5LIB}."\n";
#     print &inferpathway()."\n";
#   } 

  1;
  