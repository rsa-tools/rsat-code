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
use File::Find;
use File::Basename;

use File::Temp qw/tempfile/;
use RSAT::PathwayExtraction;
use Cwd 'abs_path';
#  my $stderr = `$command 2>&1 1>/dev/null`;
#     $stderr = &error_handling($stderr, 1);
# 
#     if ($stderr =~ /Warning/) {
# 	die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: $stderr\ncommand: $command");
#     }
   


# use Pathwayinference;

  package  PathwayExtractor_WS;
  
  ###########################################
  ## Method : infer pathways
  ## Arguments :  
  ##   $class: calling class
  ##   $seeds: seed list
  ##   $organism: name of a supported organism
  ##   $network: name of a supported network !! file name format $SOURCENAME$/metab_$SOURCENAME_v$version$_$TYPE$_network.tab (etab_MetaCyc_v141_directed_network.tab)
  ## Web service for the pathway inference tool
  ##########################################
  sub infer_pathway {
  
  

#   return "RSATPATH: $ENV{RSAT}\nREA_ROOT: $ENV{REA_ROOT}\nKWALKS_ROOT: $ENV{KWALKS_ROOT}\nCLASSPATH:  $ENV{CLASSPATH}\n";
    my ($class, $seeds,$organism,$network,$directed) = @_;
    
    my $basedir="$ENV{RSAT}/public_html/data/metabolic_networks";
    my $GERdir = "$basedir/GER_files/$organism";
    my $outputdir = "$ENV{RSAT}/public_html/tmp/";
    my $outputurl = $ENV{rsat_ws_tmp};
    
         # build random id for file uniqueness
      my $filename;
      my (undef, $file) =  File::Temp::tempfile('XXXXXXXXX', OPEN=>0);
      $file = "PathwayExtractorWS$file";#fileprefix
 
    
    open(STDERR,">","$outputdir/$file.log");
    
    # get directory from network
    my $neworkfilepattern= $network;
    $neworkfilepattern =~ s/_v.*//gi; 
    $neworkfilepattern = "$basedir/networks/$neworkfilepattern/metab_$network";
    my $networknodenames = "$neworkfilepattern"."_node_names.tab";
    my $networkfile = "$neworkfilepattern"."_network.tab";
    
    print STDERR "networkfilepattern: $neworkfilepattern\n";
    
    # check parameters
    die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: $basedir\n does not exits!") unless -e $basedir;
    die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: Organism $organism not supported!") unless -e $GERdir;
    die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: Network $network not supported!") unless -e $networkfile;
    
     $neworkfilepattern .= "/metab_$network";

    # workaround try to guess REA_ROOT and KWALKS_ROOT if unable to get it from prperties
     $ENV{"REA_ROOT"} = "$ENV{RSAT}/contrib/REA" if (not exists($ENV{"REA_ROOT"}));
     $ENV{"KWALKS_ROOT"} = "$ENV{RSAT}/contrib/kwalks/kwalks/bin" if (not exists($ENV{"KWALKS_ROOT"}));
     $ENV{"CLASSPATH"} .= ":$ENV{RSAT}/java/lib/NeAT_javatools.jar" unless ( exists($ENV{"CLASSPATH"})&& !($ENV{"CLASSPATH"}=~ /NeAT_javatools.jar/));
#     
      print STDERR "RSATPATH: $ENV{RSAT}\n";
      print STDERR "REA_ROOT: $ENV{REA_ROOT}\n";
      print STDERR "KWALKS_ROOT: $ENV{KWALKS_ROOT}\n";
      print STDERR "CLASSPATH: $ENV{CLASSPATH}\n";
    
    
    
    &RSAT::message::TimeWarn("running inferpathway:START") if ($verbose >= 2);
    SOAPAction: "running inferpathway:START";
      $outputdir = Cwd::abs_path($outputdir);
      print STDERR "outputdir: $outputdir\n";
      #my $outputfile = PathwayExtraction::InferpathwaySub(
      my $outputfile = &RSAT::PathwayExtraction::Inferpathway(
#       "NP_416523.1\tNP_416524.1\tNP_416525.1\tNP_416526.4\tNP_416527.1\tNP_416528.2\tNP_416529.1\tNP_416530.1",
      $seeds,
      0,
      $outputdir,
      "$GERdir/$organism"."-gene_refseq_ec.tab",
      $networknodenames,
      $networkfile,
      $directed,
      $outputdir, "$file",3,{'-Q' => " ",'-J' => " "});
      &RSAT::message::TimeWarn("running inferpathway:END") if ($verbose >= 2);
      &RSAT::message::TimeWarn("processing files:START in $outputdir") if ($verbose >= 2);
     
     
      &RSAT::PathwayExtraction::ProcessOutputFiles(
      $outputfile,
      $outputdir,
      "$GERdir/$organism"."-gene_refseq_ec.tab",
      $networknodenames,$directed,0);
      &RSAT::message::TimeWarn("processing files:END") if ($verbose >= 2); 
      
      my $cmd = "zip -j $outputdir/$file"."_results $outputdir/$file*";
      print STDERR $cmd."\n";
      print STDERR qx($cmd)."\n";

    my $filename = $file."_results.zip";
    my $fileurl = "$outputurl/$filename";
    open RESFILE, "<$outputdir/$filename";
    my $resultfile =  do { local $/; <RESFILE> };
    close (RESFILE);
    close STDERR;
    
    return SOAP::Data->name('response' => {'url' => $fileurl,
 				        'filename' => $filename,	
				        'file' => $resultfile});
  }
 
  sub hi{
   return "hello world\n";
  }
#   sub returnhash{
#   open RESFILE, "</home/rsat/rsa-tools//public_html/tmp/wqQQQE_results.zip";
#  my $resultfile =  do { local $/; <RESFILE> };
# #      
# my $resultfile =<RESFILE>;
#      close (RESFILE);
# #   my $resultfile = `cat /home/rsat/rsa-tools//public_html/tmp/wqQQQE_results.zip`;
#    return SOAP::Data->name('response' => {'url' => "hello world\n",
#    				        'filename' => "wqQQQE_results.zip\n",
# 				        'file' =>$resultfile});
#   }
#   
#   sub convertfile{
#     my ($class, $filename) = @_;
#     
#     open FHDL, $filename or print STDERR "unable to read file $filname";
#       my ($buffer, $data, $n);
#       while (($n = read FHDL, $buffer, 4) != 0) {
# # 	print "$n bytes read\n";
# 	$data .= $buffer ;
#       }
#       close FHDL;
#       return $data; 
#   }
#   
 ###########################################
  ## Method : get_supported_organisms
  ## Arguments :  none
  ## Web service for  getting the list of supported organisms by infer_pathways
  ##########################################
  sub get_supported_organisms{
  open(STDERR,">","$ENV{RSAT}/public_html/tmp/stdERR.log");
    my $dirtoget="$ENV{RSAT}/public_html/data/metabolic_networks/GER_files/";
   print STDERR "Directory: $dirtoget\n";
    opendir(IMD, $dirtoget) || die("Cannot open directory");
    @thefiles= readdir(IMD);
    my $ret;
    foreach my $file(@thefiles){
	next unless -d $dirtoget.$file;
# 	$file = $dirtoget.$files;
 	if ($file ne "." && $file ne ".."){
	  $ret.= $file."\n";
	  print STDERR $file."\n";
 	}
    }
    close STDERR;
    closedir(IMD); 
    return $ret;
  }
  
  ###########################################
  ## Method : get_supported_networks
  ## Arguments :  none
  ## Web service for  getting the list of supported networks by infer_pathways
  ##########################################
  sub get_supported_networks{
    my $dirtoget="$ENV{RSAT}/public_html/data/metabolic_networks/networks/";
    my @file_list =  &find_files($self,"$ENV{RSAT}/public_html/data/metabolic_networks/networks/","_network\\.tab\$");
    foreach(@file_list)
    {
      s/_network.tab$//g;
      s/metab_//g;
    }
    
    my $ret = join (',', @file_list);
    print STDERR $ret;  
    return $ret;
  }
  
  ###########################################
  ## Method : find_files
  ## Arguments : 
  ##     $class: caller class
  ##     $dirtoget: directory to starting
  ##     $pattern: regex pattern
  ##     $absolute: if true return absulute path
  ## Find recusivelly files starting in $dirtoget
  ##########################################
  sub find_files{
  my ($class,$dirtoget,$pattern,$absolute) = @_;
  open(STDERR,">","$ENV{RSAT}/public_html/tmp/stdERR.log");
    print STDERR "Directory: |$absolute|\n";
    my @file_list=();   
       File::Find::find(
       sub { 
	 if (-f && /$pattern/){
	   print STDERR $_."\t". $File::Find::name."\n" ;
	   if ($absolute){
	     push @file_list, $File::Find::name;
	   }else{
	     push @file_list, $_;
	   }
	    	   
	 }
       }, $dirtoget);
    close STDERR;
    return @file_list;
  }
  
  
#   {
#      print $ENV{RSAT}."\n";
# #     print &get_supported_organisms()."\n";
#  print &get_supported_networks()."\n";
# 
#   } 
# 
  1;
  