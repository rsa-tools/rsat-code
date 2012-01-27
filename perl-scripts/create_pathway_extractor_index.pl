#!/usr/bin/perl -w
############################################################
#
# $Id: create_pathway_extractor_index.pl,v 1.2 2012/01/27 10:39:19 jvanheld Exp $
#
############################################################

## use strict;

=pod

=head1 Extract Pathway from gene list

=head1 VERSION 1.0

=head1 DESCRIPTION

This script creates a gene index file from the results of
I<footprint-discovery>, I<infer-operons>, and pathway inference
analysis.

=head1 AUTHORS

didier.croes@ulb.ac.be

=head1 CATEGORY

Graph tool

=head1 USAGE

create_pathway_extractor_index -h -of operonfile -fd footprintdirectory [-o outputdirectory] [-v verbosity] -a gene2ec2reactions [-b gene2reaction] 

=head1 infer operon file format
; infer-operon  -v 1 -sep - -return q_info,operon,gene_nb -org Escherichia_coli_K_12_substr__MG1655_uid57779 -min_gene_nb 2 -all -dist 55
; Organism      Escherichia_coli_K_12_substr__MG1655_uid57779
; Queries       4430
; Column descriptions
;       1       query           Query string.
;       2       id              ID of the query gene.
;       3       name            Name of the query gene.
;       4       strand          Strand of the query gene.
;       5       start           Starting position of the query gene (first nucleotide of the start codon).
;       6       end             Ending position of the query gene (last nucleotide of the stop codon).
;       7       operon          Member genes of the operon.
;       8       gene_nb         number of genes in the predicted operon
#query  id      name    strand  start   end     operon  gene_nb description
NP_414543.1     NP_414543.1     thrA    D       337     2799    thrA-thrB-thrC  3       bifunctional: aspartokinase I (N-terminal); homoserine dehydrogenase I (C-terminal); fused a
spartokinase I and homoserine dehydrogenase I
NP_414544.1     NP_414544.1     thrB    D       2801    3733    thrA-thrB-thrC  3       homoserine kinase
NP_414545.1     NP_414545.1     thrC    D       3734    5020    thrA-thrB-thrC  3       threonine synthase
NP_414552.1     NP_414552.1     yaaW    R       11356   10643   yaaI-yaaW       2       conserved protein, UPF0174 family
NP_414554.1     NP_414554.1     yaaI    R       11786   11382   yaaI-yaaW       2       conserved protein, UPF0412 family
NP_414562.1     NP_414562.1     insB    R       20314   19811   insA-insB       2       IS1 transposase B

=head1  -a gene2ec2reactions [-b gene2reaction] file format
  gene_id ec_number       reaction_id     species_name    taxonomy_id     gene_name
  O22340  4.2.3.- RXN-10482       Abies grandis   46611   (4S)-limonene synthase
  O22340  4.2.3.- RXN-10483       Abies grandis   46611   (4S)-limonene synthase
  O22340  4.2.3.- RXN-10566       Abies grandis   46611  my %reactioninfos=(); (4S)-limonene synthase
  O22340  4.2.3.- RXN-10567       Abies grandis   46611   (4S)-limonene synthase
  O22340  4.2.3.- RXN-10568       Abies grandis   46611   (4S)-limonene synthase
  O22340  4.2.3.- RXN-10600       Abies grandis   46611   (4S)-limonene synthase

=head1 example:

B<pathway_extractor.pl> -i met_genes -g MetaCyc_directed_141.txt  -b METACYC_GPR_NOEC.tab -a METACYC_GPR_EC.tab -o temp

cat met_genes|B<pathway_extractor.pl> -g MetaCyc_directed_141.txt  -b METACYC_GPR_NOEC.tab -a METACYC_GPR_EC.tab -o temp


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
    local $start_time = &RSAT::util::StartScript();
    $program_version = do { my @r = (q$Revision: 1.2 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r };
#    $program_version = "0.00";

    $main::operondir;				# File name containing infer operon results
    $main::footprintsdir;				# File name containing infer operon results
    $main::selectedorganismfile;
    $main::outfile = "";			# output directory
    $main::outdir = "";			# output directory
    $main::verbose = 0;
    $main::in = STDIN;
    $main::out = STDOUT;
    
    $main::gprfile ="METACYC_GPR_EC.tab";	# GPR Gene -> EC -> REACTION annotation file path. Default (METACYC_GPR_EC.tab)
    $main::grfile="";				# GR Gene -> REACTION annotation

    ################################################################
    ## Read argument values
    &ReadArguments();

    ################################################################
    ## Check argument values
    
    
#     if (defined $main::infile{input}) {
#       ($main::in) = &OpenInputFile($main::infile{input});
   
    ################################################################
    ## Read input
#     ($main::in) = &OpenInputFile($main::infile{input});
#     while (<$main::in>) {
# 
#     }
#     

    ################################################################
    ## Print verbose
    &Verbose() if ($main::verbose);

    ################################################################
    ## Execute the command
     print "selectedorganismfile:".$selectedorganismfile ."\n"; 
    open SELECTEDORG, "<$main::selectedorganismfile";  
    
    while (my $soligne = <SELECTEDORG>) {
      next if ($soligne =~ /^;/);
      next if ($soligne =~ /^#/);
      next if ($soligne =~ /^$/);
      chomp $soligne;
     
      my @orgtax = split /\s+|\t/, $soligne;
      my $organism = "Unknown";
      my $organismid;
      my $rsatorganismid = $orgtax[0];
      my $rsattaxon=$orgtax[1];
       print "ORGANISM:" ."$rsatorganismid\t$rsattaxon"."\n";
      # my $working_dir = "";
      my $genesid;

       #load infer operon infered pathways
      print STDERR  "READING infer operon results\n";
      my %operon2filehash = ();
      my $operonpathdir = $operondir.$rsatorganismid."/pathres/";
      my $searchoperonpngcmd = "du -a $operondir$rsatorganismid|grep .png|sort -n";
      print "searchoperonpngcmd = $searchoperonpngcmd\n";
      my @operonimages= qx($searchoperonpngcmd);
      my %operon2pp = ();
      foreach  my $val (@operonimages) {
	@contents = split /\t/, $val;
	if ($contents[0] > 4){
	  my $temp= $contents[1];
	  chomp($temp);
	  $temp=~ s/$operonpathdir//;
	  @operonarray = split /_/, $temp;
	  $operonstrng=$operonarray[0];
	  chomp($operonstrng);
	  $temp=~ s/\.png$//;
	  push @{$operon2pp{$operonstrng}}, $temp;
 	  print $operonstrng."\t$temp\n";
	  
	}
      }

      my $operonresultfile = $operondir.$rsatorganismid."/".$rsatorganismid."_operons_inferred_dist55_2genes.tab"; 
      print "operonfile:".$operonresultfile."\n";
      open INFEROPERON, "<$operonresultfile";
      my %gene2operonhash = ();
      #  ; Organism      Escherichia_coli_K_12_substr__MG1655_uid57779
      while (my $ligne = <INFEROPERON>) {
# 	  print "OPERON:".$ligne;
	  chomp $ligne;
	  if ($ligne =~ /^;/){
	    my @comments = split /\s+|\t/, $ligne;
	    next;
	  };
	  next if ($ligne =~ /^#/);
	  next if ($ligne =~ /^$/);
	 
	  my @lignecp = split /\t/, $ligne;
	  if (!$organismid){
	    my $findorgcmd = "awk '/$lignecp[0]/ { print \$_; exit }' $main::gprfile";
	    
	    my @orgarray = qx($findorgcmd);
	    if ($orgarray[0]){
	      @orgarray = split /\t/,$orgarray[0];
	      $organismid = $orgarray[4];
	      print "taxid = " . $organismid."\n";;
	    }
	  }
	  push @{$gene2operonhash{$lignecp[0]}}, @lignecp;
      }
      close (INFEROPERON);
      print STDERR  "READING GPR\n";
      my @grconversiontable;
      if ($grfile){
	my $seed_converter_cmd = "awk -F'\\t+' '\$5==\"".$organismid."\" {print \$_}' $grfile";
	@grconversiontable = qx ($seed_converter_cmd) ;
	chomp(@grconversiontable);
	@grconversiontable = sort(@grconversiontable);
# 	print join( "\n",@grconversiontable) ;
      }
       
      my $seed_converter_cmd = "awk -F'\\t+' '\$5==\"".$organismid."\" {print \$_}' $gprfile";
      print $seed_converter_cmd."\n";

      my @conversiontable = qx ($seed_converter_cmd) ;
      chomp(@conversiontable);
      @conversiontable =sort(@conversiontable);
      # getting organism information
      my @tempdata = split(/\t/,$conversiontable[0]);
      $organism = $tempdata[3];
      $organismid = $tempdata[4];
#       my $groupid=( join "-",$organism,$organismid);
#       $groupid=~s/(\s|\(|\))+/_/g;
      #  merging GR and GPR in one array
      if (@grconversiontable){
	foreach  my $val (@conversiontable) {
	  @tempdata = split(/\t/,$val);
	  if (!grep(/^$tempdata[0]/, @grconversiontable)){
	    push(@grconversiontable, grep(/^$tempdata[0]/, @conversiontable));
# 	    print join( "\n",@conversiontable) ."\n"; 
	  }
	}
      #   print join( "\n",@conversiontable) ."\n"; 
	@conversiontable=@grconversiontable;
      }
      my %gene2info = ();
      my @reacinfoarray =();
#        foreach my $content(@conversiontable){
# 	my @currentarray = split(/\t/,$content);
# 	
# 	    $gene2info{$currentarray[0]}=\@currentarray;
# 	
# 	  
#        }
      foreach my $content(@conversiontable){
	my @currentarray = split(/\t/,$content);
      # 	if ($currentarray[0] eq "NP_391805.2"){
      # 	print join( "\t",@currentarray) ."\n";
      # 	}
	if ( @previousarray && !($previousarray[0] eq $currentarray[0])){ 
      # 	    print $previousarray[0]."\n";
	    my @truc = @reacinfoarray;
	    $gene2info{$previousarray[0]}=\@truc;
	    undef @reacinfoarray;
	}
	push (@reacinfoarray,\@currentarray);
	  @previousarray = @currentarray;
      }
	$gene2info{$previousarray[0]}=\@reacinfoarray;
	
      #load footprints infered pathways
      print STDERR  "load footprints infered pathways\n";
      my $footpathwaysdir = $footprintsdir."/predicted_pathways/".$rsattaxon."/".$rsatorganismid."/";
      my $searchpngcmd = "du -a $footpathwaysdir|grep .png|sort -n";
      print "DUCMD = $searchpngcmd\n";
      my @images= qx($searchpngcmd);
      my %gene2fp =();
      foreach  my $val (@images) {
	@contents = split /\t/, $val;
	if ($contents[0] > 4){ #to replace by 4 when released
	  my $temp= $contents[1];
	  chomp($temp);
	  $temp=~ s/$footpathwaysdir//;
	  $temp=~ s/\/.*//;
	  
	  push @{$gene2fp{$temp}}, $contents[1];
 	  print $temp.",";
	  
	}
      }
#   exit 0;
       
      print STDERR "Writing output";
      #        while ( my ($key, $value) = each(%gene2info) ) {
      $outdir = "$operondir$rsatorganismid/index.tab";
  open OUTFILE, ">$outdir";
      
  print OUTFILE "; infer-operon  -v 1 -sep - -return q_info,operon,gene_nb -org Escherichia_coli_K_12_substr__MG1655_uid57779 -min_gene_nb 2 -all -dist 55\n";
  print OUTFILE "; Organism\t$organismid\n";
  print OUTFILE "; Column descriptions\n";
  print OUTFILE ";\t1\t[Gene_ID]\tquery gene identifier.\n";
  print OUTFILE ";\t2\t[Gene_NAME]\tname of the query gene.\n";
  print OUTFILE ";\t3\t[Protein_name]\tname of the protein.\n";
  print OUTFILE ";\t4\t[EC numbers]\tEC Numbers.\n";
  print OUTFILE ";\t5\t[Reaction_id]\tReaction_id.\n";
  print OUTFILE ";\t6\t[op]\tthe operon predicted.\n";
  print OUTFILE ";\t7\t[opp]pathway_predicted\tthe pathway predicted.\n";
  print OUTFILE ";\t8\t[fp]phylogenetic footprints (conserved cis-regulatory element)\n";
  print OUTFILE ";\t9\t[co-fp]\tpredicted co-regulation set (genes with the same footprints)\n";
  print OUTFILE ";\t10\t[fppath]\tpathways predicted from footprints (one predicted pathway for each co-regulation neighbour)\n";
  print OUTFILE ";\t11\t[fppathproj]\tfootprint-based discovered pathway projected against METACYC\n";
  print OUTFILE ";\t12\t[Description]\tDescription\n";
  print OUTFILE "#Gene_ID\tGene_NAME\tProtein_name\tEC_numbers\tReaction_id\tOperon\topp\tfp\tco-fp\fppathproj\n";
      
      for my $key ( keys %gene2info ) {
      #         print "$key =>";
      	foreach my $content1($gene2info{$key}){
	  my @level1 = @{$content1};
	  if (@level1){
	    
	    my $operons= $gene2operonhash{$key};
	    my @newoperon;
	    my $operonfileprefix;
	    my $operonfile;
	    if ($operons){
	      @newoperon = @{$operons};
	      print "OPERON".join( "\t",@newoperon) ."\n";
	       $operonfileprefix= $operon2pp{$newoperon[6]}[0];
	       if ($operonfileprefix){
		 $operonfile=$operonpathdir.$operonfileprefix.".txt";
	       }
	    }	  
	    my $footprint= $gene2fp{$key};
	    my @newfp;
	    if ($footprint){
	      @newfp = @{$footprint};
	      print join( "\t",@newfp) ."\n";	  
	    }
# 	    
	    if ($operons||$footprint){
	      print OUTFILE "$key\t";
	      my $test1= $level1[0];
	      my @test2= @{$test1}; 
# 	      print OUTFILE "$test2[5]\t\t";
	      print OUTFILE "$level1[0][5]\t\t";
	      
	      my @reactions;
	      my @ecs = ();
	      foreach my $content1(@level1){
		my @geneinfoarray = @{$content1};	
# 		print ">>>>>>>>".join( "\t",@geneinfoarray) ."\n";
		my @grepres;
		
		if ($operons){
		  if ($operonfileprefix){
		    my $grepcmdop = "grep \"^$geneinfoarray[2]\[<>\].*Reaction\" $operonfile"; 
# 		    print $grepcmdop."\n";
		    push @grepres, qx ($grepcmdop);
		  }
		}
		my $ppfyfile = $newfp[0];
		if ($footprint){
		  my $grepcmdfp = "grep \"^$geneinfoarray[2]\[<>\].*Reaction\" $ppfyfile"; 
		  print $grepcmdfp."\n";
		  push  @grepres, qx ($grepcmdfp);
		}
		if ($grepres[0]){
		  chomp($geneinfoarray[2]);
		  chomp($geneinfoarray[1]);
		  push @reactions, $geneinfoarray[2];
		  push @ecs, $geneinfoarray[1];
	# 	print $currentarray[1]."\n";
		}
	      } 
	      print OUTFILE join("<BR>",@ecs)."\t"; 
	      print OUTFILE join("<BR>",@reactions)."\t";
	      
	      if ($operons ) {
		print OUTFILE "$newoperon[6]\t";
		if($operonfileprefix){
		  print OUTFILE "<a href=pathres/$operonfileprefix".".png>png</a> <a href=pathres/$operonfileprefix".".txt>tab</a>\t";
		}else {
		  print OUTFILE "\t";
		}
	      }else {
		 print OUTFILE "\t\t";
	      }
	      if ($footprint) {
		  my @splitpath = split /\//,$newfp[0];
		  my $ppfyfile = $splitpath[@splitpath-1];
		  $ppfyfile =~ s/\.png$//;
# 		  print ">>>>".$ppfyfile;
		  $ppfyfile= "../../footprints/predicted_pathways/$rsattaxon/$rsatorganismid/$key/$ppfyfile";
		chomp  ($ppfyfile);
# 		print ">>>>".$ppfyfile;
	      	print OUTFILE "<a href=../../footprints/$rsattaxon/$rsatorganismid/$key/$key"."_$rsatorganismid"."_$rsattaxon"."_index.html>footprints</a>\t<a href=$ppfyfile".".png>png</a> <a href=$ppfyfile".".txt>tab</a>\t\n";
	      }else {
	       print OUTFILE "\t\t\t\n";
	      }
# 		//print OUTFILE "\n";
		print "*********************\n";	

	      
	    }
	  }
	}
	
      }
      close OUTFILE;
	################################################################
	## Insert here output printing

	################################################################
	## Report execution time and close output stream
	my $exec_time = &RSAT::util::ReportExecutionTime($start_time); ## This has to be exectuted by all scripts
	print $main::out $exec_time if ($main::verbose >= 1); ## only report exec time if verbosity is specified
	close $main::out if ($main::outfile{output});
    }
   exit 0;
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

=item B<-v #>

Level of verbosity (detail in the warning messages during execution)

=cut
    if ($arg eq "-v") {
      if (&IsNatural($arguments[0])) {
	$main::verbose = shift(@arguments);
      } else {
	$main::verbose = 1;
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

=item B<-fd>

Footprint discovery results directory

=cut
    } elsif ($arg eq "-fd") {
      $main::footprintsdir= shift(@arguments);
=pod

=item B<-of>

execute gwenview to display the pathway results in png format
=cut
    } elsif ($arg eq "-od") {
     $main::operondir  = shift(@arguments);
=pod

=item	B<-s Selected organisms file>

selected organism to annotate

=cut
    } elsif ($arg eq "-s") {
      $main::selectedorganismfile = shift(@arguments);
    
=pod

=item	B<-a GPR Genes file Default (METACYC_GPR_EC.tab)>

GPR Gene -> EC -> REACTION annotation file path. Default (METACYC_GPR_EC.tab)

=cut
    } elsif ($arg eq "-a") {
      $main::gprfile = shift(@arguments);
    
=pod

=item	B<-b GR Gene -> REACTION annotation>

An gene annotation file with diredt link gene to reaction. Does not rely on the EC number annotation

=cut
    } elsif ($arg eq "-b") {
      $main::grfile = shift(@arguments);
=pod


=item	B<-o output directory>

If no output direcory is specified, the stdin directory is used. 

=cut
    } elsif ($arg eq "-o") {
      $main::outdir = shift(@arguments);
 
    } else {
      &FatalError(join("\t", "Invalid option", $arg));

    }
  }

=pod

=back

=cut

}
################################################################
## Verbose message
sub Verbose {
    print $main::out "; template ";
    &PrintArguments($main::out);
    printf $main::out "; %-22s\t%s\n", "Program version", $program_version;
    if (%main::infile) {
	print $main::out "; Input files\n";
	while (my ($key,$value) = each %main::infile) {
	  printf $main::out ";\t%-13s\t%s\n", $key, $value;
	}
    }
    if (%main::outfile) {
	print $main::out "; Output files\n";
	while (my ($key,$value) = each %main::outfile) {
	  printf $main::out ";\t%-13s\t%s\n", $key, $value;
	}
    }
}
__END__
