#!usr/bin/perl
##Programa para convertir salida TSV de RSAT feature maps a JSON jerarquico

use strict;
use warnings;
use File::Path qw(make_path);

## Modify input and output
my $input_file = "fm_input.tsv";
my $output_dir = "feature_maps_jsons";
my $output_file = "prueba_2";

## hash to keep all info
my %info_hash;

## Open TSV file
open(INPUT, '<' , $input_file);
my @feature_maps = (<INPUT>);
close(< INPUT >);

## Un contador, las primeras dos lineas me dan info general del general

my $cont = 0;
foreach(@feature_maps){
   ##Split each line  y tabs
   my @fields = split(/\t/,$_);
   my $gene = $fields[0];
   my $next = $fields[2];

   print "$gene\n";

   if($fields[2] =~ /^SEQ_START/){
      print "Start: $fields[5]\n";
      my @seq_start = split(/\-/,$fields[4]);
      $info_hash{$gene}{"SEQ_START"} = $seq_start[1];
   }
   elsif($fields[2] =~ /^SEQ_END/){
      print "End: $fields[5]\n";
      my @seq_end = split(/\-/, $fields[5]);
      $info_hash{$gene}{"SEQ_END"} = $seq_end[1];
   }

   else{
      my $seq = $fields[2];
      ##print "Seq_1: $seq\n";
      $info_hash{$gene}{"Features"}{$cont}{"Seq"} = $seq;

      my $strand = $fields[3];
      ##print "Strand: $strand\n";
      $info_hash{$gene}{"Features"}{$cont}{"Strand"} = $strand;

      my $start = $fields[4];
      ##print "Start: $start\n";
      my @start_array = split(/\-/,$start);
      $info_hash{$gene}{"Features"}{$cont}{"Start"} = $start_array[1];

      my $end = $fields[5];
      ##print "End: $end\n";
      my @end_array = split(/\-/, $end);
      $info_hash{$gene}{"Features"}{$cont}{"End"} = $end_array[1];

      my $seq_2 = $fields[6];
      ##print "Seq_2: $seq_2\n";
      $info_hash{$gene}{"Features"}{$cont}{"Seq_2"} = $seq_2;

      my $score = $fields[7];
      ##print "Score: $score\n";
      $info_hash{$gene}{"Features"}{$cont}{"Score"} = $score;

      $cont = $cont + 1;
   }
}

my $g_cont = 0;


foreach my $gene (keys %info_hash){
   my $f_cont = 0;
   foreach my $feature(keys %{$info_hash{$gene}{"Features"}}){
      $f_cont = $f_cont + 1;
   }
   $g_cont = $g_cont + 1;
   $info_hash{$gene}{"Features_cont"} = $f_cont;
}

## Print Output

## Create directory defined at the start of the program
make_path( $output_dir );

# # ## Create and open a file for each gene on the vcf
open(my $OUTPUT, '>', $output_dir."/".$output_file.".json");

print $OUTPUT "[\n"; ## Opening group of genes

my $gene_cont = 0;
foreach my $gene_name(keys %info_hash){

   print $OUTPUT "\t"."{"."\n";## Opening individual gene
   print $OUTPUT "\t\t".'"Gene_name":'.'"'.$gene_name.'"'.","."\n";
   print $OUTPUT "\t\t".'"Start":'.$info_hash{$gene_name}{'SEQ_START'}.","."\n";
   print $OUTPUT "\t\t".'"End":'.$info_hash{$gene_name}{'SEQ_END'}.","."\n";
   print $OUTPUT "\t\t".'"Features":'."\n";

   print $OUTPUT "\t\t"."["."\n";## Opening group of features

   my $feature_cont = 0;
   foreach my $feature(keys %{$info_hash{$gene_name}{"Features"}}){
      print $OUTPUT "\t\t"."{"."\n";
      print $OUTPUT "\t\t"."\t\t".'"seq":'.'"'.$info_hash{$gene_name}{"Features"}{$feature}{"Seq"}.'"'.","."\n";
      print $OUTPUT "\t\t"."\t\t".'"strand":'.'"'.$info_hash{$gene_name}{"Features"}{$feature}{"Strand"}.'"'.","."\n";
      print $OUTPUT "\t\t"."\t\t".'"start":'.$info_hash{$gene_name}{"Features"}{$feature}{"Start"}.","."\n";
      print $OUTPUT "\t\t"."\t\t".'"end":'.$info_hash{$gene_name}{"Features"}{$feature}{"End"}.","."\n";
      print $OUTPUT "\t\t"."\t\t".'"seq_2":'.'"'.$info_hash{$gene_name}{"Features"}{$feature}{"Seq_2"}.'"'.","."\n";
      print $OUTPUT "\t\t"."\t\t".'"score":'.$info_hash{$gene_name}{"Features"}{$feature}{"Score"}."\n";

      print $feature_cont. "de" .$info_hash{$gene_name}{"Features_cont"}."\n";

      my  $fe_cont = $info_hash{$gene_name}{"Features_cont"};

      if($feature_cont < $fe_cont-1){
         print $OUTPUT "\t\t"."},"."\n";
      }
      else{
         print $OUTPUT "\t\t"."}"."\n";
      }
      $feature_cont = $feature_cont + 1;
   }

   print $OUTPUT "\t\t"."]"."\n";

   if($gene_cont < $g_cont-1){ ##Closing individual genes
      print "Gene ". $gene_cont ."de". $g_cont."\n";
      print $OUTPUT "},"."\n";
   }
   else{
      print $OUTPUT "}"."\n";
   }
   $gene_cont = $gene_cont +1;
}
   print $OUTPUT "]"."\n";## Closing group of genes
