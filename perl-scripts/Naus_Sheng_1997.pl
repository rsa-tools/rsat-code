#!/usr/bin/perl

if (($ARGV[0] eq "-help") || ($ARGV[0] eq "-h")) {
  open HELP, "| more";
  print HELP <<End_of_Help;
NAME
	Naus_Sheng_1997.pl
	
	July 1997 by Jacques van Helden

CATEGORY
    statistics

DESCRIPTION
	This scripts calculates the approximations according to 
	   Naus, J.I. and Sheng, K. (1997). 
	   Matching among multiple random sequence. 
	   Bulletin of Mathematical Biology 59(3): 483-496.
USAGE
	Naus_Sheng_1997.pl -w word_length -T seq_length
		-R #_matching_seq  [-S #_seq] [-wp pattern_proba] 
				
OPTIONS
	-w #	word length
	-wp #	word probability. 
	-T #	sequence length.
	-S #	number of matching sequences. 
		If omitted, S takes the same value as R.
	-R #	number of sequences matching the word
	-bino	a binomial estimation is used instead of Naus & Sheng 
		formula.
	-v	verbose. Displays the details of calculation.
	
End_of_Help
  exit(0);
}
  
#### read arguments ####
foreach $a (0..$#ARGV) {

  if ($ARGV[$a] eq "-w") {
    $w = $ARGV[$a+1];

  } elsif ($ARGV[$a] eq "-R") {
    $R = $ARGV[$a+1];

  } elsif ($ARGV[$a] eq "-S") {
    $S = $ARGV[$a+1];

  } elsif ($ARGV[$a] eq "-T") {
    $T = $ARGV[$a+1];

  } elsif ($ARGV[$a] eq "-wp") {
    $wp = $ARGV[$a+1];

  } elsif ($ARGV[$a] eq "-v") {
    $verbose = 1;

  } elsif ($ARGV[$a] eq "-bino") {
    $bino = 1;
  }
}


unless (defined $w) {
  print "	You should define either word length.\n";
  print "	type Naus_Sheng_1997.pl -h for more info.\n";
  exit;
}


unless (defined $R) {
  print "	You should define the number of sequences.\n";
  print "	type Naus_Sheng_1997.pl -h for more info.\n";
  exit;
}

unless (defined $T) {
  print "	You should define sequence length.\n";
  print "	type Naus_Sheng_1997.pl -h for more info.\n";
  exit;
}

unless (defined $S) {
  $S = $R;
}

if ($S < $R) {
  print "	The number of matching sequences (-R) cannot be higher\n";
  print "	than the total number of sequences (-S).\n";
  print "	type Naus_Sheng_1997.pl -h for more info.\n";
  exit;
}

$m = 4;
$V = $m**$w;
$lp = 1/$m;
unless (defined $wp) {
  $wp = $lp**$w;
}


$Tprime = $T +1 - $w;


#### binomial estimate ####
if (($bino) || ($verbose)) {
  $delta_bino = 1 - `binomial -p $wp -r $Tprime -s 0`;
  $Pspec = `binomial -p $delta_bino -r $S -s $R`;
  chomp($Pspec);
  $Pspec_boe = `binomial -p $delta_bino -r $S -s $R -boe`;
  chomp($Pspec_boe);
  $Pany =1- `binomial -p $Pspec -r $V -s 0`;
}

#### Naus & Sheng approximation ####
if ((!$bino) || ($verbose)) {
  $delta = 1-exp( -$Tprime * $wp);
  $lambda = (1/$m)**($R-1);
  $Tstar = $V*(1-(($V-1)/$V)**$Tprime);
  $P8 = 1- exp(-$Tstar * $delta**($R-1) * (1 - $lambda));
  $P10 = `binomial -p $delta -r $S -s $R`;
  chomp($P10);
  $P12 = 1- (1 - `binomial -p $delta -r $S -s $R`)**$V;
}

#### print result ####
if ($bino) {
  print "$Pany\n";
} else {
  if ($R==$S) {
    print "$P8\n";
  } else {
    print "$P12\n";
  }
}

if ($verbose) {
  print <<End_verbose;

Parameters:
    w	$w	(word length)
    T	$T	(sequence length)
    R	$R	(number of sequences)
    S	$S	(matching sequences)
    
Detail of calculation:
    m	$m	(number of letters in the alphabet) 
    V	$V	(number of possible words of length $w)
    lp	$lp	(probability of each letter)
    wp	$wp	(probability of a $w letter word)
    T'	$Tprime	(number of possible word positions in a sequence)
    
Binomial estimate:
    delta	$delta_bino	(proba for a given word in 1 sequence)
    Pspec	$Pspec	(proba for a given word in $R/$S sequences)
    Pspec boe	$Pspec_boe	(proba for a given word in at least $R/$S sequences)
    Pany	$Pany	(proba for any word in $R of $S sequences)

Naus & Sheng approximations:
    T*		$Tstar	
    delta	$delta	(proba for a given word in 1 sequence)
    lambda	$lambda	
    P(8)	$P8	(proba for any word in $S/$S sequences)
    P(10)	$P10	(proba for a given word in $R/$S sequences)
    P(12)	$P12	(proba for any word in $R/$S sequences)
End_verbose
}


exit(0);



