#!/usr/bin/perl

##############################################################
###
### Statistics
###
##############################################################

package RSAT::stats;

use RSAT::GenericObject;
use RSAT::message;
use RSAT::error;
use RSAT::util;
@ISA = qw( RSAT::GenericObject );

=pod

=head1 NAME

RSAT::stats

=head1 DESCRIPTION

A method class with statistical procedures and functions.

=head1 METHODS

=over 4

=cut


=pod

=item permute(@array)

Permute the elements of an array

=cut

sub permute {
    my (@array) = @_;
    my @permuted = ();
    my $n  = scalar(@array);
    &RSAT::message::Warning(join("\t", "Permuting array of ", scalar(@array),"elements")) if ($main::verbose >= 3);
    for my $i (0..$#array) {
	my $rand = int(rand($n-$i));
	push @permuted, splice(@array, $rand, 1);
	&RSAT::message::Debug("&RSAT::stats::permute()", $i, "rand=".$rand, $#array, $#permuted)  if ($main::verbose >= 10);
    }
    return @permuted;
}

=pod

=item B<summary>

Calculate basic statistics on a list of values

=cut
sub summary {
    my (@values) = @_;
    my %stats = (n=>0,
		 sum=>0,
		 mean=>0,
		 median=>0,
		 var=>0,
		 sd=>0,
		 min=>"NA",
		 max=>"NA",
		);
    my $n = scalar(@values);
    $stats{n} = $n;
    $stats{sum} = &sum(@values);
    if ($n > 0) {

	## Calculate the mean
	$stats{mean} = $stats{sum}/$n;

	## Calculate the median
	$stats{median} = &median(@values);

	## Calculate standard deviation
	my $ssq = 0;
	foreach my $value (@values) {
	    $ssq += $value*$value;
	}
	$stats{ssq} = $ssq;
	$stats{var} = $ssq/$n - $stats{mean}*$stats{mean};
	$stats{sd} = sqrt($stats{var});

	$stats{var_est} = $n/($n-1)*$stats{var};
	$stats{sd_est} = sqrt($stats{var_est});

	## Extremes
	$stats{min} = &min(@values);
	$stats{max} = &max(@values);
    }
    return %stats;
}


##############################################################
## Return the sum of a list of numeric values
## usage:
## $sum = &RSAT::stats::sum(@value_list);
sub sum {
    my (@values) = @_;
    my $sum = 0;
    foreach my $value (@values) {
	$sum += $value;
    }
    return $sum;
}


##############################################################
## Return the mean of a list of numeric values
## usage:
## my $mean = &RSAT::stats::mean(@value_list);
sub mean {
    my (@values) = @_;
    my $mean = 0;
    my $n = scalar(@values);
    if ($n == 0) {
      $mean = "NA";
      &RSAT::message::Warning("Cannot calculate the mean of an empty list");
    } else {
      $mean = &RSAT::stats::sum(@values)/$n;
    }
    return $mean;
}


##############################################################
## Return the median of a list of numerical values
sub median {
  my $n = scalar(@_);
#  &RSAT::message::Debug("&RSAT::stats::median()", $n." values") if ($main::verbose >= 10);
  my $median = "NA";
  if ($n >= 1) {
    my @sorted_values = sort {$a <=> $b} @_;
    if ($n % 2 == 0) {
      $median = ($sorted_values[$n/2-1] + $sorted_values[$n/2])/2;
    } else {
      $median = $sorted_values[($n-1)/2];
    }
  }
  return $median;
}

## Return the maximum of a list of numeric values
## usage:
## $max_value = &RSAT::stats::max(@value_list);
sub max {
    my @sorted_values = sort {$a <=> $b} @_;
    return $sorted_values[$#sorted_values];
}


## Return the minimum of a list of numeric values
## usage:
## $min_value = &RSAT::stats::min(@value_list);
sub min {
    my @sorted_values = sort {
	$a <=> $b
	} @_;
    return $sorted_values[0];
}

## Return the minimum of a list of numeric values
## usage:
## $min_value = &RSAT::stats::min(@value_list);
sub checked_min {
    my @sorted_values = sort {
	$a <=> $b
	} @_;
    my $c;
    $c = -1;
    do {
	$c++;
    } until (($c > $#sorted_values) || (&RSAT::util::IsReal($sorted_values[$c])));
    return $sorted_values[$c];
}

## Return the minimum of a list of numeric values
## usage:
## $min_value = &RSAT::stats::min(@value_list);
sub checked_max {
    my @sorted_values = sort {
	$b <=> $a
	} @_;
    my $c;
    $c = -1;
    do {
	$c++;
    } until (($c > $#sorted_values) || (&RSAT::util::IsReal($sorted_values[$c])));
    return $sorted_values[$c];
}

## Return the sum of a list of numbers
## ignore non-numeric values
sub checked_sum {
    my @values = @_;
    my $sum = 0;
    foreach $value (@values) {
	$sum += $value if (&RSAT::util::IsReal($value));
    }
    return $sum;
}

## Return the sum of a list of numbers
## ignore non-numeric values
sub checked_avg {
    my @values = @_;
    my $sum = 0;
    my $valid_values = 0;
    foreach $value (@values) {
	if (&RSAT::util::IsReal($value)) {
	    $valid_values++;
	    $sum += $value;
	}
    }
    if ($valid_values > 0) {
	return ($sum/$valid_values);
    } else {
	&RSAT::error::FatalError( "Error: cannot calculate average of an empty list\n");
    }
}

##############################################################
### Usage:
###     &RSAT::stats::binomial($proba,$trials,$successes)
###
### This routine uses the recursive formula for the calculation of the
### binomial density function.
###
###           P(x) * p(r-x)
###  P(x+1) = -------------
###              q(x+1)
###
### We calculate everything in logarithms, and make a custom exp
### conversion to overcome the limitation of the Perl exp(x)
### function. This allows to compute small numbers as low as ~1e-300,
### instead of the classical 1e-15 limit of computation for floating
### point numbers.
###
sub binomial {
    my($proba, $trials, $succ) = @_;
    my($q) = 1 - $proba;

    if ($proba <=0) {
	&RSAT::error::FatalError( "Error : invalid probability value $proba\n");
    }

    my($logproba) = log($proba);
    my($logq) = log($q);

    my($x) = 0;
    my($logbin) = log($q)*($trials-$x);

    for $x (0..$succ-1) {
        $logbin += $logproba + log($trials - $x) - $logq - log($x+1);
    }

    my $bin = &LogToEng($logbin);
    return ($bin);
}


##############################################################
## Converts a natural logarithm (x) into e^x.
##
## This is performed by generating the output string, in order to
## cirumvent a problem with the standard Perl function exp(), which
## returns 0 for  values < e-322.
sub LogToEng {
  my ($log) = @_;
  my $base = 10;
  my $log_base = log($base);
  $log /= log(10);
  $eng = 10**(1+$log - int($log));
  $eng .= "e";
  if (int($log)-1 > 0) {
    $eng .= "+";
  }
  $eng .= int($log)-1;
  return($eng);
}

##############################################################
### usage:
###     &RSAT::stats::sum_of_binomials($proba,$trials,$from,$to)
### Calculates the sum of binomial probabilities between two values.
###
### This routine uses the recursive formula for the calculation of
### binomial:
###           P(x) * p(r-x)
###  P(x+1) = -------------
###              q(x+1)
###
sub sum_of_binomials {
    my ($proba, $trials, $from, $to) = @_;
    if ($main::verbose >= 5) {
      &RSAT::message::Debug("sum_of_binomials", "p=$proba", "t=$trials", "from=$from", "to=$to");
    }

    my $expected = $trials*$proba;

    ## Precision limits: Perl sometimes returns a very smalml negative value (-1e-15) instead of 0
    my $precision_limit = 1.0e-14;
    if (($proba < 0) && (abs($proba) < $precision_limit)) {
      &RSAT::message::Warning("RSA.stat.lib", "sum_of_binomials", "p", $proba, "Rounding proba to 0") if ($main::verbose >= 3);
      $proba = 0;
    }
    if (($proba > 1) && ($proba -1 < $precision_limit)) {
      &RSAT::message::Warning("RSA.stat.lib", "sum_of_binomials", "p", $proba, "p-1", $proba-1, "Rounding proba to 1") if ($main::verbose >= 3);
      $proba = 1;
    }

    ## Impossible probability values
    if ($proba <0) {
	&RSAT::error::FatalError("Error: invalid probability $proba (must be a positive value)\n");
    }
    if ($proba > 1) {
      &RSAT::error::FatalError("RSA.stat.lib", "sum_of_binomials",  $proba, "Invalid value for p, should be <= 1");
    }

    ($from,$to) = sort {$a <=> $b} ($from, $to);
    if ($to > $trials) {
	&RSAT::error::FatalError("&RSAT::stats::binomial()", "Successes ($to) cannot be higher than trials ($trials)\n");
    }


    ## Limit cases : the result is trivial, so there is no need to
    ## compute sum for obtaining the result
    if ($proba == 0) {
	if ($from == 0) {
	    return 1;
	} else {
	    return 0;
	}
    }
    if ($proba == 1) {
	if ($to == $trials) {
	    return 1;
	} else {
	    return 0;
	}
    }


    ## initialize
    my $q = 1 - $proba;
    my $logproba = log($proba);
    my $logq = log($q);
    my $sum_of_bin = 0;
    my $x = 0;
    my $logbin;
    my $prev_value = -1;

    $x = 0;
    $logbin = $trials * log($q);
    if ($from == 0) {
	$sum_of_bin = exp($logbin);
	$loop_start = 1;
    } else {
	$loop_start = $from;
	for $x (1..$from-1) {
	    $logbin += log($trials - $x + 1) - log($x) + $logproba - $logq;
	}
    }

    for $x ($loop_start..$to) {
	$logbin += log($trials - $x + 1) - log($x) + $logproba - $logq;
	$sum_of_bin += exp($logbin);
	last if (($x > $expected) && ($sum_of_bin <= $prev_value)); ### limit of precision
	$prev_value = $sum_of_bin;
    }

    $sum_of_bin = 1 if ($sum_of_bin > 1);
    $sum_of_bin = "ND0" if ($sum_of_bin < 0);
    return($sum_of_bin);
}


################################################################
## Usage:
##    my $p_val = &RSAT::stats::binomial_boe($proba,$trials,$successes);
##
## Calculate the probability of observing >=s successes with a
## probability $p and $r repeats, with the sum of binomials for j
## varying from s to r.
sub binomial_boe {
    my ($proba, $trials, $succ) = @_;
    if ($main::verbose >=5) {
	print join ("\t", "binomial_boe", "p=$proba", "t=$trials", "s=$succ"), "\n";
    }
    return (sum_of_binomials($proba,$trials,$succ,$trials));
}


##############################################################
## usage: &binomial_approx($proba,$trials,$successes)
##
## This is the entropy approximation for the sum of binomials.
## Note that the approximation is only valid for s/r > p.
sub binomial_approx {
    my ($proba, $trials, $successes) = @_;
    my $beta = $successes/$trials;
    my $alpha = $proba;

    if (($beta >= $alpha) &&
	($beta <1) &&
	($beta > 0)) {
	$b_approx = exp(-$trials*($beta*log(abs($beta/$alpha)) + (1-$beta)*log(abs((1-$beta)/(1-$alpha)))));
    } else {
	$b_approx = "ND";
    }
    $b_approx;
}


##############################################################
## Negative binomial distribution
##
## usage:
##    negbin($successes,$p, $k)
##    negbin2($successes,$mean, $variance)
##		  s        s    k+s
##	P(X=s) = C      * p  / q
##  	          k+s-1
##
##  	 where	s is the number of successful trials,
##		p is a real value comprized between 0 and 1
##		q = 1 + p
##
##	This distribution is aggregative, with a mean
##	     m = kp
##	and a variance
##	     v = kpq
##
##	Instead of p and k, the mean (m) and variance (v) can be
##	provided; p and k are then calculated from these parameters as
##	follows.
##
##		q = v/m
##		p = 1-q = 1-v/m
##		k = m/p = 1/(v = m)
##
##  COMPUTATION
##
##	The probability is calculated with the recursive formula.
##
##	                p(k+x)
##	    P(X=s+1) = --------*P(X=s)
##	                q(x+1)

##############################################################
## Calculate negbin with mean and variance as parameters
sub negbin2 {
    my ($x, $mean, $variance, $series) = @_;

    ### Check parameters
    &RSAT::error::FatalError("negbin2: the mean should be strictly positive") if ($mean <= 0);
    &RSAT::error::FatalError("negbin2: the variance should be positive") if ($variance < 0);
    &RSAT::error::FatalError("negbin2: the variance should be greater than the mean") if ($variance <= $mean);

    ### Convert mean and variance to p and k
    my $p = $variance/$mean -1;
    my $k = $mean/$p;
    &RSAT::message::Info("Calculating p and k from mean and variance",
	  ";\tp = ", $p, "\n",
	  ";\tk = ", $k, "\n"
	 ) if ($main::verbose >= 6);

    ### Calculate the negbin
    my @negbin = &negbin($x, $p, $k, $series);
    return ($p, $k, @negbin);

}

##############################################################
## Calculate negbin with p and k as parameters
sub negbin {
    my ($x,$p, $k, $series) = @_;

    my $q = 1 + $p;
    my $log_q = log($q);
    my $log_p = log($p);
    my $log_k = log($k);
    my @log_negbins = ();

    &RSAT::message::Info("Calculating negative binomial\n",
	  "; p = ", $p, "\n",
	  "; q = ", $q, "\n",
	  "; k = ", $k, "\n",
	  "; x = ", $x, "\n",
	 ) if ($main::verbose >= 6);

    ## Calculate proba for 0 successes
    my $log_negbin = - $log_q*$k;;
    push @log_negbins, $log_negbin if ($series);

    ## Use recursive formula
    for my $i (1..$x) {
	$log_negbin += $log_p + log($k+$i - 1) - $log_q - log($i);
	if ($series) {
	    push @log_negbins, $log_negbin; ## value for x==$i
	}
    }

    if ($series) {
	my @negbins = ();
	foreach my $log_negbin (@log_negbins) {
	    my $negbin = &LogToEng($log_negbin);
	    push @negbins, $negbin;
#	    warn join ("\t", $log_negbin, $negbin), "\n" if ($main::verbose >= 10);
	}
	return @negbins;
    } else {
	my $negbin = &LogToEng($log_negbin);
	return $negbin;
    }
}

##############################################################
## Calculate negbin with mean and variance as parameters
sub sum_of_negbin2 {
  my ($mean, $variance, $from, $to) = @_;

  ### Check parameters
  &Warning("Cannot calculate negbin with mean <= 0 and x >0") if (($mean <= 0) && ($to > 0));
  &RSAT::error::FatalError("negbin2: the mean should be strictly positive") if ($mean <= 0);
  &RSAT::error::FatalError("negbin2: the variance should be positive") if ($variance < 0);
  &RSAT::error::FatalError("negbin2: the variance should be greater than the mean") if ($variance <= $mean);

  ### Convert mean and variance to p and k
  my $p = $variance/$mean -1;
  my $k = $mean/$p;
  warn ("; Calculating p and k from mean and variance\n",
	";\tp = ", $p, "\n",
	";\tk = ", $k, "\n"
       ) if ($main::verbose >= 6);

  ### Calculate the negbin
  my $sum_of_negbin = &sum_of_negbin($p, $k, $from, $to);
  return ($p, $k, $sum_of_negbin);
}

##############################################################
## sum_of_negbin($lambda, $from, $to)
## Calculates the Negbin probability for a given interval of values
sub sum_of_negbin {
  my ($p, $k, $from,$to) = @_;
  my $sum_of_negbin;
  my $prev_sum;

  my $q = 1 + $p;
  my $log_q = log($q);
  my $log_p = log($p);
  my $log_k = log($k);
  my @log_negbins = ();

  warn ("; Calculating negative binomial\n",
	";\tp = ", $p, "\n",
	";\tq = ", $q, "\n",
	";\tk = ", $k, "\n",
	";\tfrom = ", $from, "\n",
	";\tto = ", $to, "\n",
       ) if ($main::verbose >= 10);

  ## Calculate proba for 0 successes
  my $log_negbin = - $log_q*$k;;
  $sum_of_negbin = &LogToEng($log_negbin) if ($from == 0);

  ## Use recursive formula
  for my $i (1..$to) {
    $prev_sum = $sum_of_negbin;
    $log_negbin += $log_p + log($k+$i - 1) - $log_q - log($i);
    if ($i >= $from) {
      $sum_of_negbin += &LogToEng($log_negbin);
    }
    warn join ("\t", $i, $log_negbin, &LogToEng($log_negbin), $sum_of_negbin), "\n"
      if ($main::verbose >= 10);
    last if (($sum_of_negbin > 0) && ($prev_sum == $sum_of_negbin));
  }

  return $sum_of_negbin;
}


##############################################################
## Poisson distribution
##
## Usage
##   poisson($x,$lambda)
##
## Where
##   $x is the observed number of successes
##   $lambda is the number of successes expected by chance (the mean
##           of the Poisson distribution)
##
## When the option $series is not null, the function computes the 3
## distributions (density, CDF and P-value) for the whole range from 0
## to x.
##
## note: on our sun station, this algorithm works only for m < 746
##
## Note: the distribution is computed on the basis of a recursive
## formula rather than the raw formula. This computation increases
## precision and time efficiency.
##
## Raw formula:
##    p(x) = lambda^x exp(-lambda)/x!
## Recursive formula:
##    p(x) = p(x-1) * lambda / x
##
## The implementation of this raw formula is iterative rather than
## recursive to avoid embedded function calls.
##
## Equivalent R command for validation of the whole series:
##   p=0.0001; n=500; s=120; m=p*n; data.frame(i=0:s, dpois=dpois(0:s,lambda=m),ppois=ppois(0:s,lambda=m),pval=ppois((0:s)-1,lambda=m,lower=F))
##
sub poisson {
    my ($x,$lambda, $series) = @_;
    if ($lambda <= 0) {
	&RSAT::error::FatalError($lambda." is not a valid value for the expected mean of &poisson() . Must be strictly positive.");
    }

    my $log_lambda = log($lambda);
    my $log_dpois;
    my @log_dpois = ();

    if ($x == 0) {
	$log_dpois = -$lambda;
	push @log_dpois, $log_dpois if ($series);
    } else {
	## value for x==1
	$log_dpois = -$lambda + $log_lambda;

	if ($series) {
	    push @log_dpois, -$lambda; ## value for x==0
	    push @log_dpois, $log_dpois; ## value for x==1
	}
    }

    for my $i (2..$x) {
      $log_dpois += $log_lambda - log($i);
      if ($series) {
	push @log_dpois, $log_dpois; ## value for x==$i
      }
    }


    ## Compute at once the density + cumulative density function (CDF)
    if ($series) {
      my @dpois = (); ## Poisson density
      my @ppois = (); ## Cumulative density function (CDF)

      foreach my $i (0..$#log_dpois) {
#	my $dpois = &LogToEng($log_dpois[$i]);
	$dpois = exp($log_dpois[$i]) || &LogToEng($log_dpois[$i]); ## The function &LogToEng is used here only if exp returns 0
	push @dpois, $dpois;
	$ppois[$i] = $dpois;
	$ppois[$i] += $ppois[$i-1] if ($i >= 1); ## update CDF
#	&RSAT::message::Debug($i, $dpois, $ppois[$i]) if ($main::verbose >= 10);
      }

      ################################################################
      ## Compute the P-value

      ## We first need to compute the right tail of the distribution.
      ## The Poisson distribution extends to infinite, but we can stop
      ## computing when the summed terms are below the precision
      ## limit.
      my @vpois = (); ## P-value, defined as P(X >= x)
      my $dpois = 0;
      my $last_dpois = $dpois[$#dpois];
      $i = $#dpois;
      my $after_last = 0;
      do {
	$i++;
	$log_dpois += $log_lambda - log($i);
	$dpois = exp($log_dpois) || &LogToEng($log_dpois); ## We use the tricky function &LogToEng only when the Perl exp() function fails (returns 0)
	push @dpois, $dpois;
#	&RSAT::message::Debug($i, $dpois, $last_dpois, $last_dpois - $dpois) if ($main::verbose >= #0);

	## The basic idea for stopping cumulating is that the addition
	## of a new term does not change the sum. This is however not
	## sufficient: when P-values are very small, adding the next
	## density term would not change the sum whereas adding the
	## whole remaining right tail up to infinity would change
	## it. We thus keep computing (and summing) the terms for 20
	## more steps. Checked by comparing the result with R
	if ($last_dpois + $dpois == $last_dpois) {
	  $after_last++;
	}
      } until ($after_last > 20);

      ################################################################
      ## Due to the imprecision of the Perl sum, we loose precision
      ## around 1e-285
      foreach my $j (0..$#dpois) {
	my $i = $#dpois - $j;
	$vpois[$i] = $dpois[$i];
	$vpois[$i] += $vpois[$i+1] if ($i < $#vpois);
      }

      ################################################################
      ## Truncate the arrays to return only ther equired range (the
      ## vector shad to be extended for computing the P-value).
      @dpois = @dpois[0..$x];
      @vpois = @vpois[0..$x];

      my $dpois_ref = \@dpois;
      my $ppois_ref = \@ppois;
      my $vpois_ref = \@vpois;
#      &RSAT::message::Debug($dpois_ref, $ppois_ref) if ($main::verbose >= 10);
      return ($dpois_ref, $ppois_ref, $vpois_ref);

    } else {
      my $dpois = &LogToEng($log_dpois);
      return $dpois;
    }
}

##############################################################
## sum_of_poisson($lambda, $from, $to)
## Calculates the Poisson probability for a given interval of values
sub sum_of_poisson {
  my ($lambda,$from,$to) = @_;
#    my $dpois;
  if (($lambda <= 0) && (($from > 0) || ($to > 0))){
    &Warning( "sum_of_poisson($lambda, $from, $to). Cannot calculate Poisson probability with lambda <= 0 and x > 0");
    return("NA");
  }
  my $log_lambda = log($lambda);

  my $log_dpois;
  my $sum_of_poi;
  my $start;

  $sum_of_poi = exp(-$lambda) if ($from ==0);

  if ($to ==0) {
    $log_dpois = -$lambda;
  } else {
    $log_dpois = -$lambda + $log_lambda;
    $sum_of_poi += exp($log_dpois) if ($from <= 1);
  }


  if ($from <= 1) {
    $start = 2;
  } else {
    $start = $from;
    for my $i (2..$start-1) {
      $log_dpois += $log_lambda - log($i);
    }
  }

  for my $i ($start..$to) {
    $log_dpois += $log_lambda - log($i);
    $sum_of_poi += exp($log_dpois);
    last if (($i > $lambda) && ($sum_of_poi <= $prev_value));
    $prev_value = $sum_of_poi;
  }
  return($sum_of_poi);
}


##############################################################
### Hypergeometric formula.
###
### Usage:
###     $proba = &RSAT::stats::hypergeometric($m, $n, $k, $x)
### where
###   $m = number of black balls in the urn
###   $n = total number of balls in the urn
###   $k = number of balls in the sample
###   $x = number of black balls in the sample
###
### The hypergeometric probability is
###
###             x  k-x
###            Cm Cn-m
###  P(X=x) = ---------
###               k
###              Cn
###
### For efficiency reasons, the hypergeometric is calculated with the
### recursive formula (Note that the formula is recursive, but the
### computation is iterative).
###
###                       (m - x + 1) (k - x + 1)
###  P(X = x) = P(X=x-1) ------------------------
###                         x ( n - m - k + x)
###
###
### Instead of a product, calculations are performmed with sums of
### logarithms. Then, we make a custom exp conversion to overcome the
### limitation of the Perl exp(x) function.
###
### This gives a precision of the order of e-300.
###
### This routine also calculates the sum of hypergeometrics, when an
### argument 'to' is added.
###   $p = &RSAT::stats::hypergeometric($m, $n, $k, $from, to=>$to);
###
sub hypergeometric {
    my($m, $n, $k, $x, %args) = @_;

    my $w = $n - $m; ### white balls in the urn

    ### error if too high recursion depth
    my $max_depth = 10;
    $args{depth} = $args{depth} || 0;
    if ($args{depth} > $max_depth) {
	&RSAT::error::FatalError("Hypergeometric, recursion reached depth limit $max_depth\n");
    }

    ## initialization
    my $to;
    my $proba = 0;
    my $log_proba = 0;

    ## sum of hypergeometrics
    if (defined($args{to})) {
	$to = $args{to};
    } else {
	$to = $x;
    }

    ## some verbosity
     warn join ( "\t",
  		"Hypergeometric",
  		"m=$m",
  		"w=$w",
  		"n=$n",
  		"k=$k",
  		"x=$x",
 		"depth=$args{depth}",
 		"to=".$to,
 		"check=$args{check}",
 		"prev=$args{previous_value}"
 	      ), "\n" if ($main::verbose >= 6);

    ## incompatible parameter values
    if ($k > $n) {
	&RSAT::error::FatalError("Sample ($k) cannot be larger than number of balls in the urn ($n)");
    }
    if ($m > $n) {
	&RSAT::error::FatalError("Number of black balls in the urn ($m) cannot be larger than number of balls in the urn ($n)");
    }
    if ($w > $n) {
	&RSAT::error::FatalError("Number of white balls in the urn ($m) cannot be larger than number of balls in the urn ($n)");
    }
    if ($x > $k) {
	&RSAT::error::FatalError("Number of black balls in the sample ($x) cannot be larger than sample size ($k)");
    }
    if ($x < 0) {
	&RSAT::error::FatalError( "Number of black balls in the sample ($x) must be strictly positive\n");
    }

    ## In the case of impossible events, a probability of 0 is
    ## returned, unless the routine is called with an argument
    ## &hypergeometric(..., check=>1)
    if ($x > $m) {
	if ($args{check}) {
	    &RSAT::error::FatalError("Number of black balls in the sample ($x) cannot be larger than number of black balls in the urn ($m)");
	} else {
	    return(0);
	}
    }
    if ($k - $x > $w) {
	if ($args{check}) {
	    &RSAT::error::FatalError(join( "", "Number of white balls in the sample (",
					   $k-$x,
					   ") cannot be larger than number of white balls in the urn (",
					   $w,
					   ")"));
	} else {
	    return(0);
	}
    }

#    warn "HERE WE START\n";
    ################################################################
    ## Initialization: compute the first non-null value
    if (defined($args{previous_value})) {
	#### a single recursion is sufficient
#	die;
      my $prev=$args{previous_value};
      $log_proba = log($args{previous_value});
      warn "Recursion from previous value\t", $prev, "log=".$log_proba,  "\n" if ($main::verbose >= 6);
      $start = $x;

    } elsif ($w < $k) {
      ## If the number of non-marked balls ($w) is lower than the size
      ## of the selection ($k), there will be at least $k - $m marked
      ## balls in the selection ($x >= $min_x = $k - $m) -> the proba
      ## of all values from 0 to ($min_x - 1) is null.

#       if ($m >= $k) {
# 	## If the number of marked balls ($m) is higher than or equal to the size of the selection ($k)
# 	##  it is more efficient to temporarily invert the white and black balls
# 	#	  &RSAT::error::FatalError($k, $x, $k-$x);
# 	&RSAT::message::Debug ("inverting marked and non-marked elements for computational efficiency") if ($main::verbose >= 5);
# 	my $proba = &hypergeometric($w, $n , $k, $k - $x, depth=>$args{depth}+1);
# 	if ($proba == 0) {
# 	  return (0);
# 	} else {
# 	  $log_proba = log($proba);
# 	  $start = $x+1;
# 	}

#       } else {

	$min_x = $k - $w;
	$log_proba = 0;
	# 	    for my $i (($min_x + 1)..$k) {
	# 		$log_proba += log($i)
	# 	    }
	for my $i (($n - $k + 1)..$n) {
	  $log_proba -= log($i)
	}
	# 	    for my $i (($w - $k + $min_x + 1)..$w) {
	# 		$log_proba += log($i)
	# 	    }
	for my $i (1..$min_x) {
	  $log_proba -= log($i)
	}
	for my $i (($m - $min_x + 1)..$m) {
	  $log_proba += log($i)
	}
	for my $i (($k - $min_x + 1)..$k) {
	  $log_proba += log($i)
	}
	for my $i (($w -$k + $min_x + 1)..$w) {
	  $log_proba += log($i)
	}
	$start = $min_x + 1;
#      }

    } else {
      ## calculate value for 0 successes
      $log_proba = 0;
      for my $i (($w - $k + 1)..($n - $m)) {
	$log_proba += log($i);
      }
      for my $i (($n - $k + 1)..$n) {
	$log_proba -= log($i);
      }
      $start = 1;
    }
    $proba = &LogToEng($log_proba);
    #    warn "HERE WE PASS\n";
    warn join ( "\t",
		"\tInit value",
		"depth=".$args{depth},
		"Start=".$start,
		"x=".$x,
		# 		    "to=".$to,
		"p=".$proba,
		"CDF=".$proba,
		"log(p)=".$log_proba
	      ), "\n" if ($main::verbose >= 6);

    ################################################################
    ## recursive calculation of the hypergeometric density for $x
    ## (if cumulative, $x is the first value of the sum)
    if ($start <= $x) {
      for my $i ($start..$x) {
	$log_proba += log($m - $i + 1);
	$log_proba -= log($i);
	$log_proba += log($k - $i + 1);
	$log_proba -= log($w - $k + $i);
      }
      $proba = &LogToEng($log_proba);
      warn join ( "\t",
		  "\tFrom value",
		  "depth=".$args{depth},
		  "start=".$start,
		  "x=".$x,
		  # 		    "to=".$to,
		  "p=".$proba,
		  "CDF=".$proba,
		  "log(p)=".$log_proba
		), "\n" if ($main::verbose >= 6);
    }
    #    warn "HERE WE GO\n";

    ################################################################
    ## If $to > $x (cumulative density has to be computed),
    ## pursue the recursive computation and add up the values to the sum
    for my $i (($x + 1)..$to) {
      my $last_proba = $proba; ## for the stopping condition: stop computing when adding up does not change the result anymore
      $log_proba += log($m - $i + 1);
      $log_proba -= log($i);
      $log_proba += log($k - $i + 1);
      $log_proba -= log($w - $k + $i);
      my $new_proba = &LogToEng($log_proba);
      $proba = $last_proba + $new_proba;
      warn join ( "\t",
		  "\tInterm value",
		  "depth=".$args{depth},
		  "i=".$i,
		  # 		    "to=".$to,
		  "last_CFD=".$last_proba,
		  "p=".$new_proba,
		  "CDF=".$proba,
#		  "log(p)=".$log_proba,
		), "\n" if ($main::verbose >= 6);
      last if (($proba == $last_proba) && ($proba >0)); ### beyond precision limit, it is worthless adding more elements
    }
    $proba = &min($proba, 1); ### floating point calculation errors

#     warn join ( "\t",
# 		"Hypergeometric",
# 		"m=$m",
# 		"w=$w",
# 		"n=$n",
# 		"k=$k",
# 		"x=$x",
# 		"to=$to",
# 		"depth=$args{depth}",
# 		"log(p)=$log_proba",
# 		"proba=$proba"
# 	      ), "\n" if ($main::verbose >= 6);

#    my $proba = exp($log_proba);
    return($proba);
}


##############################################################
## sum_of_hypergeometrics($m, $n, $k, $from, $to)
##
## Usage:
##   my $proba = &RSAT::stats::sum_of_hypergeometrics($m, $n, $k, $from, $to);
sub sum_of_hypergeometrics {
    my ($m, $n, $k, $from, $to) = @_;

    ## minimum number of ùmarked balls in the selection
    if ($n - $m - $k < 0) {
	$min_x = $k + $m - $n;
    } else {
	$min_x = $from;
    }
    $from = &max($from, $min_x);

    ## maximum number of marked balls in the selection
    $to = &min($to,$m);

    warn join( "\t", "sum_of_hypergeometrics", "m=$m", "n=$n", "k=$k", "from=$from", "to=$to", "min_x=$min_x"), "\n" if ($main::verbose >= 6);

    if ($min_x > $to) {
	return(0);
    } else {
	return (&hypergeometric($m,$n,$k,$from,to=>$to));
    }
}


##############################################################
### usage: &factorial($n)
sub factorial {
    my $n = $_[0];
    my $fact_n = 1;

    if ($n < 0) {
	print "	Error: invalid entry for factorial calculation\n";
	exit;
    }
    for my $j (1..$n) {
	$fact_n *= $j;
    }
    return $fact_n;
}

##############################################################
## calculates the probability according to Fisher's exact test
## Usage
## =====
## $proba = &FisherExactTest($row_nb, $col_nb, @values);
##
## where
##	$row_nb is the number of rows
##	$col_nb is th number of columns
##	@values is the list of values
## all values must be natural numbers (0,1,2, ...)
## the number of values must equal the product of col_nb by row_nb
##
## the first step is to calculate the marginal sums:
## ni+ = sum of all values from the ith row
## n+j = sum of all values from the jth column
## N = sum of all values from the table
##
## The probability is then calculated by:
##
##         PROD(ni+!)PROD(n+j!)
## proba = --------------------
##          N!PROD(PROD(nij!))
##
## the input data are reported together with all marginal sums
## by setting a global variable called $report_data to 1
sub FisherExactTest {
    my ($row_nb, $col_nb, @values) = @_;
    my $N;
    my $row;
    my $col;
    my @col_sum;
    my @row_sum;
    my $offset;
    my $i;
    my $proba;
    my $log_proba;
    my $log10_proba;
    my $val_nb = $#values + 1;
    my $left_group = 0;
    my $left_sum = 0;
    my $right_group = 0;
    my $right_sum = 0;

    ## Check parameters
    if (($row_nb < 2) || ($col_nb < 2)) {
	&RSAT::error::FatalError( ";FisherExactTest: at least 2 rows and 2 columns are required for the Fisher exact test"); # too few rows or columns
    }
    unless ($val_nb == $row_nb * $col_nb) {
	&RSAT::error::FatalError( ";FisherExactTest: invalid number of numeric values",
				  "\t".$row_nb." rows",
				  "\t".$col_nb." columns",
				  "\t".$val_nb." values");
    }
    foreach $v (@values) {
	return ";Error: $v is not a natural number" unless &IsNatural($v); # invalid values
    }


    ## calculate marginal sums ni+, n+j, and N
    $N = 0;
    for $row (1..$row_nb) {
	$col_sum[$row] = 0;
	for $col (1..$col_nb) {
	    $offset = ($row-1) * $col_nb + $col -1;
	    $col_sum[$row] += $values[$offset];
	}
	$N += $col_sum[$row];
    }
    for $col (1..$col_nb) {
	$row_sum[$col] = 0;
	for $row (1..$row_nb) {
	    $offset = ($row-1) * $col_nb + $col -1;
	    $row_sum[$col] += $values[$offset];
	}
    }

    ## calculate the probability
    $log_proba = 0;
    for $row (1..$row_nb) {

	for $i (2..$col_sum[$row]) {
	    $log_proba += log($i);
	}
    }
    for $col (1..$col_nb) {
	for $i (2..$row_sum[$col]) {
	    $log_proba += log($i);
	}
    }
    for $i (2..$N) {
	$log_proba -= log($i);
    }
    for $row (1..$row_nb) {
	for $col (1..$col_nb) {
	    $offset = ($row-1) * $col_nb + $col -1;
	    for $i (2..$values[$offset]) {
		$log_proba -= log($i);
	    }
	}
    }
    $proba = exp($log_proba);
    $log10_proba = $log_proba / log(10);

    ### data report
    if ($report_data) {
	print $out ";Fisher's exact test with $row_nb rows and $col_nb col ($val_nb values)\n";
	print $out ";DATA REPORT\n";
	print $out ";";
	for $col (1..$col_nb) {
	    print $out "\tA$col";
	}
	print $out "\tni+\n";
	for $row (1..$row_nb) {
	    print $out ";B$row";
	    for $col (1..$col_nb) {
		$offset = ($row-1) * $col_nb + $col -1;
		print $out "\t$values[$offset]";
	    }
	    print $out "\t$col_sum[$row]\n";
	}
	print $out ";n+j";
	for $col (1..$col_nb) {
	    print $out "\t$row_sum[$col]";
	}
	print $out "\t$N\n";
    }
    return $proba;
} ### end FisherExactTest



=pod

=item B<ChiSquare>

Calculate Pearson chi-square statistics for a table of numbers.

The input is an array containing the values of a data table on which
the chi2 test will be applied. Each row of the input table is
considered as a series of values. The number of columns corresponds to
the number of values per series ("classes").

For the goodness of fit test, there must be exactly two series: the
first series indicates observed values, the second expected
values. The homogeneity and independence test can be applied to tables
containing two or more series.

=head2 Usage

 $chi_square = &ChiSquare($test,$series_nb, $class_nb, $sort_by_exp, @values);

where

=over

=item I<$test> indicates the kind of  of hypothesis to test:

Supported: independence, homogeneity, goodness (of fit)

=item I<$series_nb>: number of series.

=item I<$class_nb>: number of classes

=item I<$sort_by_exp>: sort series by expected values (if parameter is > 0)

Only for goodness of fit test. Sort obs and exp series by increasing
expected values before grouping tail classes. This reduces the number
of classes to be grouped, since after sorting all the classes with exp
< 5 are on the same side of the sorted arrays.

=item I<@values>: list of values

All values must be real numbers.

The number of values must equal the product of series_nb by class_nb.

=head2 Goodness of  fit test

For the goodness of fit test, there must be exactly two series: the
first series containsq the observed frequencies, the second series the
expected frequencies.

The chi-square value is calculated by:

              (obs_j - exp_j)^2
 ChiSq = SUM  ----------------
          j      (exp_j)

where j is the index for classes.

=head2 Independence or homogeneity test

The first step is to calculate the marginal sums:

 I<ni.> = sum of all values from the ith row
 I<n.j> = sum of all values from the jth column
 I<N> = sum of all values from the table

The chi-square value is then calculated by:

                   (nij - n.j*ni./N)^2
 ChiSq = SUM (SUM -----------------)
          j    i     (n.j*ni./N)

Input data can be reported together with all marginal sums by setting
a global variable called $report_data to 1.


=head2  Applicability condition

One condition of applicability for the chi-square test is that each
class should contain a "sufficient" number of expected observations.
One commonly takes 5 as the minimal number of expected observations
per class.

If any of the expected values is 0, the chi2 value is Inf.

When the condition of acceptability is not met, the &ChiSquare() function
returns the calculated value surrounded by curly brackets, in order to
warn the user that the chi2 value is not valid.

=cut

sub ChiSquare {
    my($test, $series_nb, $class_nb,  $sort_by_exp, @values) = @_;
    my $N;
    my $row;
    my $col;
    my @col_sum;
    my @row_sum;
    my $offset;
    my $i;
    my $chi_square;
    my $val_nb = $#values + 1;
    my $obs;
    my $exp;
    my $not_valid = 0;
    my $group_tails = 1;
    my $left_group = 0;
    my $left_sum = 0;
    my $right_group = 0;
    my $right_sum = 0;
    my @observed = ();
    my @expected = ();

    &RSAT::message::Debug("RSAT::stat::ChiSquare", "test=".$test, "row_nb=".$series_nb, "col_nb=".$class_nb, scalar(@values)." values")
	if ($main::verbose >= 4);

    ### Check parameter values
    unless (($test =~ /^indep/i) || ($test =~ /^good/i) || ($test =~ /^homog/i)) {
	&RSAT::error::FatalError("Unknown test (supported: 'goodness of fit', 'homogeneity' and 'independence')");
    }
    if (($series_nb < 2) || ($class_nb < 2)) {
	&RSAT::error::FatalError("At least 2 rows and 2 columns are required"); # too few rows or columns
    }

    ## values = columns*rows
    unless ($val_nb == $series_nb * $class_nb) {
	&RSAT::error::FatalError( ";ChiSquare: invalid number of numeric values",
				  "\t".$series_nb." rows",
				  "\t".$class_nb." columns",
				  "\t".$val_nb." values");
    }

    ### All values must be real numbers
    foreach $val (@values) {
	&RSAT::error::FatalError("Invalid values: $val is not a real number") unless &RSAT::util::IsReal($val); # invalid values
    }

    ################################################################
    ## Independence of homogeneity test
    if (($test =~ /^indep/i) || ($test =~ /^homog/i)) {

	#### Calculate marginal sums ni+, n+j, and N
	$N = 0;
	for $row (1..$series_nb) {
	    $col_sum[$row] = 0;
	    for $col (1..$class_nb) {
		$offset = ($row-1) * $class_nb + $col -1;
		$col_sum[$row] += $values[$offset];
	    }
	    $N += $col_sum[$row];
	}
	for $col (1..$class_nb) {
	    $row_sum[$col] = 0;
	    for $row (1..$series_nb) {
		$offset = ($row-1) * $class_nb + $col -1;
		$row_sum[$col] += $values[$offset];
	    }
	}

	#### Calculate the chi square value
	$chi_square = 0;
	for $row (1..$series_nb) {
	    for $col (1..$class_nb) {
		$exp = $row_sum[$col] * $col_sum[$row] / $N;
		$not_valid = 1 if ($exp < 5);
		$offset = ($row-1) * $class_nb + $col -1;
		$obs = $values[$offset];
		if ($exp == 0) {
		    $chi_square = "NA";
		    &RSAT::message::Warning("&ChiSquare()", "chi2 is undefined (NA) because some expected values are null (some rows or columns contain only null values).");
		    last;
		} else {
		    $chi_square += (($obs - $exp)**2)/$exp unless ($exp == 0);
		}
	    }
	    last if ($chi_square eq "NA");
	}

	################################################################
	## Goodness of fit test
    } elsif ($test =~ /^good/i) {
	unless ($series_nb == 2) {
	    &RSAT::error::FatalError("Goodness of fit test requires esactly two lines");
	}
	$chi_square = 0;
	@observed = @values[0..$class_nb-1];
	@expected = @values[$class_nb..2*$class_nb-1];

  	if ($main::verbose >= 10) {
	    &RSAT::message::Debug("rank original", join("\t", 1..scalar(@expected)));
  	    &RSAT::message::Debug("obs original", join( "\t", @observed));
  	    &RSAT::message::Debug("exp original", join( "\t", @expected));
  	}

	#### Group tails to fulfill the condition of applicability
	if ($group_tails) {
	    &RSAT::message::Info("&ChiSquare()", "Grouping classes on left and right tails to ensure exp >= 5.") if ($main::verbose >= 3);

	    ## Sort by increasing expected values before grouping
	    ## terminal classes.  This is , in order to get the two
	    ## tails on the same side, which generally reduces the
	    ## number of classes to be grouped. I am not sure this
	    ## treatment is orthodox (I did not find it in any
	    ## textbook) but since the Chi2 test does not make any
	    ## assumption on the order of the columns, I don't see any
	    ## reason to forbid it.
	    if ($sort_by_exp) {
		&RSAT::message::Info("&ChiSquare()", "Sorting by expected value before class grouping") if ($main::verbose >= 3);
		@order = sort {$expected[$a] <=> $expected[$b]} 0..$#expected;
		@observed_sorted = @observed[@order];
		@expected_sorted = @expected[@order];
		if ($main::verbose >= 10) {
		    &RSAT::message::Debug("rank", join("\t", 1..scalar(@expected)));
		    &RSAT::message::Debug("order by incr exp", join( "\t", @order));
		    &RSAT::message::Debug("obs sorted by exp", join( "\t", @observed_sorted));
		    &RSAT::message::Debug("exp sorted by exp", join( "\t", @expected_sorted));
		}
		@observed = @observed_sorted;
		@expected = @expected_sorted;
	    }

	    ## Report expected and observed vectors before grouping
	    if ($main::verbose >= 10) {
		&RSAT::message::Debug("obs before grouping", join("\t", @observed));
		&RSAT::message::Debug("exp before grouping", join("\t", @expected));
	    }

	    ## Group left tail
	    $left_group = 0; ## Counter for the number of grouped classes on the left tail
	    while (($#expected > 1) && ($expected[0] < 5)) {
		$left_group++;
		$exp_sum = shift(@expected);
		$expected[0] += $exp_sum;
		$obs_sum = shift(@observed);
		$observed[0] += $obs_sum;
	    }

	    ## Recalculate the number of columns after left grouping
	    $class_nb = scalar(@expected);
	    &RSAT::message::Debug("Left tail collapsed", $left_group, "remaining", $class_nb) if ($main::verbose >= 3);
	    if ($main::verbose >= 10) {
		&RSAT::message::Debug("rank after left grouping", join("\t", 1..scalar(@expected)));
		&RSAT::message::Debug("obs after left grouping", join("\t", @observed));
		&RSAT::message::Debug("exp after left grouping", join("\t", @expected));
	    }

	    unless ($sort_by_exp) {
		## Group right tail
		$right_group = 0; ## Counter for the number of grouped classes on the right tail
		while (($#expected > 1) && ($expected[$#expected] < 5)) {
		    $right_group++;
		    $exp_sum = pop(@expected);
		    $expected[$#expected] += $exp_sum;
		    $obs_sum = pop(@observed);
		    $observed[$#observed] += $obs_sum;
		}

		## Recalculate the number of columns after right grouping
		$class_nb = scalar(@expected);
		&RSAT::message::Debug("Right tail collapsed", $right_group, "remaining", $class_nb) if ($main::verbose >= 3);
		if ($main::verbose >= 10) {
		    &RSAT::message::Debug("rank after right grouping", join("\t", 1..scalar(@expected)));
		    &RSAT::message::Debug("obs after right grouping", join("\t", @observed));
		    &RSAT::message::Debug("exp after right grouping", join("\t", @expected));
		}
	    }
 	}

	## Calculate the observed chi-square
	for $col (1..$class_nb) {
	    $obs = $observed[$col-1];
	    $exp = $expected[$col-1];
	    $not_valid = 1 if ($exp < 5);
	    if ($exp == 0) {
		$chi_square = "ERR_exp0";
		&RSAT::message::Warning("Cannot calculate the observed chi2 because class $col has an expected value of 0");
		last;
	    } else {
		$chi_square +=  (($obs - $exp)**2)/$exp;
	    }
	    &RSAT::message::Info($obs, $exp, $not_valid, $chi_square) if ($main::verbose >= 5);
	}
    }

    ################################################################
    ## Calculate the degrees of freedom
    my $deg_freedom = ($series_nb - 1)*($class_nb - 1);

    ################################################################
    ### data report
    if ($report_data) {
	print $out ";Pearson chi-square statistics with $series_nb rows and $class_nb col ($val_nb values)\n";
	print $out ";Test of $test\n";
	print $out ";$deg_freedom degrees of freedom\n";
	print $out ";DATA REPORT\n";
	if (($not_valid) && ($main::verbose >= 2)) {
	    print $out "; WARNING: data do not conform to the applicability conditions \n; (some classes have less than 5 expected observations).\n";
	}
	print $out ";";
	for $col (1..$class_nb) {
	    print $out "\tA$col";
	}
	if (($test =~ /^indep/i) || ($test =~ /^homog/i)) {
	    print $out "\tni+\n";
	    for $row (1..$series_nb) {
		print $out ";B$row";
		for $col (1..$class_nb) {
		    $offset = ($row-1) * $class_nb + $col -1;
		    print $out "\t$values[$offset]";
		}
		print $out "\t$col_sum[$row]\n";
	    }
	    print $out ";n+j";
	    for $col (1..$class_nb) {
		print $out "\t$row_sum[$col]";
	    }
	    print $out "\t$N\n";
	} else {
	    print $out "\n;obs";
	    for $col (1..$class_nb) {
		print $out "\t$values[$col-1]";
	    }
	    print $out "\n;exp";
	    for $col (1..$class_nb) {
		print $out "\t$values[$class_nb + $col-1]";
	    }
	    print $out "\n";
	}
    }

    my $answer = $chi_square;
    if (!&RSAT::util::IsReal($chi_square)) {
	$answer = $chi_square;
#    } elsif ($not_valid) {
#      $answer = sprintf "{%.3f}", $chi_square;
    } else {
	$answer = sprintf "%.5f", $chi_square;
    }
    &RSAT::message::Debug($answer,
			  $deg_freedom,
			  $left_group,
			  $right_group,
			  $#expected,
			  $expected[0],
			  $expected[$#expected],
	) if ($main::verbose >= 6);

    &RSAT::message::Debug("chi2=",$answer, "df=".$deg_freedom, "left_group=".$left_group, "right_group=".$right_group, "obs:".scalar(@observed), "exp:".scalar(@expected))
	if ($main::verbose >= 5);

    return ($answer, $deg_freedom, $left_group, $right_group, \@observed, \@expected);
} ### end ChiSquare


##############################################################
##
## calculates log-likelihood statistics for a table of numbers
##
## Usage
## =====
## $log_likelihood = &LogLikelihood($row_nb, $col_nb, @values);
##
## where
##	$row_nb is the number of rows
##	$col_nb is th number of columns
##	@values is the list of values
## all values must be real numbers
## the number of values must equal the product of col_nb by row_nb
##
## the first step is to calculate the marginal sums:
## ni+ = sum of all values from the ith row
## n+j = sum of all values from the jth column
## N = sum of all values from the table
##
## The log-likelihood value is then calculated by:
##
##
## LogLikelihood = N*log(N) + SUM SUM nij*log(nij)
##                             i   j 
##
##                          - SUM ni+*log(ni+) - SUM n+j*log(n+j)
##                             i                  j
##
##
## the input data are reported together with all marginal sums
## by setting a global variable called $report_data to 1
sub LogLikelihood {
    my ($row_nb, $col_nb, @values) = @_;
    my $N;
    my $row;
    my $col;
    my @col_sum;
    my @row_sum;
    my $offset;
    my $i;
    my $log_likelihood;
    my $val_nb = $#values + 1;
    my $obs;
    my $deg_freedom = ($row_nb - 1)*($col_nb - 1);

    if (($row_nb < 2) || ($col_nb < 2)) {
	return ";Error: at least 2 rows and 2 columns are required"; # too few rows or columns
    }
    unless ($val_nb == $row_nb * $col_nb) {
	return ";Error: invalid number of numeric values"; # inappropriate number of values provided
    }
    foreach $v (@values) {
	return ";Error: $v is not a real number" unless &RSAT::util::IsReal($v); # invalid values
    }

#    foreach $v (@values) {
#	return ";Error: ($v) all values must be strictly positive" unless ($v > 0); # invalid values
#    }

    ## calculate marginal sums ni+, n+j, and N
    $N = 0;
    for $row (1..$row_nb) {
	$col_sum[$row] = 0;
	for $col (1..$col_nb) {
	    $offset = ($row-1) * $col_nb + $col -1;
	    $col_sum[$row] += $values[$offset];
	}
	$N += $col_sum[$row];
    }
    for $col (1..$col_nb) {
	$row_sum[$col] = 0;
	for $row (1..$row_nb) {
	    $offset = ($row-1) * $col_nb + $col -1;
	    $row_sum[$col] += $values[$offset];
	}
    }

    ## calculate the log likelihood value
    $log_likelihood = - $N * log($N);
    for $row (1..$row_nb) {
        for $col (1..$col_nb) {
	    $offset = ($row-1) * $col_nb + $col -1;
	    $obs = $values[$offset];
	    $log_likelihood -= $obs * log($obs) unless ($obs <= 0);
        }
    }
    for $row (1..$row_nb) {
	$log_likelihood += $col_sum[$row] * log($col_sum[$row]) unless ($col_sum[$row] <= 0);
    }
    for $col (1..$col_nb) {
	$log_likelihood += $row_sum[$col] * log($row_sum[$col]) unless ($row_sum[$col] <= 0);;
    }
    $log_likelihood = 0 if ($log_likelihood > 0);

    ### data report
    if ($report_data) {
	print $out ";Log-likelihood statistics with $row_nb rows and $col_nb col ($val_nb values)\n";
	print $out ";DATA REPORT\n";
	print $out ";";
	for $col (1..$col_nb) {
	    print $out "\tA$col";
	}
	print $out "\tni+\n";
	for $row (1..$row_nb) {
	    print $out ";B$row";
	    for $col (1..$col_nb) {
		$offset = ($row-1) * $col_nb + $col -1;
		print $out "\t$values[$offset]";
	    }
	    print $out "\t$col_sum[$row]\n";
	}
	print $out ";n+j";
	for $col (1..$col_nb) {
	    print $out "\t$row_sum[$col]";
	}
	print $out "\t$N\n";
	print $out ";$deg_freedom degrees of freedom\n";
	print $out ";-2log(L) = ",-2*$log_likelihood,"\n";
    }
    return $log_likelihood;
} ### end LogLikelihood


1;
