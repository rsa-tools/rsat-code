################################################################
#### This package is calculating the Area Under the Curve from ROC curve
#### informations.
####
#### This has been implemented by Rekin's Janky taking inspiration from the
#### routine impemented by Alejandra Marina Rivera
################################################################

package auc;

# require Exporter;
# require AutoLoader;
# @ISA = qw( Exporter AutoLoader );
# @EXPORT = qw(test);

################################################################
#### Calculate geometrically the Area Under the Curves (AUC)
####
#### The total area is calculated by summing the area for each adjacent points
#### of the curve, (x1,y1) and (x2,y2), and the base. This is done by
#### calculating the area of the rectangle between (x1,0),(x2,0),(x1,y1) and
#### (x2,y1), then the area of the triangle corresponding to (x1,y1),(x2,y1)
#### and (x2,y2).
################################################################

sub calc_geometric{
    my ($Sn, $FPR) = @_;

    my ($base, # distance between 2 points in x axis
	$x1, $x2, # x points
	$y1, $y2, # y points
	$height_t, # height of the triangle
	$areaR, # area of rectangle
	$areaT, # area of triangle
	);

    my $AUC=0 ;	
    for my $i (0..$#{$FPR}-1) { # foreach element of the x array (FPR)
	
	$x1=$FPR->[$i]; 
	$x2=$FPR->[$i+1];
	$base= $x2 - $x1;
	$y1=$Sn->[$i];
	$y2=$Sn->[$i+1];
	$height_t=abs($y2-$y1);
	$areaR= $base * $y2;
	$areaT= ($base * $height_t) / 2 ;
	$AUC += ($areaR + $areaT);
#	&OPTOOLS::util::Warning(1,join("\t", "x1",$x1,"x2",$x2,"AUC",$AUC));
    }
    return ($AUC);
}

################################################################
#### Calculate the Area Under the Curves (AUC) using Normal distribution
#### assumption
################################################################

sub calc_normal{
    my ($Sn, $FPR) = @_;
    my ($X_mean, $X_var) = mean_and_var($FPR);
    my ($Y_mean, $Y_var) = mean_and_var($Sn);
		
    use Math::CDF;
    my $AUC=0 ;
    my $Z=($X_mean - $Y_mean)/sqrt($X_var + $Y_var);
    $AUC=1-&Math::CDF::pnorm($Z);  
    return ($AUC);
}
		

sub mean_and_var {
	my $array = shift;
	
	my $n = $#{$array};

	my $mean = 0;
	$mean += $_ foreach @{$array};
	$mean /= $n;
	
	my $var = 0;
	$var += ($_ - $mean) **2 foreach @{$array};
	$var /= $n;
	
	return ($mean,$var);
}

1;
