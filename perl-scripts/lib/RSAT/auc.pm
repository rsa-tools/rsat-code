################################################################
#### This package is calculating the Area Under the Curve from ROC curve
#### informations.
####
#### This has been implemented by Rekin's Janky taking inspiration from the
#### routine impemented by Alejandra Marina Rivera
################################################################

package RSAT::auc;



use RSAT::GenericObject;
use RSAT::message;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );



################################################################
#### Calculate geometrically the Area Under the Curves (AUC)
####
#### The total area is calculated by summing the area for each adjacent points
#### of the curve, (x1,y1) and (x2,y2), and the base. This is done by
#### calculating the area of the rectangle between (x1,0),(x2,0),(x1,y1) and
#### (x2,y1), then the area of the triangle corresponding to (x1,y1),(x2,y1)
#### and (x2,y2).
#### 
#### This method is also called the trapezoïdal method
################################################################

sub calc_geometric{
    &RSAT::message::Warning(join("\t","Calculating the AUCg (using geometric method)")) if ($main::verbose >= 2);
    my $AUC=&calc_geometric_in_total_area(@_);	
    return($AUC);
}

sub calc_geometric_all{
    &RSAT::message::Warning(join("\t","Calculating the AUCg (using geometric method)")) if ($main::verbose >= 2);
#    my $AUCg_local=&calc_geometric_in_local_area(@_);	
    my $AUCg_total=&calc_geometric_in_total_area(@_);	
    my $AUCg_total_extended=&calc_geometric_in_total_area_with_extremes(@_);	
#    return($AUCg_total,$AUCg_total_extended,$AUCg_local);
    return($AUCg_total,$AUCg_total_extended);
}

sub calc_geometric_in_total_area{
    my ($Sn, $FPR) = @_;
    &RSAT::message::Warning(join("\t","Calculating the AUCgt : AUCg expressed as a ratio of the total area (the total area is estimated as the area between (0,0) and (1,1) )")) if ($main::verbose >= 2);
    my $AUC=0 ;	
    my $total_area = 1;

    for my $i (0..$#{$FPR}-1) { # foreach element of the x array (FPR)
	$x1=$FPR->[$i]; 
	$x2=$FPR->[$i+1];
	$y1=$Sn->[$i];
	$y2=$Sn->[$i+1];
	$AUC += &area($x1,$x2,$y1,$y2);
    }

    $AUC = $AUC/$total_area;
    return ($AUC);
}

sub calc_geometric_in_total_area_with_extremes{
    my ($Sn, $FPR) = @_;
    &RSAT::message::Warning(join("\t","Calculating the AUCgtx : AUCg expressed as a ratio of the total area (the total area is estimated as the area between (0,0) and (1,1), those points are added if not existing)")) if ($main::verbose >= 2);
    my $AUC=0 ;	
    my $total_area = 1;

    unless (($FPR->[0]==0)&&($Sn->[0]==0)){
	&RSAT::message::Warning(join("\t","\tAdding the point (FPR,TPR)=(0,0)")) if ($main::verbose >= 2);
	$AUC += &area(O,$FPR->[0],0,$Sn->[0]);
    }

    for my $i (0..$#{$FPR}-1) { # foreach element of the x array (FPR)
	$x1=$FPR->[$i]; 
	$x2=$FPR->[$i+1];
	$y1=$Sn->[$i];
	$y2=$Sn->[$i+1];
	$AUC += &area($x1,$x2,$y1,$y2);
    }

    # Add (TP,FP) = (1,1)
    unless (($FPR->[$#{$FPR}-1]==1)&&($Sn->[$#{$FPR}-1]==1)){
	&RSAT::message::Warning(join("\t","\tAdding the point (FPR,TPR)=(1,1)")) if ($main::verbose >= 2);
	$AUC += &area($FPR->[$#{$FPR}],1,$Sn->[$#{$FPR}],1);
    }

    $AUC = $AUC/$total_area;
    return ($AUC);
}

sub calc_geometric_in_local_area{
    my ($Sn, $FPR) = @_;
    &RSAT::message::Warning(join("\t","Calculating the AUCgl : AUCg expressed as a ratio of the local area (the total area is estimated as the area between the extremes values (min,max) of the FPR and Sn)")) if ($main::verbose >= 2);
    my $AUC=0 ;	

    my $maxFPR=0;    
    my $minFPR=0;
    my $maxSn=0;    
    my $minSn=0;
    for my $i (0..$#{$FPR}-1) { # foreach element of the x array (FPR)
	
	$x1=$FPR->[$i]; 
	$x2=$FPR->[$i+1];
	$y1=$Sn->[$i];
	$y2=$Sn->[$i+1];
	$AUC += &area($x1,$x2,$y1,$y2);
#	&RSAT::message::Warning(join("\t",$AUC));
	$maxFPR = $x2 if ($x2 > $maxFPR);
	$maxSn = $y2 if ($y2 > $maxSn);
	$minFPR = $x1 if ($x1 < $minFPR);
	$minSn = $y1 if ($y1 < $minSn);
    }

    if  ($main::verbose >= 2){
	&RSAT::message::Warning(join("\t","\tFPR","min",$minFPR,"max",$maxFPR));
	&RSAT::message::Warning(join("\t","\tSn","min",$minSn,"max",$maxSn));
    }

    $total_area = ($maxFPR - $minFPR)*($maxSn - $minSn);
    $AUC = $AUC/$total_area;
    return ($AUC);
}

################################################################
#### Calculate the trapezoidal area given the coordinates
################################################################

sub area{
    my ($x1,$x2,$y1,$y2)=@_;
    my $base= $x2 - $x1;  # distance between 2 points in x axis
    my $height_t=abs($y2-$y1); # height of the triangle
    my $areaR= $base * $y1; # area of rectangle
    my $areaT= ($base * $height_t) / 2 ;# area of triangle
    my $AUC = ($areaR + $areaT);
#    &RSAT::message::Warning(join("\t",$x1,$y1,$x2,$y2,$AUC));
    return ($AUC);
}

################################################################
#### Calculate the Area Under the Curves (AUC) using Normal distribution
#### assumption
################################################################

sub calc_normal{
    my ($Sn, $FPR) = @_;
    &RSAT::message::Warning(join("\t","Calculating the AUCn (using Normal distribution assumption)")) if ($main::verbose >= 2);
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

################################################################
#### Calculate the Area Under the Curves (AUC) using Mann-Whitney U
################################################################
# to do in R :
#    x1 = x[y==1]; n1 = length(x1); 
#    x2 = x[y==0]; n2 = length(x2);
#     r = rank(c(x1,x2))  
#   auc = (sum(r[1:n1]) - n1*(n1+1)/2) / (n1*n2) 

sub calc_MannWhitney{
    my ($Sn, $FPR) = @_;
    my $AUC=0 ;	    

    # Add (TP,FP)=(0,0)
    unless (($FPR->[0]==0)&&($Sn->[0]==0)){
	&RSAT::message::Warning(join("\t","\tAdding the point (FPR,TPR)=(0,0)")) if ($main::verbose >= 2);
	unshift(@$FPR,0.000);  
	unshift(@$Sn,0.000);
     }
    
    # Add (TP,FP) = (1,1)
    unless (($FPR->[$#{$FPR}-1]==1)&&($Sn->[$#{$FPR}-1]==1)){
	&RSAT::message::Warning(join("\t","\tAdding the point (FPR,TPR)=(1,1)")) if ($main::verbose >= 2);
	push(@$FPR,1.000);  
	push(@$Sn,1.000);
    }

    my $n1 = scalar(@{$Sn});  # number of elements in Sn
    my $n2 = scalar(@{$FPR}); # number of elements in FPR

    ################################################################
    #  Use the ranks to calculate the Mann-Whitney statistic 

    my %freq=();
    my %min_rank=();
    my %max_rank=();
    my $r=0;
    my $R = 0; # sum of the ranks of sensitivities

    foreach my $val (sort {$a<=>$b} (((@{$FPR}),(@{$Sn})))) {
	$r++;
	$min_rank{$val}=$r unless ($freq{$val});
	$max_rank{$val}=$r;
	$freq{$val}++;
#	&RSAT::message::Warning(join("\t",$r,"value",$val,$freq{$val},"min",$min_rank{$val},"max",$max_rank{$val})) if ($main::verbose >= 10);
    }
    
    my %rank=();
    foreach my $val (sort keys %freq){
	my $min = $min_rank{$val};
	my $max = $max_rank{$val};
	my $fq = $freq{$val};
	$rank{$val}=$min +(($max-$min)/2);
#	&RSAT::message::Warning(join("\t",$val,"rank",$rank{$val})) if ($main::verbose >= 10);
    }

    for my $i (0..$#{$Sn}) { # foreach element of the Sn array
	$R+= $rank{$Sn->[$i]};
#	&RSAT::message::Warning(join("\t","Sn",$Sn->[$i],"rank",$rank{$Sn->[$i]},"R",$R)) if ($main::verbose >= 5);
    }
    my $U=($R - $n1*($n1+1)/2) / ($n1*$n2);
    &RSAT::message::Warning(join("\t","AUC",$U)) if ($main::verbose >= 2);
    return ($U);
}


1;
