################################################################
#### This package is calculating the Area Under the Curve from ROC curve
#### informations.
####
#### This has been implemented by Rekin's Janky taking inspiration from the
#### routine implemented by Alejandra Marina Rivera
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
    my $AUCg_local=&calc_geometric_in_local_area(@_);	
    my $AUCg_total=&calc_geometric_in_total_area(@_);	
    my $AUCg_total_extended=&calc_geometric_in_total_area_with_extremes(@_);	
    return($AUCg_local,$AUCg_total,$AUCg_total_extended);
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

# sub calc_MannWhitney{
#     my ($Sn, $FPR) = @_;


#     my ($base, # distance between 2 points in x axis
# 	$x1, $x2, # x points
# 	$y1, $y2, # y points
# 	$height_t, # height of the triangle
# 	$areaR, # area of rectangle
# 	$areaT, # area of triangle
# 	);

#     my $AUC=0 ;

#     #Before that add the points (TP,FP) =(0,0) and (TP,FP) = (1,1)
#     unshift(@$FPR,0);  
#     unshift(@$Sn,0);
#     push(@$FPR,1);  
#     push(@$Sn,1);
    
#     #Use the ranks to calculate the Mann-Whitney statistic 
# 	#U = Npos*Nneg + Nneg(Nneg+1)/2 - R, 
    
# for my $i (0..$#{$FPR}-1) { # foreach element of the x array (FPR)
# 	my $U = $Npos*$Nneg + $Nneg($Nneg+1)/2 - $R, 
# 	$x1=$FPR->[$i]; 
# 	$x2=$FPR->[$i+1];
# 	$base= $x2 - $x1;
# 	$y1=$Sn->[$i];
# 	$y2=$Sn->[$i+1];
# 	$height_t=abs($y2-$y1);
# 	$areaR= $base * $y2;
# 	$areaT= ($base * $height_t) / 2 ;
# 	$AUC += ($areaR + $areaT);
# #	&OPTOOLS::util::Warning(1,join("\t", "x1",$x1,"x2",$x2,"AUC",$AUC));
#     }
#     return ($AUC);

# # https://list.scms.waikato.ac.nz/pipermail/wekalist/2004-January/002113.html
# #
# # 1. extract TP and FP rates (forming the ROC curve), keep them in threshold
# # (= P) order

# # 2. Use the nonparametric AUC (Area Under the Curve) approximation with the
# # Mann-Whitney rank sum test as follows. Before that add the points (TP,FP) =
# # (0,0) and (TP,FP) = (1,1) of the
# # ROC (Weka leaves them out of the curve) and link them to the most extreme
# # threshold values (mostly 0 and 1 if it
# # scales between 0 and 1).

# # a. The threshold order (0 to 1) determines the rank (assign an average rankorder
# # over obeservations with the same threshold to correct for ties)

# # b. Use the ranks to calculate the Mann-Whitney statistic 
# # U = Npos*Nneg + Nneg(Nneg+1)/2 - R, 
# # with 
# #         R = the ranksum of the negative sample,
# #         Npos= number of positives in the sample and
# #         Nneg= the number of negatives in the sample

# # c. AUC = U/(Npos + Nneg)
# #    ---
# # 3. The standard error for AUC can be approximated Using the Hanley MacNeil
# # test following:
# #        --------------
# #         /* For clarity, the same symbols are used as in the Hanley-McNeil
# # paper. */
# #         n_A    = Npos;
# #         n_N    = Nneg;
# #         theta  = AUC;
# #         //
# #         theta2 = theta * theta;
# #         Q1     = theta / (2 - theta);
# #         Q2     = 2 * theta2 / (1 + theta);
# #         se2    = (theta * (1 - theta) + (n_A - 1) * (Q1 - theta2) + (n_N -
# # 1) * (Q2
# #         - theta2)) / (n_A * n_N);
# #         //
# #         SE_auc = squareRoot (se2);
# #         Assumption:
# #         0.0 > AUC < 1.0 , Npos > 0 and Npos > 0
	
# # The threshold values for every observation can also be obtained from the
# # classification function if made external (like for weka logistic regression), 
# # and the TP and FP fraction from the cummulative frequency for negatives and 
# # positives when sorted by the classification function values.
# }

1;
