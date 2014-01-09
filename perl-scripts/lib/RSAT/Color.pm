###############################################################
#
# A class for util handling
#
package RSAT::Color;
use RSAT::GenericObject;
use RSAT::message;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

=pod

=head1 NAME

    RSAT::Color

=head1 DESCRIPTION

Color specifications for graphics.

=head1 METHODS

=cut



################################################################

=pod

=item RGBColorPalette

Redutns an array of colors, each defined as an (R, G, B) array, where
R, G and B are Integer values between 0 and 255.

Usage: 
  my @palette = &RSAT::Color::RGBColorPalette();
  foreach my $c (@palette) {
    my $color_ref = $palette[$c];
    my ($R,$G,$B, $name) = @{$color_ref};
    print join (",", $R,$G,$B, $name);
  }

=cut
sub RGBColorPalette {
  my @white = qw(255 255 255 white);
  my @black = qw(0 0 0 black);
  my @gray_032 = qw(32 32 32 gray_032);
  my @gray_064 = qw(64 64 64 gray_064);
  my @gray_096 = qw(96 96 96 gray_096);
  my @gray_125 = qw(125 125 125 gray_125);
  my @gray_150 = qw(150 150 150 gray_150);
  my @gray_200 = qw(200 200 200 gray_200);
  my @gray_225 = qw(225 225 225 gray_225);
  my @yellow = qw(255 255 0 yellow);
  my @yellow_225 = qw(225 225 0 yellow_225);
  my @yellow_200 = qw(200 200 0 yellow_200);
  my @yellow_127 = qw(127 127 0 yellow_127);
  my @yellow_light = qw(255 255 150 yellow_light);
  my @yellow_pale = qw(255 255 200 yellow_pale);
  my @green = qw(0 255 0 green);
  my @green_200 = qw(0 200 0 green_200);
  my @green_175 = qw(0 175 0 green_175);
  my @green_127 = qw(0 127 0 green_127);
  my @green_096 = qw(0 96 0 green_096);
  my @green_064 = qw(0 64 0 green_064);
  my @cyan = qw(0 255 255 cyan);
  my @cyan_200 = qw(0 200 200 cyan_200);
  my @cyan_127 = qw(0 127 127 cyan_127);
  my @cyan_096 = qw(0 96 96 cyan_096);
  my @blue = qw(0 0 255 blue);
  my @blue_191 = qw(0 0 191 blue_191);
  my @blue_127 = qw(0 0 127 blue_127);
  my @blue_096 = qw(0 0 96 blue_096);
  my @blue_064 = qw(0 0 64 blue_064);
  my @magenta = qw(255 0 255 magenta);
  my @magenta_191 = qw(191 0 191 magenta_191);
  my @red = qw(255 0 0 red);
  my @red_191 = qw(191 0 0 red_191);
  my @red_127 = qw(127 0 0 red_127);
  my @red_096 = qw(96 0 0 red_096);
  my @red_064 = qw(64 0 0 red_064);
  my @pink = qw(255 80 180 pink);
  my @orange = qw(255 100 0 orange);
  my @violet = qw(120 0 200 violet);
  my @brown = qw(100 31 31 brown);
  my @pistache = qw(100 225 150 pistache);
  my @violet_pale = qw(230  215  255 violet_pale);
  my @pink_pale = qw(255 230 210 pink_pale);
  my @champagne = qw(255 240 200 champagne);
  my @pistache_pale = qw(200 255 200 pistache_pale);

  my @palette = ();
  push @palette, \@blue;
  push @palette, \@red;
  push @palette, \@green_175;
  push @palette, \@pink;
  push @palette, \@cyan_200;
  push @palette, \@orange;
  push @palette, \@violet;
  push @palette, \@gray_125;
  push @palette, \@brown;
  push @palette, \@black;
  push @palette, \@yellow_225;
  push @palette, \@cyan;
  push @palette, \@pistache;
  push @palette, \@magenta;
  push @palette, \@green;
  push @palette, \@green_127;
  push @palette, \@gray_150;
  push @palette, \@red_191;
  push @palette, \@blue_127;
  push @palette, \@cyan_096;
  push @palette, \@blue_096;
  push @palette, \@red_096;
  push @palette, \@blue_064;
  push @palette, \@gray_064;

  &RSAT::message::Info("Defined RGB palette with ".scalar(@palette)." colors") if ($main::verbose >= 4);
  if ($main::verbose >= 4) {
    foreach my $c (0..$#palette) {
      my $color_ref = $palette[$c];
      my ($R,$G,$B, $name) = @{$color_ref};
      &RSAT::message::Debug("Color palette", $c, $R,$G,$B,$name) if ($main::verbose >= 0);
    }
  }
  return @palette;
}


return 1;


__END__


