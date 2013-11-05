###############################################################
#
# Wrapper to use PostScript::Simple with the same functions as GD
#
package RSAT::PostscriptWrapper;

use PostScript::Simple;

=pod

=head1 NAME

    RSAT::PostscriptWrapper

=head1 DESCRIPTION

A wrapper to export images as postscripts, whilst using the GD
interface. GD functions are substituted for calls to
PostScript::Simple functions. 

Warning: this is not a full wrapper: I only defined some basic
fonctionalities of GD, which were needed by feature-map.

Since the authors of PostScript::Simple envisage to implement their
own compatibility between GD and PostScript, this is only a temporary
fix.

=over

=cut


################################################################
=pod

=item new()

Create a new image. This object is a wrapper around a
PostScript::Simple object. All arguments are passed to
PostScript::Simple creator.

=cut

sub new {
    my ($class, $graph_x_size, $graph_y_size, %args) = @_;
    my $self = bless {
    }, $class;
    my $landscape = 0;
#    if ($graph_x_size > $graph_y_size) { $landscape=1 ; }

    $self::image = new PostScript::Simple(xsize=>$graph_x_size, 
					  ysize=>$graph_y_size,
					  colour=>1,
					  eps=>1,
#					  landscape=>$landscape,
#					  clip=>0,
					  direction=>"RightUp",
					  coordorigin=>'LeftTop', ## Warning: This option seems not to work. This is problematic since the default is not the same as for GD.
					  units=>"bp",
					  %args);
    $self::image->setlinewidth(1);
    $self::ysize = $graph_y_size;

    return $self;
}


################################################################
=pod

=item rectangle()

Draw a rectangle

=cut

sub rectangle {
    my ($self, $left, $top, $right, $bottom, $color) = @_;
    $top = $self->convert_y($top);
    $bottom = $self->convert_y($bottom);
    $self->setcolour($color);
    $self::image->box({filled=>0}, $left, $top, $right, $bottom)
}



################################################################
=pod

=item filledRectangle()

Draw a filled rectangle

=cut

sub filledRectangle {
    my ($self, $left, $top, $right, $bottom, $color) = @_;
    $top = $self->convert_y($top);
    $bottom = $self->convert_y($bottom);
    $self->setcolour($color);
    $self::image->box({filled=>1}, $left, $top, $right, $bottom)
}



################################################################
=pod

=item line()

Draw a line.

=cut

sub line {
    my ($self, $left, $top, $right, $bottom, $color) = @_;
    $top = $self->convert_y($top);
    $bottom = $self->convert_y($bottom);
    $self->setcolour($color);
    $self::image->line($left, $top, $right, $bottom)
}


################################################################
=pod

=item arc()

Draw an arc or a circle.

=cut

sub arc {
    my ($self, $x,$y,$width, $height,$start,$end,$color, $filled) = @_;
    $y = $self->convert_y($y);
    my $radius = ($width+$height)/4;
    $self->setcolour($color);
    $self::image->arc({filled=>$filled}, $x,$y, $radius, $start, $end);
}


################################################################
=pod 

=item filledArc()

Draw a filled arc or a circle.

=cut

sub filledArc {
    my ($self, $x,$y,$width, $height,$start,$end,$color) = @_;
    $self->arc($x,$y,$width, $height,$start,$end,$color, 1);
}


################################################################
=pod 

=item polygon()

Draw a polygon. 

=cut

sub polygon {
    my ($self, $filled, $color, @points) = @_;
    warn join ("\t", ";", "drawing polygon", $filled, $color), "\n" if ($main::verbose >= 10);
    for my $i (1..$#points) {
	next unless ($i%2);
	$points[$i] = $self->convert_y($points[$i]);
    }
    $self->setcolour($color);
    $self::image->polygon({filled=>$filled}, @points);
}


################################################################
=pod 

=item filledPolygon()

Draw a filled polygon. 

=cut

sub filledPolygon {
    my ($self, $color, @points) = @_;
    $self->polygon(1, $color, @points);
}


################################################################
=pod

=item fillToBorder()

NOT YET IMPLEMENTED

=cut

sub fillToBorder {
    my ($x,$y,$border_color,$color);
    warn "fillToBorder() Not yet implemented for postscript format\n" if ($main::verbose >= 1);
}



################################################################
=pod

=item string()

Draw a string text

=cut

sub string {
    my ($self, $font, $x,$y, $string, $color) = @_;
    $y = $self->convert_y($y);
    warn join( "\t", ";", "text", $font, $x, $y, $string, $color), "\n" 
	if ($main::verbose >= 10);

    my ($font_name, $font_size) = split /\s+/, $font;
    $font_name = "Arial" unless ($font_name);
    $font_size = 10 unless ($font_size);

    $y -= $font_size; ## Correct differences between GD and ps font position

    $self::image->setfont($font_name, $font_size);
    $self::image->setcolour(0,0,0);
    $self::image->text($x,$y, $string);
}

################################################################
=pod

=item colorAllocate()

Define a color and associate it to a name

=cut

sub colorAllocate {
    my ($self, $red, $green, $blue) = @_;
    my $name = join " ", $red, $green, $blue;
    warn join( "\t", ";", "allocating color", $red, $green, $blue), "\n" 
	if ($main::verbose >= 10);
    return $name;
}


################################################################
#  =pod

=item convert_y()

Convert Y positions (in GD, 0 is the top, in ps, 0 is the bottom). 

=cut

sub convert_y {
    my ($self, $y) = @_;
    
    my $height = $self::ysize;
    if ($height > 0) {
	$y =  $height - $y;
    } else {
	warn "ysize is not defined. The image will be mirrored vertically\n";
    }
    return $y;
}

################################################################
## Calls for methods from PostScript::Simple

################################################################
=pod

=item setcolour()

Set the current color.

=cut

sub setcolour {
    my ($self, $colour) = @_;
    my ($red,$green,$blue) = split /\s+/, $colour;

    warn join("\t", ";", "setting colour", $red, $green, $blue), "\n" if ($main::verbose >= 10);
    $self::image->setcolour($red,$green,$blue);
}

################################################################
=pod

=item output()

Print the output in a ps file

=cut

sub output {
    my ($self, $outputfile) = @_;
    if ($outputfile) {
      $self::image->output($outputfile);
    } else {
      $page = $self::image->_builddocument($self, "title");
      foreach $i (@$page) {
	if (ref($i) eq "SCALAR") {
	  print STDOUT $$i;
	} else {
	  print STDOUT $i;
	}
      }
    }
}

return 1;

__END__

=back
