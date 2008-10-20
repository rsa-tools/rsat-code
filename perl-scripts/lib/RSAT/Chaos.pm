###############################################################
#
# Class Chaos.pm
#
# developed by Morgane Thomas-Chollier
#
package RSAT::Chaos;

use RSAT::GenericObject;
use RSAT::error;
use Data::Dumper;
@ISA = qw( RSAT::GenericObject );

### class attributes

##coords
local $coord_x;
local $coord_y;

##declare default coord adress
local $xh ;
local $xb;
local $yl;
local $yr ;

=pod

=head1 NAME

    RSAT::Chaos

=head1 DESCRIPTION

Class for handling Chaos Game Representation

=cut


################################################################
=pod

=item B<new()>

Create a new Chaos object.

=cut
sub new {
    my ($class, %args) = @_;
    my $object = bless {
	}, $class;
    return $object;
}

################################################################
=pod

=item B<init_mapping_table()>

Create a new Chaos object.

=cut
sub init_mapping_table {
    my ($self) = @_;
   	my %letter_coord = ();
   	
   	$self -> set_attribute("start_flag",0);

	##initialise default coord values
	$$xh= 0 ;
	$$xb= 1 ;
	$$yl= 0 ;
	$$yr= 1 ;
	
	## universal coord for each letter (mapping)
	$letter_coord{"a"}->{"X"} = $xb;
	$letter_coord{"a"}->{"Y"} = $yl;
	$letter_coord{"c"}->{"X"} = $xh;
	$letter_coord{"c"}->{"Y"} = $yl;
	$letter_coord{"g"}->{"X"} = $xh;
	$letter_coord{"g"}->{"Y"} = $yr;
	$letter_coord{"t"}->{"X"} = $xb;
	$letter_coord{"t"}->{"Y"} = $yr;
   
   $self->set_hash_attribute("letter_coord", %letter_coord);
}

################################################################
=pod

=item B<init_mapping_table()>

Initialise chaos data matrix

=cut
sub init_chaos_table {
    my ($self,$w_length) = @_;

   	## initialise chaos data matrix
	my @chaos =();
	for (my $x = 0; $x < (2**$w_length); $x++) {
		for (my $y = 0; $y < (2**$w_length); $y++) {
			$chaos[$x][$y] = 0;
		}
	}	
	$self->set_array_attribute("chaos_table", @chaos);
}

################################################################
## Calculate the coordinate (x,y) of a given word
sub getCoordinate {
	my ($self,$word) = @_;
	my @chars = split(//, $word);
	my $letter = pop(@chars);
	
	my $new_word =  join('',@chars);

	## handle first the last letter
	if ($self -> {start_flag}) {
		$coord_x = $self->{letter_coord}->{$letter}->{"X"};	
		$coord_y = $self->{letter_coord}->{$letter}->{"Y"};
		$self -> {start_flag} = 0;
	
	## process remaining letters
	} else {	
		$$xh = $$coord_x*2;
		$$xb = $$xh +1;
		$$yl = $$coord_y *2;
		$$yr = $$yl +1;
				
		$coord_x = $self->{letter_coord}->{$letter}->{"X"};
		$coord_y = $self->{letter_coord}->{$letter}->{"Y"};	
		}
		
		## stop condition
		unless ($new_word eq "") {
			($coord_x,$coord_y)  = $self-> getCoordinate($new_word);		
		} else {
			return ($self->{letter_coord}->{$letter}->{"X"},$self->{letter_coord}->{$letter}->{"Y"});
		}
}

################################################################
=pod

=item B<get_chaos_table>

Return the filled chaos table. If necessary, run
$self->fill_chaos()

=cut
sub get_chaos_table {
  my ($self,$return_value,$oligo_freq_ref) = @_;
  unless ($self->get_attribute("chaos_calculated")) {
    $self->fill_chaos($return_value,$oligo_freq_ref);
  }
  return $self->get_attribute("chaos_table");
}

################################################################
=pod

=item B<fill_chaos>

Fill the chaos table.

=cut
sub fill_chaos {
  my ($self,$return_value,$oligo_freq_ref) = @_;
  
  my %oligo_freq = %$oligo_freq_ref;
  foreach my $word (keys(%oligo_freq)){
		$self -> {start_flag} = 1;	
		$word = lc($word);

		##re-initialise default value
		$$xh= 0 ;
		$$xb= 1 ;
		$$yl= 0 ;
		$$yr= 1 ;
		
		##get the coord (x,y)
		($coord_x,$coord_y) = $self-> getCoordinate($word);	
		 &RSAT::message::Info("Coord of word ", $word, ": ($$coord_x,$$coord_y) ") if ($main::verbose >= 2);
		##store the frequency in the chaos table
		if ($return_value eq "word") {
			$self->{chaos_table} -> [$$coord_x][$$coord_y] = $word;
		}else{
			$self->{chaos_table} -> [$$coord_x][$$coord_y] = $oligo_freq{$word};
		}
	}
	$self->set_attribute("chaos_calculated",1);
}


return 1;

__END__

