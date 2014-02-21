###############################################################
#
# Class Classification
#
package RSAT::Classification;

use RSAT::GenericObject;
use RSAT::util;
use RSAT::Family;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

### class attributes

=pod

=head1 NAME

    RSAT::Classification

=head1 DESCRIPTION

This class stores a classification, i.e. a set of classes, each
containing a set of objects.

The classification can be a partition (i.e. each objects belogs to one
and only one class) or contain multiple assignations.

The class is equipped with methods for importing and exporting a
classification from various formats.

Hierarchical classifications can be loaded, but the inter-class
hierarchy is not treated in the current version.

=cut


################################################################

=pod

=item new()

Create a new Classification.

=cut
sub new {
    my ($class, %args) = @_;
    my $object = bless {
	}, $class;
    return $object;
}


################################################################

=pod

=item B<read_from_file>

Read the classification from a text file. 


 Title    : read_from_file
 Usage    : $classificaion->->read_from_file($input_file, $input_format)
 Function : Read the Classification from a text file. 
 Supported formats: mcl, tab, profiles, ermg, rnsc

=cut

sub read_from_file {
    my ($self, $input_file, $input_format, $name_file, %args) = @_;
    if ($input_format eq "mcl") {
	$self->read_mcl($input_file, %args);
    } elsif ($input_format eq "tab") {
	$self->read_tab($input_file, %args);
    } elsif ($input_format eq "profiles") {
	$self->read_profiles($input_file, %args);
    } elsif ($input_format eq "ermg") {
	$self->read_ermg($input_file, %args);
    } elsif ($input_format eq "rnsc") {
	$self->read_rnsc($input_file, $name_file, %args);
    } elsif ($input_format eq "mcode") {
        $self->read_mcode($input_file, $name_file, %args);
    } else {
	&RSAT::error::FatalError(join ("\t", "Classification::read_from_file", $input_format, "is not a supported input format"));
    }
}


################################################################

=pod

=item B<read_mcl>

Read the classification from a mcl text file.


 Title    : read_mcl
 Usage    : $classificaion->->read_mcl($input_file)
 Function : Read the Classification from a mcl text file.

=cut

sub read_mcl {
    my ($self, $input_file) = @_;
    ($main::in) = &RSAT::util::OpenInputFile($input_file);

    ## Load the classification
    my $class_number = 0;
    while (<$main::in>) {
	next unless (/\S/);
	chomp;
	$class_number++;
	my $class_name = "cl_".$class_number;
	my @labels = split(/\s+/);
	my $class = new RSAT::Family(name=>$class_name);
	&RSAT::message::Info("Reading", $class, $class_number, $class_name) if ($main::verbose >= 3);
	foreach my $label (@labels) {
	    $class->new_member($label);
	}
	$self->push_attribute("classes", $class);
    }
}
################################################################

=pod

=item B<read_mcode>

Read the classification from a MCODE output text file (exported from Cytoscape plugin).


 Title    : read_rnsc
 Usage    : $classificaion->->read_rnsc($input_file, $input_file_names)
 Function : Read the Classification from a rnsc text file.
 Argument : $input_file -> the output of RNSC
          : $input_file_names -> names of the nodes
 

=cut

sub read_mcode {
    
    my ($self, $input_file, $input_file_names) = @_;
    ($main::in) = &RSAT::util::OpenInputFile($input_file);
    ## Load the classification
    while (my $ligne = <$main::in>) {
      next if ($ligne !~ /^[0-9]+\t/);
      chomp $ligne;
      my @lignecp = split /\t/, $ligne;
      my $class_number = $lignecp[0];
      my $nodes_id_list = $lignecp[4];
      my @labels = split /, /, $nodes_id_list;
      
      if ((scalar @labels) > 1) { 
        my $class_name = "cl_".$class_number;
        my $class = new RSAT::Family(name=>$class_name);
        &RSAT::message::Info("Reading", $class, $class_number, $class_name) if ($main::verbose >= 3);
        foreach my $label (@labels) {
          if ($label ne "-1") {
            my $member_name = $label;
            $member_name = $id_names{$label} if defined($id_names{$label});
            $class->new_member($member_name);
          }
        }
        $self->push_attribute("classes", $class);
      }
  }
}
################################################################

=pod

=item B<read_rnsc>

Read the classification from a RNSC output text file.


 Title    : read_rnsc
 Usage    : $classificaion->->read_rnsc($input_file, $input_file_names)
 Function : Read the Classification from a rnsc text file.
 Argument : $input_file -> the output of RNSC
          : $input_file_names -> names of the nodes
 

=cut

sub read_rnsc {
    my ($self, $input_file, $input_file_names) = @_;
    ($main::in) = &RSAT::util::OpenInputFile($input_file);
    ($main::names) = &RSAT::util::OpenInputFile($input_file_names) if defined($input_file_names);
    # Read the nodes names
    my %id_names = ();
    if (defined $input_file_names) {
      while (my $ligne = <$main::names>) {
        next if ($ligne =~ /^;/);
        next if ($ligne =~ /^\#/);
        next unless ($ligne =~ /\S/);
        chomp $ligne;
        my @lignecp = split /\t/, $ligne;
        $id_names{$lignecp[0]} = $lignecp[1]; 
      }
    }
    ## Load the classification
    my $class_number = 0;
    while (my $ligne = <$main::in>) {
      chomp $ligne;
      my @labels = split / /, $ligne;
      if ((scalar @labels) > 1) { 
        $class_number++;
        my $class_name = "cl_".$class_number;
        my $class = new RSAT::Family(name=>$class_name);
        &RSAT::message::Info("Reading", $class, $class_number, $class_name) if ($main::verbose >= 3);
        foreach my $label (@labels) {
          if ($label ne "-1") {
            my $member_name = $label;
            $member_name = $id_names{$label} if defined($id_names{$label});
            $class->new_member($member_name);
          }
        }
        $self->push_attribute("classes", $class);
      }
  }
}
################################################################

=pod

=item B<read_tab>

Read the classification from a tab-delimited text file.  Each row
describes one class membership.  By default, the first column
indicates the member, the second column the class. 

 Title    : read_tab
 Usage    : $classificaion->->read_tab($input_file)
 Function : Read the Classification from a tab-delimited text file.

The columns can be specified by specifying additional arguments. For
instance, if the input is a tab-delimited file with members in the
third and classes in the seventh column :

 Usage :  $classificaion->->read_tab($input_file,member_column=>3,class_column=>7)

A score column can optionally be specified with the argument score_column. 

 Usage :  $classification->->read_tab($input_file,score_column=3)

=cut
sub read_tab {
    my ($self, $input_file, %args) = @_;
    my ($in) = &RSAT::util::OpenInputFile($input_file);
    my %class = (); ## classes indexed by name
    $member_column = $args{member_column} || 1;
    $class_column = $args{class_column} || 2;
    $score_column = $args{score_column} || 0;
    
    ## Verbosity
    &RSAT::message::TimeWarn("Reading classification from tab file", $input_file) if ($main::verbose >= 2);
    &RSAT::message::Info(join("\n\t",
			      "Columns",
			      "member column: ".$member_column,
			      "class column: ".$class_column,
			      "score column: ".$score_column,
			     )) if ($main::verbose >= 3);

    ## Load the classification
    my $line = 0;
    while (<$in>) {
	$line++;
	next if (/^;/);
	next if (/^\#/);
	next unless (/\S/);
	chomp();
	s/\r$//;
	my @fields = split /\t/;
	
	### class member
	my $member_name = &RSAT::util::trim($fields[$member_column -1]);
	if ($args{toupper}) {
	  ## Case-insensitive comparison
	  $member_name = uc($member_name);
	}
	unless ($member_name =~ /\S/) {
	    &RSAT::message::Warning(join ("\t", "Error class file", 
					  $class_file,  "line", 
					  $line, "member not specified")) if ($main::verbose >= 1);
	    next;
	}
	
	### class name
	my $class_name = &RSAT::util::trim($fields[$class_column - 1]);
	if ($args{toupper}) {
	  ## Case-insensitive comparison
	  $class_name = uc($class_name);
	}
	unless ($class_name) {
	    &RSAT::message::Warning(join("\t", "Error class file", 
					 $class_file,  "line", 
					 $line, "class not specified")) if ($main::verbose >= 1);
	    next;
	}
	
	#### create a new class if required
	unless ($class{$class_name}) {
	    $class{$class_name} = new RSAT::Family(name=>$class_name);
	    $self->push_attribute("classes", $class{$class_name});
	}
	
	## Add the member to the class
	if ($score_column) {
	    ### read score
	    local $score = &RSAT::util::trim($fields[$score_column - 1]);
#	    unless (&RSAT::util::IsReal($score)) {
#		&RSAT::error::FatalError(join("\t", $score, "Invalid score (must be a Real number).", 
#					      "class file", $class_file,  
#					      "line", $line,
#					      "member", $member_name,
#					      "class", $class_name,
#					     ));
#	    }
	    
	    ## Check threshold on score
	    next if ((defined($args{min_score})) && ($score < $args{min_score}));
	    $class{$class_name}->new_member($member_name, 0, score=>$score);
#		&RSAT::message::Debug("member score", $class_name, $member_name, $score) if ($main::verbose >= 0);
	} else {
	    $class{$class_name}->new_member($member_name);
	}
	&RSAT::message::Debug( join ("\t",  $class_name,
				     $member_name,
				     $name{$id}) ) if ($main::verbose >= 5);
	
	if (($line % 10000 == 0) &&
	    ($main::verbose >= 2)) {
	  if (defined($input_file)) {
	    &RSAT::message::TimeWarn( "Reading classes from file", $input_file, "lines read",  $line);
	  } else {
	    &RSAT::message::TimeWarn( "Reading classes from", "STDIN", "lines read",  $line);
	  }
	}
#	&RSAT::message::Debug("Scores", $class_name, $class{$class_name}->get_attribute("scores")) if ($main::verbose >= 0);
	
    }
    close $in if ($input_file);
}

################################################################

=pod

=item B<read_profiles>

Read the classification from a profile file.  Each row
describes one member, each column one class. 

 Title    : read_profiles
 Usage    : $classificaion->->read_profiles($input_file)
 Function : Read the Classification from a profile file.

Rows starting with ";" are ignored.

The first row is suppose to contain the header (description of
classes), describing the member names. This row must start with a "#'.

The first word off each following row contains the class name. 

=cut

sub read_profiles {
    my ($self, $input_file, %args) = @_;
    my ($in) = &RSAT::util::OpenInputFile($input_file);

    my $null = 0;
    if (defined($args{null})) {
	$null = $args{null};
    }

    my %class = (); ## classes indexed by name
    my $got_header = 0;

    ## Verbosity
    if ($main::verbose >= 1) {
	&RSAT::message::TimeWarn("Reading classification from profile file", $input_file) if ($main::verbose >= 2);
    }


    ## Load the classification
    my $line = 0;
    my @class_names = ();
    while (<$in>) {
      $line++;
      next if (/^;/); ## Skip comment lines
      next unless (/\S/); ## Skip empty lines
      s/\r//; ## Replace DOS-specific carriage return characters by Unix carriage return
      chomp(); ## Suppres carriage return


      ## Read class names from the header line, and create all the classes
      unless ($got_header) {
	if (/^\#/) {
	  $got_header = 1;
	  @class_names = split /\t/;
	    shift @class_names;
	    &RSAT::message::Debug("line", $line, "Header", "class names",join( "; ", @class_names)) if ($main::verbose >= 5);
	    foreach my $c (0..$#class_names) {
	      my $class_name = &RSAT::util::trim($class_names[$c]);

	      ### Check class name
	      unless ($class_name) {
		&RSAT::error::FatalError("Profile file", $class_file,
					 "line", $line,
					 "class nb", $c+1,
					 "file column", $c+2,
					 "empty class name in the header");
	    }

	      ### Create a new class if required
	    if ($class{$class_name}) {
	      &RSAT::error::FatalError("The header contains several columns with the same name (".$classname.")");
	    } else {
	      $class{$class_name} = new RSAT::Family(name=>$class_name);
	      $self->push_attribute("classes", $class{$class_name});
	    }
	  }
	  next;
	} else {
	  &RSAT::error::FatalError("Reading classes from profile file. The first non-comment line should contain the header and start with #");
	}
      }

	my @fields = split /\t/;
	### class members
	$member_name = shift @fields;
	&RSAT::message::Info(join("\t", "Treating row", $line, $member_name)) if ($main::verbose >= 5);
	unless ($member_name =~ /\S/) {
	    &RSAT::message::Warning(join ("\t", "Error class file", 
					  $class_file,  "line", 
					  $line, "member not specified")) if ($main::verbose >= 1);
	    next;
	}
	foreach my $class_name (@class_names) {
	    my $score = shift (@fields);
	    if (($score) && ($score ne $null)) {
		unless (&RSAT::util::IsReal($score)) {
		    &RSAT::error::FatalError(join("\t", $score, "Invalid score (must be a Real number).", 
						  "class file", $class_file,  
						  "line", $line,
						  "member", $member_name,
						  "class", $class_name,
						 ));
		}
		## Check threshold on score
		next if ((defined($args{min_score})) && ($score < $args{min_score}));
		&RSAT::message::Info(join("\t", "New member", $line, $class_name, $member_name, $score)) if ($main::verbose >= 4);
		$class{$class_name}->new_member($member_name, 0, score=>$score);
	    }
	}
	
    }
    close $in if ($input_file);

}


################################################################

=pod

=item B<read_ermg>

Read the classification from the output file of ERMG algorithm
developed by Stephane Robin (submitted).

 Title    : read_ermg
 Usage    : $classificaion->->read_ermg($input_file)
 Function : Read the Classification from an output file of the ERMG program.

=cut

sub read_ermg {
    my ($self, $input_file, %args) = @_;
    my ($in) = &RSAT::util::OpenInputFile($input_file);
    my %class = (); ## classes
    my $all_scores = 0;
    if (defined($args{all_scores})) {
	$all_scores = 1;
    }
    
    ## Verbosity
    &RSAT::message::TimeWarn(join("\t", 
				  "Reading classification from ERMG file", 
				  $input_file) ) if ($main::verbose >= 2);

    ## Load the classification
    my $line = 0;
    my $started = 0;
    my $clusters = 0;
    while (<$in>) {
	$line++;
	next unless (/\S/);
	chomp();
	s/\r$//;
	
	if (/^MAP \+ Posterior probabilities/) {
	    $started = 1;
	    &RSAT::message::Debug("Starting to read class elements", "line", $line) if ($main::verbose >= 5);
	    next;
	} elsif (/Number of groups: Q = (\d+)/) {
	    ## Number of clusters
	    $clusters = $1;
	    unless ((&RSAT::util::IsNatural($clusters)) && ($clusters >= 1)) {
		&FatalError(join("\t", "read_ermg", $clusters, "Invalid number of clustrers, must be a strictly positive natural number."));
	    }
	    
	    #### create all the classes
	    for my $cl (1..$clusters) {
		&RSAT::message::Info(join("\t", "Creating class for cluster", $cl)) if ($main::verbose >= 4);
		my $class = new RSAT::Family(name=>$cl);
		$self->push_attribute("classes", $class);
		$classes{$cl} = $class;
	    }
	    next;
	}
	next unless ($started);

	my @fields = split /\s+/;
	
	### class member
	$member_name = shift(@fields);
	unless ($member_name =~ /\S/) {
	    &RSAT::message::Warning(join ("\t", "Error class file", $class_file,  
					  "line", $line, "element not specified")) if ($main::verbose >= 1);
	    next;
	}
	
	### class name
	$class_name = shift(@fields);
	unless ((&RSAT::util::IsNatural($class_name)) && ($class_name >= 0)) {
	    &RSAT::message::Warning(join("\t", "Error class file", $class_file,  
					 "line", $line, "class not specified")) if ($main::verbose >= 1);
	    next;
	}
	if ($class_name > $clusters) {
	    &RSAT::message::Warning(join("\t", "Error class file", $class_file,  
					 "line", $line, "class number is larger than the number of clusters")) if ($main::verbose >= 1);
	    next;	    
	}
	

	## Use the posterior probabiliy as score
	for my $cl (1..$clusters) {
	    my $score = shift @fields;

	    ## Check threshold on score
	    next if ((defined($args{min_score})) && ($score < $args{min_score}));

	    if (($cl == $class_name) || ($all_scores)) {
		## Add the member to the class
		$classes{$cl}->new_member($member_name, 0, score=>$score);
		&RSAT::message::Info(join("\t", "Element", $member_name, 
					  "Class", $class_name,
					  "Score", $score
					 )) if ($main::verbose >= 4);
	    }
	}
    }
    close $in if ($input_file);

}



################################################################

=pod

=item B<to_text>

Convert the classification into a string for export to a text file.


 Title    : to_text
 Usage    : $classificaion->->to_text($output_format)
 Supported formats: tab

=cut

sub to_text {
    my ($self, $output_format, %args) = @_;
    if ($output_format eq "tab") {
	$self->to_tab(%args);
    } elsif ($output_format eq "mcl") {
	$self->to_mcl(%args);
    } elsif ($output_format eq "profiles") {
	$self->to_profiles(%args);
    } else {
	&RSAT::error::FatalError(join ("\t", "Classification::to_text", $output_format, "is not a supported output format"));
    }
}

################################################################

=pod

=item B<to_tab>

Convert the classification to a tab-delimited string. 
One row per element-class association. 
First column contains the element.
Second column contains the class.

 Title    : to_tab
 Usage    : $classificaion->->to_tab()

=cut

sub to_tab {
    my ($self) = @_;
    my @classes = $self->get_attribute("classes");
    my $string = "";
    my $some_scores = 0;

    foreach my $class (@classes) {
	my %scores = $class->get_attribute("scores");
	my $class_name = $class->get_attribute("name");
	foreach my $member ($class->get_members()) {
	    $string .= join ("\t", $member, $class_name);
	    if ($scores{$member}) {
		$some_scores = 1;
		$string .= "\t".$scores{$member};
	    }
	    $string .= "\n";
	}
    }

    my $header = join ("\t", "#element", "class");
    if ($some_scores) {
	$header .= "\tscore";
    };
    $header .= "\n";
    $string = $header.$string;
    return($string);
}

################################################################

=pod

=item B<to_profiles>

Convert the classification to a profile file. 
One row per element. 
One column per class.

 Title    : to_profiles
 Usage    : $classificaion->->to_tab()

=cut

sub to_profiles {
    my ($self, %args) = @_;
    my @classes = $self->get_attribute("classes");
    my $null = 0;
    if (defined($args{null})) {
	$null = $args{null};
    }

    ## Index class-member associations
    my %cross_table = ();
    my @class_names = ();
    my $class_name;
    my $member;
    foreach my $class (@classes) {
	my %scores = $class->get_attribute("scores");
#	&RSAT::message::Debug( "SCORES", "class", $class_name, "members", %scores);
	$class_name = $class->get_attribute("name");
	push @class_names, $class_name;
	foreach $member ($class->get_members()) {
	    if (defined($scores{$member})) {
		$cross_tab{$member}{$class_name} = $scores{$member};		
	    } else {
		$cross_tab{$member}{$class_name}++;
	    }
	}
    }
    my @members = sort(keys %cross_tab);

    ## Header
    my $string = join ("\t", "#", @class_names);
    $string .= "\n";


    ## Print profiles
    foreach $member (@members) {
	$string .= $member;
	foreach $class_name (@class_names) {
	    $string .= "\t";
	    if ($cross_tab{$member}{$class_name}) {
		$string .= $cross_tab{$member}{$class_name};
	    } else {
		$string .= $null;
	    }
	}
	$string .= "\n";
    }
    return $string;
}

################################################################

=pod

=item B<to_mcl>

Convert the classification to the mcl cluster format. MCL is a
graph-based clustering algorithm developed by Stijn Van Dongen
(http://micans.org/mcl/).


One row per cluster. Each row provides a list of elements, separated
by tab characters. 

 Title    : to_mcl
 Usage    : $classificaion->->to_mcl()

=cut

sub to_mcl {
    my ($self) = @_;
    my @classes = $self->get_attribute("classes");
    my $string = "";
    foreach my $class (@classes) {
	my $class_name = $class->get_attribute("name");
#	&RSAT::message::Debug("Printing", $class, $class_name) if ($main::verbose >= 5);
	my @members = $class->get_members();
	$string .= join ("\t", @members);
	$string .= "\n";
    }
    return $string;
}

return 1;

__END__

