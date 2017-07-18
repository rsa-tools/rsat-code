##############################################################
#
# Class MatrixReader
#
package RSAT::MatrixReader;

use RSAT::GenericObject;
use RSAT::matrix;
use RSAT::error;
use RSAT::feature;
@ISA = qw( RSAT::GenericObject );

### class attributes


=pod

=head1 NAME

RSAT::MatrixReader

=head1 DESCRIPTION

A class for reading position-specific scoring matrices from different
formats.

=cut


################################################################
## Class variables
%supported_input_format = (
			   'alignace'=> 1,
			   'assembly'=>1,
			   'cb'=>1,
			   'clustal'=>1,
			   'cluster-buster'=>1,
			   'consensus'=>1,
			   'sequences'=>1,
			   'feature'=>1,
			   'footprintdb'=>1,
			   'gibbs'=> 1,
			   'infogibbs'=>1,
			   'info-gibbs'=>1,
			   'jaspar'=>1,
			   'homer'=>1,
			   'mscan'=>1,
			   'meme'=>1,
			   'meme_block'=>1,
			   'motifsampler'=>1,
			   'stamp'=>1,
			   'stamp-transfac'=>1,
			   'tab'=>1,
			   'tf'=>1,
			   'transfac'=>1,
			   'cis-bp'=>1,
			   'uniprobe'=>1,
                           'yeastract'=>1,
			   'encode'=>1,
			  );
$supported_input_formats = join ",", sort(keys %supported_input_format);

=pod

=item List supported input formats

List available input formats.

=cut
sub ListInputMatrixFormats {
  return(%supported_input_format);
}
sub ListMatrixFormats { ## FOR BACKWARD COMPATIBILITY
  return(%supported_input_format);
}

=pod

=item readFromFile($file, $format)

Read a matrix from a file.

Supported arguments.

=over

=item  top=>X

Only return the X top matrices of the file.

Example: 
 my @matrices =
    &RSAT::MatrixReader::readFromFile($matrix_file, $input_format, top=>3);

=item Other arguments

Other arguments are passed to some format-specific readers (tab,
cluster-buster).

=back

=cut
sub readFromFile {
    my ($file, $format, %args) = @_;
    $format = lc($format);
    $format =~ s/^cb$/cluster-buster/;
    $format =~ s/^tf$/transfac/;
    $format =~ s/^info-gibbs$/infogibbs/;

    my @matrices = ();

    if (($format eq "consensus") || ($format =~ /^wc/i)) {
	@matrices = _readFromConsensusFile($file);
    } elsif ($format eq "transfac") {
	@matrices = _readFromTRANSFACFile($file, "transfac");
    } elsif ($format eq "cis-bp") {
	@matrices = _readFromTRANSFACFile($file, "cis-bp");
    } elsif ($format eq "yeastract"){
	@matrices = _readFromTRANSFACFile($file, "yeastract");
    
#    } elsif ($format eq "stamp-previous") {
    } elsif ($format eq "stamp") { ## This is the format of STAMP demo matrices
	@matrices = _readFromSTAMPFile($file);
    } elsif ($format eq "stamp-transfac") { ## STAMP exports files in a transfac-like format
	@matrices = _readFromTRANSFACFile($file, "stamp"); 
    } elsif ($format eq "footprintdb") { ## footprintDB exports files in a transfac-like format
	@matrices = _readFromTRANSFACFile($file, "stamp"); 
    } elsif ($format eq "infogibbs") {
	@matrices = _readFromInfoGibbsFile($file);
    } elsif ($format eq "assembly") {
	@matrices = _readFromAssemblyFile($file);
    } elsif ($format eq "gibbs") {
	@matrices = _readFromGibbsFile($file);
    } elsif ($format eq "alignace") {
	@matrices = _readFromAlignACEFile($file);
    } elsif ($format eq "tab") {
	@matrices = _readFromTabFile($file, %args);
    } elsif ($format eq "cluster-buster") {
	@matrices = _readFromClusterBusterFile($file, %args);
    }elsif ($format eq "encode") {
	@matrices = _readFromEncodeFile($file, %args);
    } elsif (($format eq "jaspar") || ($format eq "mscan")) {
	@matrices = _readFromJasparFile($file, %args);
    } elsif ($format eq "homer") {
	@matrices = _readFromHomerFile($file, %args);
    } elsif ($format eq "uniprobe") {
	@matrices = _readFromUniprobeFile($file, %args);
    } elsif ($format eq "motifsampler") {
	@matrices = _readFromMotifSamplerFile($file);
    } elsif ($format eq "meme_block") {
	@matrices = _readFromMEMEFile($file);
	} elsif ($format eq "meme") {
	@matrices = _readFromMEMEFile_2015($file, %args);
    } elsif ($format eq "feature") {
	@matrices = _readFromFeatureFile($file);
      } elsif ($format eq "sequences") {
	@matrices = _readFromSeq($file, %args);
    } else {
	&main::FatalError("&RSAT::matrix::readFromFile", "Invalid format for reading matrix\t$format");
    }

    ################################################################
    ## Check that there was at least one matrix in the file
    if (scalar(@matrices) == 0) {
      ## In case not a single matrix has been read, check if file is
      ## empty.
      $msg_file = &RSAT::util::hide_RSAT_path($file);
      if (-z $file) {
	&RSAT::message::Warning("Matrix file is empty (file size is zero)", $msg_file);
      } else {
	&RSAT::message::Warning("Matrix file",  $msg_file, "does not contain any matrix in", $format , "format. Please check format. ");
      }
    }  else {
      &RSAT::message::Info("Read ".scalar(@matrices)." matrices from file ", $file) if ($main::verbose >= 3);
    }


    ################################################################
    ## Assign or re-assign general parameters.
    ##
    ## Depending on the input format, some of these have already been
    ## assigned from the input file but we prefer to reassign them for
    ## the sake of consistency.

    &RSAT::message::Info("Assigning matrix parameters") if ($main::verbose >= 3);
    my $matrix_nb = 0;
    foreach my $matrix (@matrices) {

      ## Some programs require sorted alphabet (A,C,G,T). We re-order
      ## the rows alphabetically.
      $matrix->sort_rows();

      ## Reassign matrix numbers
      $matrix_nb++;
      $matrix->set_parameter("matrix.nb", $matrix_nb);

      ## If a prefix is specified, use it in the name, ID and accession number
      if ($args{prefix}) {
	my $prefix = $args{prefix};

	if (scalar(@matrices) > 1) {
	  $m++;
	  $prefix .= "_m".$matrix_nb;
	}

	$matrix->force_attribute("accession", $prefix);
	$matrix->force_attribute("AC", $prefix);
	$matrix->force_attribute("name", $prefix);
	$matrix->force_attribute("id", $prefix);

      }

      ## The -prefix option deletes the AC of the matrix.  I (Jaime)
      ## think is more informative (and required for matrix-clustering)
      ## to save the old AC. So now the -prefix output is like this:
      ## ${prefix}_${old_AC} If a prefix is specified, use it in the
      ## name, ID and accession number The new ID is the prefix
      ## indictated + the current motif id.
      if ($args{prefix_id}) {

	my $prefix_id = $args{prefix_id};

	# my $old_AC = $matrix->get_attribute("accession");
	$old_AC = $matrix->get_attribute("accession");

	$prefix_id = $prefix_id."_m".$matrix_nb."_".$old_AC;

	$matrix->force_attribute("accession", $prefix_id);
	$matrix->force_attribute("AC", $prefix_id);
#	$matrix->force_attribute("ID", $prefix_id);
	# $matrix->force_attribute("id", $prefix_id);
	# $matrix->force_attribute("name", $prefix_id);
      }

      ## Check that each matrix contains at least one row and one col
      if ($matrix->nrow() > 0) {
	if ($matrix->ncol() > 0) {
	  &RSAT::message::Info("Matrix read",
			       "nrow = ".$matrix->nrow(),
			       "ncol = ".$matrix->ncol(),
			       "prior : ".join (" ", $matrix->getPrior()),
			      ) if ($main::verbose >= 5);

	  ## Count number of sites per matrix
	  my $site_nb = scalar($matrix->get_attribute("sequences"));
	  if ($site_nb) {
	    $matrix->set_parameter("sites", $site_nb);
	  }

	} else {
	  &RSAT::message::Warning("The matrix", $matrix_nb, $matrix->get_attribute("id"), "contains 0 columns.");
	}
      } else {
	&RSAT::message::Warning("The matrix", $matrix_nb, $matrix->get_attribute("id"), "contains 0 rows.");
      }

      ## Replace undefined values by 0
      $matrix->treat_null_values();

      ## Compute MAP per sites
      if ((my $map = $matrix->get_attribute("MAP")) && (my $sites = $matrix->get_attribute("sites"))) {
	$matrix->set_parameter("MAP.per.site", $map/$sites);
      }

      ## Suppress invalid characters from the matrix ID
      my $id = $matrix->get_attribute("id");
      if ($id) {
#	$id =~ s/[\(\)\/]/_/g;
	$id = &clean_id($id);
	$matrix->force_attribute("id", $id);
      }

      ## Suppress invalid characters from the matrix accession, and ensure that matrices always have an accession
      my $ac = $matrix->get_attribute("accession");
      if ($ac) {
#	  $ac =~ s/[\(\)\/]/_/g;
	$ac = &clean_id($ac);
      } else {
	  $ac = &clean_id($id);
      }
      $matrix->force_attribute("accession", $ac);
      $matrix->force_attribute("AC", $ac);
    }

    ## Only retain the N top matrices if requested
    my $top = 0;
    if (defined($args{top})) {
      $top = $args{top};
      if ((&RSAT::util::IsNatural($top)) && ($top > 0)) {
	my $matrix_nb = scalar(@matrices);
	if ($matrix_nb > $top) {
	  &RSAT::message::Info("Retaining", $top, "top matrices among", (scalar(@matrices))) if ($main::verbose >= 1);
	  foreach my $m (($top+1)..$matrix_nb) {
	    pop @matrices;
	  }
	}
      }
    }


    ## Skip the N first matrices if requested
    my $skip = 0;
    if (defined($args{skip})) {
      $skip = $args{skip};
      &RSAT::message::Info("Skipping", $skip, "top matrices among", (scalar(@matrices))) if ($main::verbose >= 1);
      if ((&RSAT::util::IsNatural($skip)) && ($skip > 0)) {
	for my $m (1..$skip) {
	  shift(@matrices);
	}
      }
    }



    ## Select matrices specified as list of IDs
    my @selected_ids = ();
    if (defined($args{selected_ids})) {
      @selected_ids = @{$args{selected_ids}};
      &RSAT::message::Info("Selected IDs", join (",", @selected_ids)) if ($main::verbose >= 2);
      if (scalar(@selected_ids) > 0) {
	## Index selected IDs in a hash table
	my %selected_id = ();
	foreach my $id (@selected_ids) {
	  $selected_id{lc($id)} = 1;
	}

	## Select matrices
	my @retained_matrices = ();
	foreach my $matrix (@matrices) {
	  my $matrix_id = $matrix->get_attribute("id");
	  push @retained_matrices, $matrix if ($selected_id{lc($matrix_id)});
	}
	@matrices = @retained_matrices;
      }
    }

    ## Select matrices specified as list of ACs
    my @selected_acs = ();
    if (defined($args{selected_acs})) {
      @selected_acs = @{$args{selected_acs}};
      &RSAT::message::Info("Selected ACs", join (",", @selected_acs)) if ($main::verbose >= 2);
      if (scalar(@selected_acs) > 0) {
	## Index selected ACs in a hash table
	my %selected_ac = ();
	foreach my $ac (@selected_acs) {
	  $selected_ac{lc($ac)} = 1;
	}

	## Select matrices
	my @retained_matrices = ();
	foreach my $matrix (@matrices) {
	  my $matrix_ac = $matrix->get_attribute("accession");
	  push @retained_matrices, $matrix if ($selected_ac{lc($matrix_ac)});
	}
	@matrices = @retained_matrices;
      }
    }


    ## Select matrices specified as list of names
    my @selected_names = ();
    if (defined($args{selected_names})) {
      @selected_names = @{$args{selected_names}};
      &RSAT::message::Info("Selected names", join (",", @selected_names)) if ($main::verbose >= 2);
      if (scalar(@selected_names) > 0) {
	## Index selected names in a hash table
	my %selected_name = ();
	foreach my $name (@selected_names) {
	  $selected_name{lc($name)} = 1;
	}

	## Select matrices
	my @retained_matrices = ();
	foreach my $matrix (@matrices) {
	  my $matrix_name = $matrix->get_attribute("name");
	  push @retained_matrices, $matrix if ($selected_name{lc($matrix_name)});
	}
	@matrices = @retained_matrices;
      }
    }


    ## Report remaining matrices after selection
    if ((($skip + $top) > 0)
	|| (scalar(@selected_acs) > 0)
	|| (scalar(@selected_ids) > 0)
	|| (scalar(@selected_names) > 0)
	) {
      &RSAT::message::Info("Number of matrices after selection",  scalar(@matrices)) if ($main::verbose >= 1);
      if ($main::verbose >= 2) {
	my @retained_acs = ();
	my @retained_names = ();
	foreach my $matrix (@matrices) {
#	  push @retained_acs, $matrix->get_attribute("accession");
#	  push @retained_names, $matrix->get_attribute("name");
	  push @retained_message, $matrix->get_attribute("accession")." (".$matrix->get_attribute("name").")";
	}
#	&RSAT::message::Info("Retained matrix ACs",  join ",", (@retained_acs));
#	&RSAT::message::Info("Retained matrix names",  join ",", (@retained_names));
	&RSAT::message::Info("Retained matrices",  join "; ", (@retained_message));
      }
    }

    ## Return the matrices
    return @matrices;
}


=pod

Suppress problematic characters from the ID (parentheses, spaces, pipe, ...)

=cut
sub clean_id {
  my ($id) = @_;
  $id =~ s/[\(\)\/\\\s\|]/_/g;
  return $id;
}

=pod

=item InitializeEquiPriors

Initialize prior residue frequencies as equiprobable alphabet.

=cut
sub InitializeEquiPriors {
  my @matrices = @_;
  &RSAT::message::Info("Initializing equiprobable priors for all matrices") if ($main::verbose >= 3);
  foreach my $matrix (@matrices) {
    my @alphabet = $matrix->getAlphabet();
    if (scalar(@alphabet) == 0) {
#      @alphabet = qw(a c g t);
      $matrix->set_alphabet_for_type();
    } else {
      $matrix->setAlphabet(@alphabet);
    }

#    $matrix->force_attribute("nrow", scalar(@alphabet));
    my %tmp_prior = ();
    my $prior = 1/scalar(@alphabet);
    foreach my $residue (@alphabet) {
      $tmp_prior{$residue} = $prior;
      #	&RSAT::message::Debug("initial prior", $residue, $prior) if ($main::verbose >= 10);
    }
    $matrix->setPrior(%tmp_prior);
    if ($main::verbose >= 5) {
      &RSAT::message::Debug("Read matrix with alphabet", join(":", $matrix->getAlphabet()));
      &RSAT::message::Debug("Initialized prior as equiprobable", join(":", $matrix->getPrior()));
      &RSAT::message::Debug("Matrix size", $matrix->nrow()." rows",  $matrix->ncol()." columns");
    }
  }
}


################################################################
## read matrices from matrix file list paths

sub readMatrixFileList {
    my ($mlist, $input_dir,$input_format) = @_;
    my @matrix_files = ();
    my @matrices = ();
    while (<$mlist>) {
	next if (/'^;'/);		# skip comment lines
	next if (/'^#'/);		# skip header lines
	next if (/'^--'/);	# skip mysql-type comment lines
	next unless (/\S/);	# skip empty lines
	my @fields = split /\s+/;
	my $matrix_file = $fields[0];
	push @matrix_files, $matrix_file;

    }
    close $mlist;
    &RSAT::message::Info("Read matrix list from file", $main::infile{mlist2}, scalar(@matrix_files), "matrices") 
      if ($main::verbose >= 2);

    if (scalar(@matrix_files >= 1)) {
      foreach my $matrix_file (@matrix_files) {
	my @matrices_from_file = &readFromFile($matrix_file, $input_format);
	foreach my $matrix (@matrices_from_file) {
	  my ($matrix_name) = &RSAT::util::ShortFileName($matrix_file);
	  $matrix_name =~ s/\.\w+$//; ## suppress the extension from the file name
	  unless (defined($matrix->get_attribute("name"))){
	    $matrix->set_attribute("name", $matrix_name);
	  }
	  push @matrices, $matrix;
	}
      }
    }else{
	&RSAT::error::FatalError("The matrix ist must contain at least one matrix file path."); 
    }
    return @matrices;
}


=pod

=item _readFromTRANSFACFile($file, $format)

Read a matrix from a TRANSFAC file. This method is called by the method 
C<readFromFile($file, "TRANSFAC")>.

=cut
sub _readFromTRANSFACFile {
  my ($file, $format) = @_;
  &RSAT::message::Info ("Reading matrix from TRANSFAC file", $file) if ($main::verbose >= 3);

  ## open input stream
  my ($in, $dir) = &main::OpenInputFile($file);
  my $current_matrix_nb = 0;
  my @matrices = ();
  my $matrix;
  my $command = "";
  my $ncol = 0;
  my $transfac_consensus = "";

  my $short_file_name;
  if ($file) {
    $short_file_name = &RSAT::util::ShortFileName($file);
  } else {
    $short_file_name = "input";
  }

  my %prior = ();
  my $l = 0;
  while (<$in>) {
    $l++;
    next unless (/\S/);
    next if (/^;/);
    s/\r//;
    chomp();
    my $version = "";

    ## Read the command line
    if (/^VV\s+/) {
      $version = $';		# '
      &RSAT::message::Warning("TRANSFAC file version", $version) if ($main::verbose >= 3);

      ## empty field separator
    } elsif (/^XX/) {

      ## Start a new matrix (one TRANSFAC file contains several matrices)
    } elsif (/^AC\s*(.*)/) {
      my $accession = &clean_id($1);
      $current_matrix_nb++;
      $transfac_consensus = "";
      $matrix = new RSAT::matrix();
      $matrix->set_parameter("program", "transfac");
      $matrix->set_parameter("matrix.nb", $current_matrix_nb);
      push @matrices, $matrix;
      unless ($accession) {
	&RSAT::message::Info("Empty accession, set to", $accession) if ($main::verbose >= 3);
	$accession = $short_file_name."_".$current_matrix_nb;
      }
      &RSAT::message::Info("Setting AC, accession, id and name", $accession) if ($main::verbose >= 5);
      $matrix->set_parameter("accession", $accession);
      $matrix->set_parameter("AC", $accession);

      ## TRANSFAC Accession number corresponds to our matrix ID
      $matrix->set_parameter("id", $accession);

      ## TRANSFAC ID corresponds to our matrix name. We temporarily
      ## set the name to AC, since some badly formated fomres
      ## (e.g. coming from STAMP) do not contain the ID field.
      $matrix->set_parameter("name", $accession);
      $matrix->set_parameter("version", $version);
      $ncol = 0;
#      next;

    } elsif ((/^NA\s*(.*)/) && ($format eq "yeastract") ) {
      my $accession = &clean_id($1);
      $current_matrix_nb++;
      $transfac_consensus = "";
      $matrix = new RSAT::matrix();
      $matrix->set_parameter("program", "yeastract");
      $matrix->set_parameter("matrix.nb", $current_matrix_nb);
      push @matrices, $matrix;
      unless ($accession) {
	&RSAT::message::Info("Empty accession, set to", $accession) if ($main::verbose >= 3);
	$accession = $short_file_name."_".$current_matrix_nb;
      }
      &RSAT::message::Info("Setting AC, accession, id and name", $accession) if ($main::verbose >= 5);
      $matrix->set_parameter("accession", $accession);
      $matrix->set_parameter("AC", $accession);

      ## TRANSFAC Accession number corresponds to our matrix ID
      $matrix->set_parameter("id", $accession);

      ## TRANSFAC ID corresponds to our matrix name. We temporarily
      ## set the name to AC, since some badly formated fomres
      ## (e.g. coming from STAMP) do not contain the ID field.
      $matrix->set_parameter("name", $accession);
      $matrix->set_parameter("version", $version);
      $ncol = 0;
#      next;
  
      ## Read prior alphabet from the matrix header (PO line)
      ## Equiprobable alphabet

    } elsif ((/^PO\s+/)  || (/^P0\s+/) || (/^Pos\s+/)) { ## 2009/11/03 JvH fixed a bug, in previous versions I used P0 (zero) instead of PO (big "o")
	## CIS-BP database matrices are similar to transfac but do not contain an AC or ID line.
	## Intialize CIS-BP matrix
	if  ($format eq "cis-bp" ){
	    $current_matrix_nb++;
	    $transfac_consensus = "";
	    $matrix = new RSAT::matrix();
	    $matrix->set_parameter("program", "cis-bp");
	    $matrix->set_parameter("matrix.nb", $current_matrix_nb);
	    push @matrices, $matrix;
	    &RSAT::message::Info("Empty accession, set to", $accession) if ($main::verbose >= 3) ;
	    $accession = $short_file_name."_".$current_matrix_nb;

	    ## TRANSFAC ID corresponds to our matrix name. 
	    ## CIS-BP does not contain AC, ID or name.
      
	    &RSAT::message::Info("Setting AC, accession, id and name", $accession) if ($main::verbose >= 5);
	    $matrix->set_parameter("accession", $accession);
	    $matrix->set_parameter("AC", $accession);
	    $matrix->set_parameter("id", $accession);
	    $matrix->set_parameter("name", $accession);
	    $matrix->set_parameter("version", $version);
	    $ncol = 0;

	}
      ## Check for a likely error: user has entered only the PO and
      ## frequencies, without AC and ID.
      if ($current_matrix_nb < 1) {
	my $msg_file = &RSAT::util::hide_RSAT_path($file);
	&RSAT::error::FatalError("Wrongly formatted matrix file ".$msg_file, "In TRANSFAC format, matrices should start with an accession line ('AC  ').");
      }

      my $header = $'; #'
      $header = RSAT::util::trim($header);

      ## Alphabet is parsed from the TRANSFAC matrix header (PO row)
      my @alphabet = split /\s+/, $header;
#      $matrix->setAlphabet_lc(@alphabet);
      $matrix->setAlphabet(@alphabet);
#      $matrix->set_alphabet_for_type(); ## TO TEST (JvH): do we prefer to set the alphabet based on TRANSFAC column header or according to the user-specified type ? 
      &RSAT::message::Debug("Alphabet", join(";",@alphabet)) if ($main::verbose >= 5);

      ## Check that prior has been specified
      unless ($matrix->get_attribute("prior_specified")) {
	foreach my $letter (@alphabet) {
	  $prior{lc($letter)} = 1/scalar(@alphabet) 
	    unless (defined($prior{$letter}));;
	}
	$matrix->setPrior(%prior);
      }

	
      ## Other matrix parameters

    } elsif ($matrix) {
      ## Sites used to build the matrix
      if (/^BS\s+/) {
	my $bs = $'; #'
	my ($site_sequence, $site_id) = split(/\s*;\s*/, $bs);
	#      my $site_sequence = $1;
	#      my $site_id = $2;
	if ($site_sequence) {
	  $matrix->push_attribute("sequences", $site_sequence);
	  if ($site_id) {
	    $matrix->push_attribute("site_ids", $site_id);
	  }
	}
	&RSAT::message::Info("TRANSFAC site", $site_id, $site_sequence) if ($main::verbose >= 5);
#      &RSAT::message::Debug("line", $l, "site", $site_sequence, $site_id, $bs) if ($main::verbose >= 10);

      ## Count column of the matrix file (row in TRANSFAC format)
      } elsif (/^(\d+)\s+/) {
	my $values = $'; #'
	$values = &RSAT::util::trim($values);
	my @fields = split /\s+/, $values;
	my $consensus_residue= "";
	if ($fields[$#fields] =~ /[A-Z]/i) {
	  $consensus_residue = pop @fields;
	  $transfac_consensus .= $consensus_residue;
	}
	&RSAT::message::Debug("line ".$l, "adding column", join (":", @fields)) if ($main::verbose >= 5);
	$matrix->addColumn(@fields);
	$ncol++;
	$matrix->force_attribute("ncol", $ncol);

	## TRANSFAC identifier corresponds to our name
      } elsif ((/^ID\s+(\S+)/)
	       || ((/^NA\s+(\S+)/) && ($format eq "stamp")) ## For TRANSFAC-like format exported by STAMP
	       || ((/^NA\s+(\S+)/) && ($format eq "footprintdb")) ## For TRANSFAC-like format exported by footprintDB
	      ) {

	## TRANSFAC identifier corresponds to our matrix name
	## Our ID is in the TRANSFAC field AC (accession)
	$matrix->set_parameter("name", $1);

	&RSAT::message::Info("TRANSFAC",
			     "AC -> id", $matrix->get_attribute("id"),
			     "ID -> name", $matrix->get_attribute("name")
			    ) if ($main::verbose >= 5);

	## Bound factor
      } elsif (/^BF\s+/) {
	my $factor_description = $';
	my $factor_id = "";
	$matrix->push_attribute("binding_factor_desc", $factor_description); #'
	if ($factor_description =~ /^(T\d+)/) {
	  $factor_id = $1;
	  $matrix->push_attribute("binding_factor", $factor_id); #'
	}
	&RSAT::message::Info("TRANSFAC binding factor", $factor_id, $factor_description) if ($main::verbose >= 5);
	$matrix->set_parameter("binding_factors", join(";", $matrix->get_attribute("binding_factor")));
	&RSAT::message::Info("TRANSFAC binding factors", $matrix->get_attribute("binding_factors")) if ($main::verbose >= 5);

	## Short factor description
      } elsif (/^SD\s+/) {
	$matrix->push_attribute("short_foactor_description", $'); #'
	&RSAT::message::Info("TRANSFAC short factor desc", $factor_id, $factor_description) if ($main::verbose >= 5);

	## Statistical basis
      } elsif (/^BA\s+/) {
	$matrix->set_parameter("statistical_basis", $'); #'

	## Matrix description
      } elsif (/^DE\s+/) {
	$matrix->set_parameter("description", $'); #'

	## Store the consensus at the end of the matrix
      } elsif (/^\/\//) {
	$matrix->set_parameter("transfac_consensus", $transfac_consensus);

	## Row containing other field
      } elsif (/^(\S\S)\s+(.*)/) {
	my $field = $1;
	my $value = $2;
	&RSAT::message::Warning("Not parsed", $field, $value) if ($main::verbose >= 5);

      } else {
	&RSAT::message::Warning("Skipped invalid row", $_);
      }
    }
  }
  close $in if ($file);

  &RSAT::message::Debug("&RSAT::MatrixReader::readFromTRANSFAC()", "Finished reading TRANSFAC matrices") 
      if ($main::verbose >= 5);

  return @matrices;

}


=pod

=item _readFromSTAMPFile($file)

Read a matrix from a STAMP file. This method is called by the method
C<readFromFile($file, "STAMP")>.

=cut
sub _readFromSTAMPFile {
  my ($file) = @_;
  &RSAT::message::Info ("Reading matrix from STAMP file", $file) if ($main::verbose >= 3);

  ## open input stream
#  my $in = STDIN;
  my ($in, $dir) = &main::OpenInputFile($file);
#  if ($file) {
#    open INPUT, $file;
#    $in = INPUT;
#  }
  my $current_matrix_nb = 0;
  my @matrices = ();
  my $matrix = "";
  my $command = "";
  my $ncol = 0;
  my $comment_nb = 0;
  my @pre_matrix_comments = ();

  my %prior = ();
  my $l = 0;
  my $last_field_name = 0; ## 1 if the last read field was "NA"

  while (<$in>) {
    $l++;
    next if (/^;/);
    s/\r//;
    chomp();
    my $version = "";

    ## Empty row is the matrix separator
    unless (/\S/) {
      @pre_matrix_comments = ();
      $matrix = "";
#      &RSAT::message::Debug("End of matrix", $current_matrix_nb, $matrix) if ($main::verbose >= 10);
      next;
    }

    ## in STAMP, the XX is used to put comments (contrary to TRANSFAC where it is a field separator)
    if (/^XX/) {
      if (/^XX\s+(\S+.+)/) {
#	&RSAT::message::Debug("STAMP comment", $current_matrix_nb, $matrix) if ($main::verbose >= 10);
	if ($matrix) { ## Ignore the XX preceding the DE field, since the matrix is not yet created
	  $comment_nb++;
	  $matrix->set_parameter("STAMP_comment_".$comment_nb, $1);
	} else {
	  push @pre_matrix_comments, $1;
	}
      }

      ## Start a new matrix (one STAMP file contains several matrices)
    } elsif ((/^(DE)\s+(.+)/) || (/^(NA)\s+(.+)/)) {
      ## Problem: STAMP has 2 formats (so called "TRANSFAC" and the
      ## "TRANSFAC-like", but actually the TRANSFAC is not properly speaking TRANSFAC.
      ##
      ## In early versions, there was no AC field, whereas this is an
      ## essential field for TRANSFAC).
      ##
      ## The beginning of a record is marked EITHER by "NA"
      ## ("TRANSFAC" format) OR by "DE" ("TANSFAC-like" format).
      ## We have to fiddle around to circumvent this ambiguity.
      ##
      ## Apparently in 2013 there is an AC field. I adapted
      ## _readFromTRANSFACFile() to support the new STAMP format. This
      ## one is maintained for backwards compatibility.

      my $field = $1;
      my $accession = &clean_id($2); ## Name and description are both acceptable as AC

      ## STAMP uses the description field as accession number
      $comment_nb = 0;
      &RSAT::message::Info("STAMP accession", $accession) if ($main::verbose >= 5);
      $current_matrix_nb++;
      $matrix = new RSAT::matrix();
      $matrix->set_parameter("matrix.nb", $current_matrix_nb);
      push @matrices, $matrix;


      my $name = "";
      my $description = "";
      if ($field eq "NA") {
	## If field "name" is defined, we are in the so-called
	## "TRANSFAC" format. The "name" is then used as accession +
	## as ID + as name
	$last_field_name = 1;
	$name = $2;
	$matrix->set_parameter("name", $name);
	$matrix->set_parameter("accession", $accession);
	$matrix->set_parameter("AC", $accession);
	$matrix->set_parameter("id", $accession);
      } else {
	## Field DE is always taken as description
	$description = $2;
	$matrix->set_parameter("description", $description);

	unless ($last_field_name) {
	  ## If no name has been specified earlier, we are in the
	  ## "TRANSFAC-like" format. In this case the description is
	  ## used as name, accession and ID.
	  $matrix->set_parameter("name", $accession);
	  $matrix->set_parameter("accession", $accession);
	  $matrix->set_parameter("AC", $accession);
	  $matrix->set_parameter("id", $accession);
	}
	$last_field_name = 0;
      }

      $matrix->set_parameter("program", "STAMP");
      $matrix->set_parameter("version", $version);
      $ncol = 0;
      while (my $comment = shift(@pre_matrix_comments)) {
	$comment_nb++;
	$matrix->set_parameter("STAMP_comment_".$comment_nb, $comment);
      }

      ## STAMP has no header to specify the alphabet. Columns are supposed to contain A,C,G,T respectively.
      # my @alphabet = qw(A C G T);
      # $matrix->setAlphabet_lc(@alphabet);
      $matrix->set_alphabet_for_type();
      &RSAT::message::Debug("Alphabet", join(";",@alphabet)) if ($main::verbose >= 5);

      ## Equiprobable alphabet
      ## Check that prior has been specified
      unless ($matrix->get_attribute("prior_specified")) {
	foreach my $letter (@alphabet) {
	  $prior{lc($letter)} = 1/scalar(@alphabet) 
	    unless (defined($prior{$letter}));;
	}
	$matrix->setPrior(%prior);
      }

      ## Other matrix parameters
    } elsif ($matrix) {
      ## Sites used to build the matrix
      if (/^BS\s+/) {
	my $bs = $'; #'
	my ($site_sequence, $site_id) = split(/\s*;\s*/, $bs);
	#      my $site_sequence = $1;
	#      my $site_id = $2;
	if ($site_sequence) {
	  $matrix->push_attribute("sequences", $site_sequence);
	  if ($site_id) {
	    $matrix->push_attribute("site_ids", $site_id);
	  }
	}
	&RSAT::message::Info("TRANSFAC site", $site_id, $site_sequence) if ($main::verbose >= 5);
#      &RSAT::message::Debug("line", $l, "site", $site_sequence, $site_id, $bs) if ($main::verbose >= 10);

      ## Count column of the matrix file (row in TRANSFAC format)
      } elsif (/^(\d+)\s+/) {
	my $values = $'; #'
	$values = &RSAT::util::trim($values);
	my @fields = split /\s+/, $values;
	my $consensus_residue= "";
	if ($fields[$#fields] =~ /[A-Z]/i) {
	  $consensus_residue = pop @fields;
	  $transfac_consensus .= $consensus_residue;
	}
	&RSAT::message::Debug("line ".$l, "adding column", join (":", @fields)) if ($main::verbose >= 5);
	$matrix->addColumn(@fields);
	$ncol++;
	$matrix->force_attribute("ncol", $ncol);

	## Rows containing other fields
      } elsif (/^(\S\S)\s+(.*)/) {
	my $field = $1;
	my $value = $2;
	&RSAT::message::Warning("Not parsed", $field, $value) if ($main::verbose >= 5);


      } else {
	&RSAT::message::Warning("Skipped invalid row", $_);
      }
    }
  }
  close $in if ($file);

  return @matrices;

}


=pod

=item _readFromInfoGibbsFile($file)

Read a matrix from a result file from I<info-gibbs> (implementation by
Matthieu Defrance, 2008). This method is called by the method
C<readFromFile($file, "InfoGibbs")>.

=cut
sub _readFromInfoGibbsFile {
    my ($file, %args) = @_;
    &RSAT::message::Info ("Reading matrices from info-gibbs file", $file) if ($main::verbose >= 3);

    ## open input stream
    my ($in, $dir) = &main::OpenInputFile($file);

    ## read header
    if ($args{header}) {
	$header = <$in>;
	$header =~ s/\r//;
	chomp ($header);
	$header =~ s/\s+/\t/g;
	@header = split "\t", $header;
	$matrix->push_attribute("header", @header);
    }

    ################################################################
    ## Initialize the matrix list
    my @matrices = ();
    my $matrix = new RSAT::matrix();
    $matrix->set_parameter("program", "infogibbs");
    $matrix->set_parameter("matrix.nb", $current_matrix_nb);
    push @matrices, $matrix;
    my $current_matrix_nb = 1;
    #    my $id = $file."_".$current_matrix_nb;
    my $id_prefix = $file || "matrix";
    my $id = &clean_id($id_prefix."_".$current_matrix_nb);
    $matrix->set_attribute("AC", $id);
    $matrix->set_attribute("id", $id);
    my $l = 0;
    my $matrix_found = 0;
    my $new_matrix = 0;
    my $no_motif = 1;
    my $read_seqs=0;
    while ($line = <$in>) {
      $l++;
      next unless ($line =~ /\S/); ## Skip empty lines
      chomp($line); ## Suppress newline
      $line =~ s/\r//; ## Suppress carriage return
      $line =~ s/(^.)\|/$1\t\|/; ## Add missing tab after residue
      $line_aux = $line;
      @aux_line=split (/ +/,$line_aux);

      if ($line =~ /^; random see/ ) {
	  $matrix->force_attribute("random_seed", pop (@aux_line) );
	  next;
      } elsif ($line =~ /^; number of runs/ ){
	  $matrix->set_attribute("num_runs", pop (@aux_line) );
	  next;
      } elsif ($line =~ /^; number of iterations/){
	  $matrix->set_attribute("num_iterations", pop (@aux_line) ) ;
	  next;
      } elsif ($line =~ /^; sequences/){
	  $matrix->set_attribute("nb_seq" , pop (@aux_line) );  
	  next;
      } elsif ($line =~ /^; total size in bp/){
	  $matrix->set_attribute("total_size_bp",  pop (@aux_line) ) ;
	  next;
      } elsif ($line =~ /^; expected motif occurrences/){
	  $matrix->set_attribute("exp_motif_occ",  pop (@aux_line) ) ;  
	  next;
      }elsif ($line =~ /^; avg.llr/){
	  $matrix->set_attribute("avg_llr", pop (@aux_line) ) ;  
	  next;
      }elsif ($line =~ /^; avg.ic/){
	  $matrix->set_attribute("avg_ic",  pop (@aux_line) ) ;  
	  next;
      }  elsif ($line =~ /^; log likelihood ratio/){
	  $matrix->set_attribute("llr",  pop (@aux_line) ) ;  
	  next;
      }  elsif ($line =~ /^; information content/){
	  $matrix->set_attribute("ic",  pop (@aux_line) ) ;  
	  next;
      } elsif ($line =~ /^; seq/ ) {
      	  $no_motif = 0;
       } elsif( ($line =~ /^;.+-i/) ){
	  $line =~ s/; //;
	  $matrix->set_attribute("command", $line );
	  next;
      }
      if ($line =~ /^; seq	strand	pos	site/ ){
	  &RSAT::message::Info("Reading sequences from infogibbs ") if ($main::verbose >= 5);
	  $read_seqs=1; 
	 # $no_motif = 0;0
	#  <STDIN>;
	  next;
      }

      next if ( ($line =~ /^;/) && ($no_motif)) ; # skip comment lines

      if ($read_seqs && ($line =~ /^;\s+\d/) ){
	  $line =~ s/\s+/\t/g; ## Replace spaces by tabulation
	  #print $line."\n";
	  local @fields2 = split /\t/, $line;
	  local $site_sequence= $fields2[4];
	  $matrix->push_attribute("sequences", $site_sequence);
	  &RSAT::message::Debug("site", $site_sequence) if ($main::verbose >= 5);
	 # <STDIN>;
	  next;
      }
      else {
	  $no_motif=1;
      }
      $line =~ s/\s+/\t/g; ## Replace spaces by tabulation
      $line =~ s/\[//g; ## Suppress [ and ] (present in the tab format of Jaspar and Pazar databases)
      $line =~ s/\]//g; ## Suppress [ and ] (present in the tab format of Jaspar and Pazar databases)
      $line =~ s/://g; ## Suppress : (present in the tab format of Uniprobe databases)
      #die $line;
      ## Create a new matrix if required
      if  ($line =~ /^\/\//) {
      	$new_matrix = 0; # tgis is to track the end of file...
	$no_motif=1;
	$read_seqs=0;
	$matrix = new RSAT::matrix();
	$matrix->set_parameter("program", "tab");
	push @matrices, $matrix;
	$current_matrix_nb++;
	$id = &clean_id($id_prefix."_".$current_matrix_nb);
	$matrix->set_attribute("AC", $id);
	$matrix->set_attribute("id", $id);
	&RSAT::message::Info("line", $l, "new matrix", $current_matrix_nb) if ($main::verbose >= 5);
	next;
      }

      if ($line =~ /^\s*(\S+)\s+/) {
	  next if ($line =~ /^;/);
	  $new_matrix = 1;
	  $matrix_found = 1; ## There is at least one matrix row in the file
	  my @fields = split /\t/, $line;


	## residue associated to the row
	my $residue = lc(shift @fields);

	## skip the | between residue and numbers
	shift @fields unless &RSAT::util::IsReal($fields[0]);
	$matrix->addIndexedRow($residue, @fields);
      }
    }
    close $in if ($file);

    ## Initialize prior as equiprobable alphabet
    if ($matrix_found) {
      if ($new_matrix == 0) {
	# eliminate empty matrix at the end
	pop(@matrices);
	$current_matrix_nb--;
      }
      &InitializeEquiPriors(@matrices);

#       foreach my $matrix (@matrices) {
# 	my @alphabet = $matrix->getAlphabet();
# 	my %tmp_prior = ();
# 	my $prior = 1/scalar(@alphabet);
# 	foreach my $residue (@alphabet) {
# 	  $tmp_prior{$residue} = $prior;
# 	  #	&RSAT::message::Debug("initial prior", $residue, $prior) if ($main::verbose >= 10);
# 	}
# 	$matrix->setPrior(%tmp_prior);
# 	if ($main::verbose >= 5) {
# 	  &RSAT::message::Debug("Read matrix with alphabet", join(":", $matrix->getAlphabet()));
# 	  &RSAT::message::Debug("Initialized prior as equiprobable", join(":", $matrix->getPrior()));
# 	  &RSAT::message::Debug("Matrix size", $matrix->nrow()." rows",  $matrix->ncol()." columns");
# 	}
#       }
    } else {
      @matrices = ();
    }

    return (@matrices);
}


# =pod

# =item _readFromOldInfoGibbsFile($file)

# Read a matrix from a result file from InfoGibbs (implementation by
# Gregory Gathy, 2007). This method is called by the method
# C<readFromFile($file, "InfoGibbs")>.

# This format was a customized version of the TRANSFAC format, developed
# for the Master thesis of Gregory Gathy. The program is not supported
# anymore, it has been replaced by Matthieu Defrance's implementation
# I<info-gibbs>.

# =cut
# sub _readFromOldInfoGibbsFile {
#   my ($file) = @_;
#   &RSAT::message::Info ("Reading matrices from InfoGibbs file", $file) if ($main::verbose >= 3);

#   ## open input stream
#   my ($in, $dir) = &main::OpenInputFile($file);
# #  my ($in) = STDIN;
# #  if ($file) {
# #    open INPUT, $file;
# #    $in = INPUT;
# #  }
#   my $current_matrix_nb = 0;
#   my @matrices = ();
#   my $matrix;
#   my $version = "";
#   my $command = "";
#   my $ncol = 0;
#   my $infogibbs_consensus = "";
#   my %select_type = ("final"=>1);

#   my %prior = ();
#   my $l = 0;
#   while (<$in>) {
#     $l++;
#     chomp();

# #    &RSAT::message::Debug($l, $_) if ($main::verbose >= 10);
#     next unless (/\S/);
#     s/\r//;
#     my $version = "";

#     ## Read the command line
#     if (/^VV\s+/) {
#       ## InfoGibbs version
#       $version = $';		# '
#       &RSAT::message::Info("InfoGibbs version", $version) if ($main::verbose >= 5);

#     } elsif (/^CM\s+/) {
#       ## InfoGibbs command
#       $command = $';		# '
#       &RSAT::message::Info("InfoGibbs command", $command) if ($main::verbose >= 5);

#     } elsif (/^PR\s+/) {
#       my $prior_error = 0;
#       ## InfoGibbs command
#       my $prior_line = $';		# '
#       my @fields = split /;\s*/, $prior_line;
#       my %new_prior = ();
#       foreach my $field (@fields) {
# 	if ($field =~ /([A-Z]):(\S+)/i) {
# 	  my $residue = lc($1);
# 	  my $prior = $2;
# 	  if (&RSAT::util::IsReal($prior)) {
# 	    $new_prior{$residue} = $prior;
# 	  } else {
# 	    &RSAT::message::Warning("InfoGibbs file", "line ".$l, "Invalid prior specification", $prior_line);
# 	    $prior_error = 1;
# 	  }
# 	} else {
# 	  &RSAT::message::Warning("InfoGibbs file", "line ".$l, "Invalid prior specification", $prior_line);
# 	    $prior_error = 1;
# 	}
#       }
#       unless ($prior_error) {
# 	%prior = %new_prior;
# #	&RSAT::message::Debug("New prior", join (" ", %prior)) if ($main::verbose >= 10);
#       }

#       ## Start a new matrix (an InfoGibbs file contains several matrices)
#     } elsif (/^AC\s+(\S+)/) {
#       my $accession = &clean_id($1);
#       &RSAT::message::Info("New matrix", $accession) if ($main::verbose >= 2);
#       $current_matrix_nb++;
#       $matrix = new RSAT::matrix();
#       $matrix->set_parameter("accession", "IG.".$accession);
#       $matrix->set_parameter("id", $accession);
#       $matrix->set_parameter("program", "InfoGibbs");
#       $matrix->set_parameter("version", $version);
#       $matrix->set_parameter("command", $command);
#       $matrix->set_parameter("matrix.nb", $current_matrix_nb);
#       if (scalar(keys(%prior)) > 0) {
# 	$matrix->setPrior(%prior);
# #	&RSAT::message::Debug("Prior", join (" ", %prior)) if ($main::verbose >= 5);
#       }
#       push @matrices, $matrix;
#       $ncol = 0;
#       $infogibbs_consensus = "";

#       &RSAT::message::Info("Parsing matrix",  $current_matrix_nb, $matrix->get_attribute("accession")) 
# 	if ($main::verbose >= 2);
#       next;

#       ## Parameters for the current matrix
#     } elsif ($matrix) {

#       ## Read prior alphabet from the matrix header (PO line)
#       ## Equiprobable alphabet.
#       if ((/^PO\s+/) || (/^P0\s+/)) {
# 	my $header = $'; #'
# 	$header = &RSAT::util::trim($header);

# 	## Alphabet is parsed from the InfoGibbs matrix header (PO row)
# 	my @alphabet = split /\s+/, $header;
# 	$matrix->setAlphabet_lc(@alphabet);

# 	## Check that prior has been specified
# 	unless ($matrix->get_attribute("prior_specified")) {
# 	  foreach my $letter (@alphabet) {
# 	    $prior{lc($letter)} = 1/scalar(@alphabet) 
# 	      unless (defined($prior{$letter}));;
# 	  }
# 	  $matrix->setPrior(%prior);
# #	  &RSAT::message::Debug("Prior", join (" ", %prior)) if ($main::verbose >= 5);
# 	}

# 	## Count column of the matrix file (row in TRANSFAC/InfoGibbs format)
#       } elsif (/^(\d+)\s+/) {
# 	my $values = $'; #'
# 	$values = &RSAT::util::trim($values);
# 	my @fields = split /\s+/, $values;
# 	my $consensus_residue= "";
# 	if ($fields[$#fields] =~ /[A-Z]/i) {
# 	  $consensus_residue = pop @fields;
# 	  $infogibbs_consensus .= $consensus_residue;
# 	}
# 	$matrix->addColumn(@fields);
# 	$ncol++;
# 	$matrix->force_attribute("ncol", $ncol);
# #	&RSAT::message::Debug("line ".$l, "adding column", $ncol, "counts", join (":", @fields))
# #	  if ($main::verbose >= 10);

# 	## Sites used to build the matrix
#       } elsif (/^BS\s+/)  {
# 	my $bs = $'; #'
# 	my ($site_sequence, $site_id) = split(/\s*;\s*/, $bs);
# 	#      my $site_sequence = $1;
# 	#      my $site_id = $2;
# 	if ($site_sequence) {
# 	  $matrix->push_attribute("sequences", $site_sequence);
# 	  if ($site_id) {
# 	    $matrix->push_attribute("site_ids", $site_id);
# 	  }
# 	}
# #	&RSAT::message::Debug("line", $l, "site", $site_sequence, $site_id, $bs) if ($main::verbose >= 10);

# 	## Information content computed by InfoGibbs
#       } elsif (/^IC\s+/) {
# 	$matrix->set_parameter("IC", $'); #'

# 	## Consensus score computed by InfoGibbs
#       } elsif (/^CS\s+/) {
# 	$matrix->set_parameter("CS", $'); #'

# 	## Log-Likelihood computed by InfoGibbs
#       } elsif (/^LL\s+/) {
# 	$matrix->set_parameter("LL", $'); #'

# 	## TRANSFAC matrix parameters (not used by InfoGibbs, but maintained for compatibility)
#       } elsif (/^XX/) {
# 	## field separator
# 	next;

#       } elsif (/^ID\s+/) {
# 	$matrix->force_attribute("name", $'); #'
# 	$matrix->force_attribute("id", $'); #'

#       } elsif (/^BF\s+/) {
# 	$matrix->set_parameter("binding_factor", $'); #'

# 	#    } elsif (/^SD\s+/) {
# 	#      $matrix->set_parameter("short_foactor_description", $'); #'

#       } elsif ((/^BA\s+/)   && ($matrix)){
# 	$matrix->set_parameter("statistical_basis", $'); #'

#       } elsif ((/^DE\s+/)   && ($matrix)){
# 	$matrix->set_parameter("description", $'); #'

# 	## Matrix type
#       } elsif ((/^TY\s+/)   && ($matrix)){
# 	$matrix->set_parameter("type", $'); #'

#       } elsif (/^\/\//) {
# 	if ($matrix) {
# 	  $matrix->set_parameter("infogibbs_consensus", $infogibbs_consensus);
# 	}

# 	## Unknown field
#       } elsif (/^(\S\S)\s+(.*)/) {
# 	my $field = $1;
# 	my $value = $2;
# 	&RSAT::message::Warning("Unknown field, not parsed", "line ".$l, $field, $value) if ($main::verbose >= 5);

#       } else {
# 	&RSAT::message::Warning("skipped invalid row", "line ".$l, $_);
#       }
#     }

#   }
#   close $in if ($file);

#   my @selected_matrices = ();
#   foreach my $matrix (@matrices) {
#     my $type = $matrix->get_attribute("type");
#     if (($type) && ($select_type{$type})) {
#       push @selected_matrices, $matrix;
#     }
#   }

#   return @selected_matrices;

# }


=pod

=item _readFromAlignACEFile($file)

Read a matrix from an AlignACE file. This method is called by the method 
C<readFromFile($file, "AlignACE")>.

=cut
sub _readFromAlignACEFile {
  my ($file) = @_;
  my @matrices = ();

  my ($in, $dir) = &main::OpenInputFile($file);
  my $l = 0;
  my $matrix;
  my $matrix_nb;
  my $site_nb;
  my $expect = "";
  my $gcback = 0.5;
  my $minpass = "";
  my $seed = "";
  my $numcols = "";
  my $undersample = "";
  my $oversample = "";
  my %seq_id = ();
  while (<$in>) {
    $l++;
    chomp();
    next unless (/\S/); ## Skip empty lines
    if (/^Motif (\d+)/) {
      $matrix_nb = $1;
      $site_nb = 0;
      $matrix = new RSAT::matrix();
      $matrix->set_parameter("id", "Motif_".$matrix_nb);
      $matrix->set_parameter("accession", "Motif_".$matrix_nb);
      $matrix->set_parameter("matrix.nb", $matrix_nb);
      $matrix->set_parameter("program", "AlignACE");
      $matrix->set_parameter("command", $AlignACE_command);
      $matrix->set_parameter("program.version", $AlignACE_version);
      $matrix->set_parameter("program.release", $AlignACE_date);
      $matrix->set_parameter("seed", $seed);
      $matrix->set_parameter("alignace.minpass", $minpass);
      $matrix->set_parameter("alignace.expect", $expect);
      $matrix->set_parameter("alignace.undersample", $undersample);
      $matrix->set_parameter("alignace.oversample", $oversample);
      &RSAT::message::Info("Starting to read matrix", $matrix_nb) if ($main::verbose >= 5);

      ## Default nucletodide alphabet
#      $matrix->setAlphabet_lc("a","c","g","t");
      $matrix->set_alphabet_for_type();
      my $atback = 1-$gcback;
      $matrix->setPrior(a=>$atback/2, c=>$gcback/2,t=>$atback/2, g=>$gcback/2);
      $matrix->force_attribute("nrow",4);
      push @matrices, $matrix;
      $in_matrix = 1;
    } elsif (/^AlignACE (\d+\.\d+)\s+(\S+)/) {
      $AlignACE_version = $1;
      $AlignACE_date = $2;
    } elsif (/AlignACE/) {
      $AlignACE_command = $_;
    } elsif (/gcback = \s+(\S+)/) {
      $gcback = $1;
    } elsif (/seed\s+=\s+(\S+)/) {
      $seed = $1;
    } elsif (/expect\s+=\s+(\S+)/) {
      $expect = $1;
    } elsif (/minpass\s+=\s+(\S+)/) {
      $minpass = $1;
    } elsif (/oversample\s+=\s+(\S+)/) {
      $oversample = $1;
    } elsif (/undersample\s+=\s+(\S+)/) {
      $undersample = $1;
    } elsif (/#(\d+)\s+(\S+)/) {
      $seq_id[$1] = $2;
    } elsif ($in_matrix) {
      if (/MAP Score:\s(\S+)/i) {
	$matrix->set_parameter("MAP", $1);
      } elsif (/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
	$site_nb++;
	my $site_seq = $1;
	my $seq_nb = $2;
	my $site_pos = $3;
	my $site_strand = $4;
	my $site_id = join ("_", "mtx".$matrix_nb, "site".$site_nb, "seq".$seq_nb, $seq_id[$seq_nb], $site_pos, $site_strand);
	$matrix->add_site(lc($site_seq), id=>$site_id, score=>1);
      }
    }
  }
  close $in if ($file);

  ## Compute matrix consensus
  foreach my $matrix (@matrices) {
      ## Replace undefined values by 0
    $matrix->treat_null_values();
    $matrix->calcConsensus();
    my $consensus = $matrix->get_attribute("consensus.IUPAC");
#    my $ac = $matrix->get_attribute("accession");
    $matrix->force_attribute("id", $consensus);
    $matrix->force_attribute("name", $consensus);
  }

  return @matrices;
}


=pod

=item _readFromGibbsFile($file)

Read a matrix from a gibbs file. This method is called by the method 
C<readFromFile($file, "gibbs")>.

By default, the matrix is parsed from the sites.

The variable $parse_model allows to parse the matrix model exported by
the gibbs sampler. I prefer to avoid this because this model is
exported in percentages, and with some fuzzy rounding.

=cut
sub _readFromGibbsFile {
    my ($file) = @_;

    my $parse_model = 0; ## boolean: indicates whether or not to parse matrices from motif models
    my $initial_matrices = 0;  ## boolean: indicates whether or not to export the initial matrices resulting from the optimization
    my $final_matrices = 1; ## boolean: indicates whether or not to export the final matrices

    ## open input stream
   my ($in, $dir) = &main::OpenInputFile($file);
#   my $in = STDIN;
#     if ($file) {
#	open INPUT, $file;
#	$in = INPUT;
#    }
    $in_matrix = 0;

    my $matrix;
    my @matrices = ();
    my @matrix = ();
    my @alphabet = ();
    my $ncol = 0;
    my $nrow = 0;
    my $m = 0;
    my $gibbs_command = "";
    my $seed;

    while (<$in>) {
      $l++;

      chomp();

      ## Suppress DOS-type newline characters
      s/\r//;

      ## Empty rows indicate the end of a matrix
      unless (/\S/) {
	if ($in_matrix) {
	  $in_matrix = 0;
	}
	next;
      }

      if (/^gibbs /) {
	$gibbs_command = $_;
#	&RSAT::message::Debug("line ".$l, "gibbs command", $gibbs_command) if ($main::verbose >= 10);

      } elsif (/seed: (\S+)/) {
	$seed = $1;
#	&RSAT::message::Debug("line ".$l, "seed", $seed) if ($main::verbose >= 10);

      } elsif (/^\s*(\d+)\-(\d+)\s+(\d+)\s+([a-z]*)\s+([A-Z]+)\s+([a-z]+)\s+(\d+)\s*(\S*)/) {

	if ($in_matrix) {
#	  &RSAT::message::Debug("line ".$l, "Parsing one matrix row") if ($main::verbose >= 10);
	} else {
#	  &RSAT::message::Debug("line ".$l, "Starting to read a matrix") if ($main::verbose >= 10);
	  $matrix = new RSAT::matrix();
	  $matrix->set_parameter("program", "gibbs");
	  $matrix->set_parameter("command", $gibbs_command);
	  $matrix->set_parameter("seed", $seed);
	  push @matrices, $matrix;
	  $in_matrix = 1;
	  # default nucletodide alphabet
#	  $matrix->setAlphabet_lc("a","c","g","t");
	  $matrix->set_alphabet_for_type();
#	  $matrix->force_attribute("nrow",4);
	}

	my $seq_nb = $1;
	my $site_nb=$2;
	my $start=$3;
	my $left_flank=$4;
	my $site_seq=$5;
	my $right_flank=$4;
	my $end=$7;
	my $score=$8; $score =~ s/\(//; $score =~ s/\)//;
	my $site_id;
	if ($score eq "") {
	  $matrix->set_parameter("gibbs.type", "initial");
	  $site_id = join ("_", $seq_nb, $site_nb, $start, $end);
	} else {
	  $matrix->set_parameter("gibbs.type", "final");
	  $site_id = join ("_", $seq_nb, $site_nb, $start, $end, $score);
	}

	$matrix->add_site(lc($site_seq), id=>$site_id, score=>1);

      } elsif ((/^Motif model/) && ($parse_model)) {
	&RSAT::message::Debug("line ".$l, "Creating a new model matrix") if ($main::verbose >= 5);
#      if (/^\s*MOTIF\s+(\S+)/) {
	$matrix = new RSAT::matrix();
	$matrix->set_parameter("program", "gibbs");
	$matrix->set_parameter("command", $gibbs_command);
	$matrix->set_parameter("seed", $seed);
	push @matrices, $matrix;
	&RSAT::message::Debug("Starting to read a motif") if ($main::verbose >= 5);
	$in_matrix = 1;
	# default nucletodide alphabet
#	$matrix->setAlphabet_lc("a","c","g","t");
	$matrix->set_alphabet_for_type();
	next;

      } elsif ((/model map = (\S+); betaprior map = (\S+)/) && ($in_matrix)) {
	$matrix->set_parameter("gibbs.model.map", $1);
	$matrix->set_parameter("gibbs.betaprior.map", $2);
#	&RSAT::message::Warning("gibbs matrix", $matrix,
#				"model map", $matrix->get_attribute("gibbs.model.map"),
#				"betaprior map", $matrix->get_attribute("gibbs.betaprior.map"));

      } elsif ((/^\s*MAP = (\S+)/) && ($matrix)) {
	$matrix->set_parameter("MAP", $1);

      } elsif ((/^\s*NetMAP = (\S+)/) && ($matrix)) {
	$matrix->set_parameter("NetMAP", $1);

      } elsif ((/^\s*sites: MAP = (\S+)/) && ($matrix)) {
	$matrix->set_parameter("sites.MAP", $1);

      } elsif ((/^\s*Initial MAP = (\S+)/) && ($matrix)) {
	$matrix->set_parameter("initial.MAP", $1);

      } elsif (($in_matrix) && ($parse_model)) {
	if (/^\s*POS/) {
	  ## Header line, indicating the alphabet
	  s/\r//;
	  chomp;
	  @header = split " +";
	  @alphabet = @header[1..$#header-1];
#	  $matrix->setAlphabet_lc(@alphabet);
	  $matrix->setAlphabet(@alphabet); ## Here we use the alphabet defined in the file rather than the user-specified type
#	  $matrix->set_alphabet_for_type();
	  &RSAT::message::Debug("Alphabet", join(":", @alphabet)) if ($main::verbose >= 5);

	} elsif (/^\s*\d+\s+/) {
	  ## Add a column to the matrix (gibbs rows correspond to our columns)
	  s/\r//;
	  chomp;
	  s/^\s+//;
	  @fields = split " +";
#	  &RSAT::message::Debug("col", $ncol, "fields", scalar(@fields), @fields) if ($main::verbose >= 10);
	  @values = @fields[1..$#header-1];
	  $nrow = scalar(@values);
	  foreach my $v (0..$#values) {
	    $values[$v] =~ s/^\.$/0/;
	    $matrix[$ncol][$v] = $values[$v];
	  }
#	  &RSAT::message::Debug("col", $ncol, "values", scalar(@values), @values) if ($main::verbose >= 10);
	  $ncol++;

	} elsif (/site/) {
	  $in_matrix = 0;
	  $matrix->force_attribute("nrow", $nrow);
	  $matrix->force_attribute("ncol", $ncol);
	  $matrix->setMatrix ($nrow, $ncol, @matrix);
	  @matrix = ();
	  $nrow = 0;
	  $ncol = 0;
	  next;
	}
      }
    }
    close $in if ($file);


    ## Delete the pre-matrices (first half of the matrices)
    unless ($initial_matrices) {
      for my $m (1..(scalar(@matrices)/2)) {
	shift @matrices;
      }
    }
    unless ($final_matrices) {
      for my $m (1..(scalar(@matrices)/2)) {
	pop @matrices;
      }
    }

    return @matrices;
}


=pod

=item _readFromConsensusFile($file)

Read a matrix from a consensus file. This method is called by the
method C<readFromFile($file, "consensus")>.

=cut
sub _readFromConsensusFile {
  my ($file) = @_;
  &RSAT::message::Info ("Reading matrix from consensus file", $file) if ($main::verbose >= 3);
  #    ($in, $dir) = &main::OpenInputFile($file);

  ## open input stream
#  my $in = STDIN;
  my ($in, $dir) = &main::OpenInputFile($file);
  my $current_matrix_nb = 0;
  my @matrices = ();
  my $matrix;
  my $command = "";

  my $prefix = "matrix";
  if ($file) {
    $prefix = &RSAT::util::ShortFileName($file);
  }

  my %prior = ();
  my $l = 0;
  my $final_cycle = 0;
  while (<$in>) {
    $l++;
#    &RSAT::message::Debug("line", $l, $_) if ($main::verbose >= 5);
    next unless (/\S/);
    s/\r//;
    chomp();


    ## Read the command line
    if (/COMMAND LINE: /) {
      $command = $';		# '

    } elsif (/THE LIST OF MATRICES FROM FINAL CYCLE/) {
      $final_cycle = 1;

      ## Start a new matrix (one consensus file contains several matrices)
    } elsif ((/MATRIX\s(\d+)/) && ($final_cycle)) {
      $current_matrix_nb = $1;
      $matrix = new RSAT::matrix();
      push @matrices, $matrix;
      $matrix->set_parameter("program", "consensus");
      $matrix->set_parameter("id", $prefix.$current_matrix_nb);
      $matrix->set_parameter("matrix.nb", $current_matrix_nb);
      $matrix->set_parameter("command", $command);
      $matrix->setPrior(%prior);
      next;

      ## Read prior frequency for one residue in the consensus header
    } elsif (/letter\s+\d+:\s+(\S+).+prior frequency =\s+(\S+)/) {
      my $letter = lc($1);
      my $prior = $2;
      &RSAT::message::Info ("Prior from consensus file", $letter, $prior) if ($main::verbose >= 5);
      $prior{$letter} = $prior;

    } elsif ($current_matrix_nb >= 1) {

      ## Matrix content (counts) for one residue
      if (/^\s*(\S+)\s+\|/) {
	my @fields = split /\s+/, $_;
	## residue associated to the row
	my $residue = lc(shift @fields);
	$residue =~ s/\s*\|//;

	## skip the | between residue and numbers
	shift @fields unless &RSAT::util::IsReal($fields[0]);
	$matrix->addIndexedRow($residue, @fields);


#	&RSAT::message::Debug("&_readFromConsensusFile", $residue, "alphabet", join(":", $matrix->getAlphabet()), join ", ", @fields)
#	  if ($main::verbose >= 10);

	## Sites used to build the matrix
      } elsif (/(\d+)\|(\d+)\s*\:\s*(-){0,1}(\d+)\/(\d+)\s+(\S+)/) {
	my $site_nb = $1;
	my $site_cycle = $2;
	my $site_strand = $3;
	my $site_seq_nb = $4;
	my $site_pos = $5;
	my $site_sequence = $6;
	$matrix->push_attribute("sequences", $site_sequence);
	&RSAT::message::Debug("line", $l, "site", $site_sequence) if ($main::verbose >= 10);

	## Other matrix parameters
      } elsif (/number of sequences = (\d+)/) {
	$matrix->set_parameter("sites", $1); 
      } elsif (/unadjusted information = (\S+)/) {
	$matrix->set_parameter("cons.unadjusted.information", $1); 
      } elsif (/sample size adjusted information = (\S+)/) {
	$matrix->set_parameter("cons.adjusted.information", $1); 
      } elsif (/ln\(p\-value\) = (\S+)   p\-value = (\S+)/) {
	$matrix->set_parameter("cons.ln.Pval", $1); 
	$matrix->set_parameter("cons.Pval", $2); 
      } elsif (/ln\(e\-value\) = (\S+)   e\-value = (\S+)/) {
	$matrix->set_parameter("cons.ln.Eval", $1);
	$matrix->set_parameter("cons.Eval", $2); 
      } elsif (/ln\(expected frequency\) = (\S+)   expected frequency = (\S+)/) {
	$matrix->set_parameter("cons.ln.exp", $1); 
	$matrix->set_parameter("cons.exp", $2); 
      }
    }
  }

  ##Check if there was at least one matrix obtained after final cycle
  unless ($final_cycle) {
    &RSAT::message::Warning("This file does not contain the \"FINAL CYCLE\" header") if ($main::verbose >= 2);
  } 
  unless (scalar(@matrices) > 0) {
    &RSAT::message::Warning("This file does not contain any final cycle matrix") if ($main::verbose >= 2);
  }
  close $in if ($file);

  &RSAT::message::Debug("matrices read", scalar(@matrices)) if ($main::verbose >= 3);

  &RSAT::message::Debug("matrices after sorting", scalar(@matrices)) if ($main::verbose >= 3);

  return @matrices;
}


=pod

=item _readFromAssemblyFile($file)

Read a matrix from the output file of pattern-assembly. This method is
called by the method C<readFromFile($file, "assembly")>.

=cut
sub _readFromAssemblyFile {
  my ($file) = @_;
  &RSAT::message::Info ("Reading matrix from pattern-assembly file", $file) if ($main::verbose >= 3);

  ## open input stream
  my ($in, $dir) = &main::OpenInputFile($file);

  my $current_matrix_nb = 0;
  my @matrices = ();
  my $matrix;
  my $command = "";

  my $l = 0;
  while (my $line = <$in>) {
    $l++;
    next unless ($line =~  /\S/);
    chomp($line);
    &RSAT::message::Debug(" _readFromAssemblyFile", $l, $line) if ($main::verbose >= 5);

    ## Read the command line
    if ($line =~ /; *assembly # (\d+)*\s+(.*)seed:\s+(\S+)/) {
      $current_matrix_nb = $1;
      my $seed = $3;
      my $cluster = 0;
      if ($2 =~ /cluster\s+#\s+(\d+)/) {
	$cluster = $1;
      }
      $matrix = new RSAT::matrix();
      $matrix->set_parameter("program", "pattern-assembly");
      $matrix->set_parameter("matrix.nb", $current_matrix_nb);
      my $matrix_id;
      if ($cluster) {
	$matrix_id = "assembly_".$current_matrix_nb."_cluster_".$cluster;
      } else {
	$matrix_id =  "assembly_".$current_matrix_nb;
      }
      $matrix->set_parameter("id", $matrix_id);
      my $matrix_name = $seed;
      $matrix->set_parameter("name", $matrix_name);
      push @matrices, $matrix;
#      $matrix->setAlphabet_lc("A","C","G","T");
      $matrix->set_alphabet_for_type();
#      $matrix->force_attribute("nrow", 4);
      $matrix->set_parameter("asmb.seed", $seed);
      &RSAT::message::Debug("New matrix from assembly", $current_matrix_nb."/".scalar(@matrices), "seed", $seed) if ($main::verbose >= 5);

    } elsif ($line =~ /^(\S+)\t(\S+)\s+(\S+)\s+isol/) {
      $current_matrix_nb++;
      my $pattern = $1; 
      my $pattern_rc = $2; 
      my $score = $3;
      $matrix = _from_isolated($pattern, $pattern_rc, $score, @matrices);
      push @matrices, $matrix;

    } elsif ($line =~ /^(\S+)\s+(\S+)\s+isol/) {
      $current_matrix_nb++;
      my $pattern = $1; 
      my $score = $2;
      $matrix = _from_isolated($pattern, "", $score, @matrices);
      push @matrices, $matrix;


      ## Consensus from a 2-strand assembly
    } elsif ($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+best consensus/) {
      $matrix->force_attribute("name", $1);
      $matrix->set_attribute("consensus.assembly", $1);
      $matrix->set_attribute("consensus.assembly.rc", $2);
      $matrix->set_attribute("assembly.top.score", $3);
#      &RSAT::message::Debug("Consensus for matrix", $current_matrix_nb, $1) if ($main::verbose >= 5);

      ## Consensus from a 1-strand assembly
    } elsif ($line =~ /^(\S+)\s+(\S+)\s+best consensus/) {
      $matrix->set_attribute("consensus.assembly", $1);
      $matrix->set_attribute("assembly.top.score", $3);
#      &RSAT::message::Debug("Consensus for matrix", $current_matrix_nb, $1) if ($main::verbose >= 5);

    } elsif ($line =~ /^;/) {
      next;

      ## New site from a two-strands assembly
    } elsif ($line =~ /^(\S+)\t(\S+)\s+(\S+)/) {
      my $pattern = $1; 
      my $pattern_rc = $2;
      my $score = $3;
      my $pattern_id = $pattern."|";
      $pattern =~ s/\./n/g;
      $pattern_rc =~ s/\./n/g;
#      &RSAT::message::Debug("ASSEMBLY LINE", $l, $pattern, $pattern_rc, $score) if ($main::verbose >= 10);
      $matrix->add_site(lc($pattern), score=>$score, id=>$pattern_id, max_score=>1);

      ## New site from a single-strand assembly
    } elsif ($line =~ /^(\S+)\s+(\S+)/) {
      my $pattern = $1; 
      my $score = $2;
      $pattern =~ s/\./n/g;
      &RSAT::message::Debug("ASSEMBLY LINE", $l, $pattern, $pattern_rc, $score) if ($main::verbose >= 5);
      $matrix->add_site(lc($pattern), score=>$score, id=>$pattern, max_score=>1);

    } else {
      &RSAT::message::Warning("&RSAT::Matrixreader::_readFromAssemblyFile", "line", $l, "not parsed", $_) if ($main::verbose >= 5);
    }
  }
  close $in if ($file);

  return @matrices;
}


=pod

=item _from_isolated($pattern, $pattern_rc, $score, @matrices)

Create a matrix from an isolated pattern string.

=cut
sub _from_isolated {
  my ($pattern, $pattern_rc, $score, @matrices) = @_;
  unless (&RSAT::util::IsReal($score)) {
    &RSAT::message::Warning($score, "is not a valid score value for pattern", $pattern) if ($main::verbose >= 5);
    $score = 1;
  }
  my $pattern_id = $pattern;
  $pattern =~ s/\./n/g;
  if ($pattern_rc) {
    $pattern_rc =~ s/\./n/g;
    $pattern_id = $pattern."|".$pattern_rc;
  }
  $matrix = new RSAT::matrix();
  $matrix->set_parameter("program", "pattern-assembly");
  $matrix->set_parameter("matrix.nb", $current_matrix_nb);
  # $matrix->setAlphabet_lc("A","C","G","T");
  $matrix->set_alphabet_for_type();
#  $matrix->force_attribute("nrow", 4);
  $matrix->set_parameter("asmb.seed", $pattern);
  $matrix->set_attribute("asmb.consensus", $pattern);
  $matrix->set_attribute("asmb.consensus.rc", $pattern_rc);
  $matrix->set_attribute("asmb..top.score", $score);
  $matrix->add_site(lc($pattern), score=>$score, id=>$pattern_id,, max_score=>1);
  &RSAT::message::Debug("New matrix from isolated pattern", $current_matrix_nb."/".scalar(@matrices), "seed", $seed) if ($main::verbose >= 5);
  return $matrix;
}


=pod

=item _readFromTabFile($file)

Read a matrix from a tab-delimited file. This method is called by the
method C<readFromFile($file, "tab")>.

=cut
sub _readFromTabFile {
    my ($file, %args) = @_;
    &RSAT::message::Info(join("\t", "Reading matrix from tab file\t",$file)) if ($main::verbose >= 3);

    ## Open input stream
    my ($in, $dir) = &main::OpenInputFile($file);

    ## Read header
    if ($args{header}) {
	$header = <$in>;
	$header =~ s/\r//;
	chomp ($header);
	$header =~ s/\s+/\t/g;
	@header = split "\t", $header;
	$matrix->push_attribute("header", @header);
    }

    ################################################################
    ## Initialize the matrix list
    my @matrices = ();
    my $matrix = new RSAT::matrix();
    $matrix->set_parameter("program", "tab");
    $matrix->set_parameter("matrix.nb", $current_matrix_nb);
    push @matrices, $matrix;
    my $current_matrix_nb = 1;
    #    my $id = $file."_".$current_matrix_nb;
    my $no_path_file=$file;
    $no_path_file=~s/.+\/+//;
    my $id_prefix = $no_path_file|| "matrix";
    my $id = &clean_id($id_prefix."_".$current_matrix_nb);
    $matrix->set_attribute("AC", $id);
    $matrix->set_attribute("id", $id);
    my $l = 0;
    my $matrix_found = 0;
    my $new_matrix = 0;
    while ($line = <$in>) {
      $l++;
      next unless ($line =~ /\S/); ## Skip empty lines
      chomp($line); ## Suppress newline
      $line =~ s/\r//; ## Suppress carriage return
      $line =~ s/(^.)\|/$1\t\|/; ## Add tab after residue if missing
      $line =~ s/\s+/\t/g; ## Replace spaces by tabulation
      next if ($line =~ /^;/) ; # skip comment lines
      $line =~ s/\[//g; ## Suppress [ and ] (present in the tab format of Jaspar and Pazar databases)
      $line =~ s/\]//g; ## Suppress [ and ] (present in the tab format of Jaspar and Pazar databases)
      $line =~ s/://g; ## Suppress : (present in the tab format of Uniprobe databases)

      ## Create a new matrix if required
      if  ($line =~ /\/\//) {
      	$new_matrix = 0; # tgis is to track the end of file...
	$matrix = new RSAT::matrix();
	$matrix->set_parameter("program", "tab");
	push @matrices, $matrix;
	$current_matrix_nb++;
	$id = &clean_id($id_prefix."_".$current_matrix_nb);
	$matrix->set_attribute("AC", $id);
	$matrix->set_attribute("id", $id);
	&RSAT::message::Info("line", $l, "new matrix", $current_matrix_nb) if ($main::verbose >= 5);
	next;
      }

      ## Detect error in the input format specification (user selected
      ## tab for TRANSFAC-formatted matrix)
      if (($line =~ /^po/i) 
	  || ($line =~ /^ac/i)
	  || ($line =~ /^id/i)) {
	  &RSAT::error::FatalError("Matrix does not seem to be in tab-delimited format. Seems to be a transfac-formatted matrix.");
      }

      if ($line =~ /^\s*(\S+)\s+/) {
	$new_matrix = 1;
	$matrix_found = 1; ## There is at least one matrix row in the file
	my @fields = split /\t/, $line;

	## residue associated to the row
	my $residue = lc(shift @fields);

	## skip the | between residue and numbers
	shift @fields unless &RSAT::util::IsReal($fields[0]);
	$matrix->addIndexedRow($residue, @fields);
      }
    }
    close $in if ($file);

    ## Initialize prior as equiprobable alphabet
    if ($matrix_found) {
    	if ($new_matrix == 0) {
	    # eliminate empty matrix at the end
	    pop(@matrices);
	    $current_matrix_nb--;
	}
	&InitializeEquiPriors(@matrices);
    } else {
      @matrices = ();
    }

    return (@matrices);
}


=pod

=item _readFromClusterBusterFile($file)

Read a matrix from a file in ClusterBuster format (files with
extension .cb). This method is called by the method
C<readFromFile($file, "cluster-buster")>.

=cut
sub _readFromClusterBusterFile {
    my ($file, %args) = @_;
    &RSAT::message::Info(join("\t", "Reading matrix from ClusterBuster file\t",$file)) if ($main::verbose >= 3);


    ## open input stream
    my ($in, $dir) = &main::OpenInputFile($file);

    ## Initialize the matrix list
    my @matrices = ();
    my @alphabet = qw(a c g t); ## cluster-buster does not explicitly indicate the alphabet, it is always supposed to represent DNA matrices
    my $matrix;
    my $current_matrix_nb = 1;
    my $l = 0;
    my $ncol = 0;
    while ($line = <$in>) {
      $l++;
      next unless ($line =~ /\S/); ## Skip empty lines
      chomp($line); ## Suppress newline
      $line =~ s/\r//; ## Suppress carriage return
      $line =~ s/\s+/\t/g; ## Replace spaces by tabulation
      next if ($line =~ /^;/) ; # skip comment lines
      #	&RSAT::message::Debug("line", $l, $line) if ($main::verbose >= 10);
      ## Create a new matrix if required
      if  ($line =~ /^\>(\S*)/) {
	my $name = $1;
	if ($line =~ /\/name=(\S*)/) { $name = $1;}
	$matrix = new RSAT::matrix();
	$matrix->set_parameter("program", "clusterbuster");
	$ncol = 0;
	if ($name) {
	  $matrix->set_attribute("name", $name);
	  $matrix->set_attribute("AC", $name);
	  $matrix->set_attribute("accession", $name);
	}
	# my @alphabet = qw(a c g t);
	# $matrix->setAlphabet_lc(@alphabet);
	$matrix->set_alphabet_for_type();
#	$matrix->force_attribute("nrow", 4);
	push @matrices, $matrix;
	$current_matrix_nb++;
	&RSAT::message::Info("line", $l, "new matrix", $current_matrix_nb, $name) if ($main::verbose >= 5);
	next;
      }

      if ($line =~ /^\s*(\S+)\s+/) {
	$line = &main::trim($line);
	my @fields = split /\t/, $line;
#	&RSAT::message::Info("line", $l, "adding column", join(";", @fields)) if ($main::verbose >= 10);
	$matrix->addColumn(@fields);
	$ncol++;
	$matrix->force_attribute("ncol", $ncol);
      }
    }
    close $in if ($file);

    ## Initialize prior as equiprobable alphabet
    &InitializeEquiPriors(@matrices);

    return (@matrices);
}
=pod

=item B<_readFromEncodeFile($file)>

Read a matrix from a file in Encode format. This method is called by the method
C<readFromFile($file, "encode")>.

Source: 
  http://compbio.mit.edu/encode-motifs/

Reference:

    Pouya Kheradpour and Manolis Kellis (2013). Systematic discovery
    and characterization of regulatory motifs in ENCODE TF binding
    experiments Nucleic Acids Research, 2013 December 13,.
    doi:10.1093/nar/gkt1249 

=cut
sub _readFromEncodeFile {
    my ($file) = @_;
    &RSAT::message::Info(join("\t", "Reading matrix from Encode file\t",$file)) if ($main::verbose >= 3);


    ## open input stream
    my ($in, $dir) = &main::OpenInputFile($file);

    ## Initialize the matrix list
    my @matrices = ();
    my @alphabet = qw(a c g t); ## encode does not explicitly indicate the alphabet, it is always supposed to represent DNA matrices
    my $matrix;
    my $current_matrix_nb = 1;
    my $l = 0;
    my $ncol = 0;
    while ($line = <$in>) {
      $l++;
      next unless ($line =~ /\S/); ## Skip empty lines
      chomp($line); ## Suppress newline
      $line =~ s/\r//; ## Suppress carriage return
      $line =~ s/\s+/\t/g; ## Replace spaces by tabulation
      next if ($line =~ /^;/) ; # skip comment lines
      #	&RSAT::message::Debug("line", $l, $line) if ($main::verbose >= 10);
      ## Create a new matrix if required
      if  ($line =~ /^\>(\S*)/) {
	  ## First line is compossed of two elements, first the TF ID, then separated by space
	  ## a description that some times includes the pograms used to discover the motif
	  my ($name,$comment)=split(" ",$line);
	  $name=~s/^>//;
	  $comment=~s/#/_/g;
	  #print join ("+++" ,$name,$comment);
	  #die "BOOM";
	  #if ($line =~ /\/name=(\S*)/) { $name = $1;}
	  $matrix = new RSAT::matrix();
	  $matrix->set_parameter("program", "encode");
	  $ncol = 0;
	  if ($name) {
	      $matrix->set_attribute("name", $name);
	      $matrix->set_attribute("AC", $name);
	      $matrix->set_attribute("accession", $name);
	  }
	  if ($comment){
	      $matrix->set_parameter("description", $comment);
	  }
	  # my @alphabet = qw(a c g t);
	  # $matrix->setAlphabet_lc(@alphabet);
	  $matrix->set_alphabet_for_type();
	  #	$matrix->force_attribute("nrow", 4);
	  push @matrices, $matrix;
	  $current_matrix_nb++;
	  &RSAT::message::Info("line", $l, "new matrix", $current_matrix_nb, $name) if ($main::verbose >= 5);
	  next;
      }

      if ($line =~ /^\s*(\S+)\s+/) {
	$line = &main::trim($line);
	my @fields = split /\t/, $line;
	## First character in the column is the consesus
	shift (@fields);
#	&RSAT::message::Info("line", $l, "adding column", join(";", @fields)) if ($main::verbose >= 10);
	$matrix->addColumn(@fields);
	$ncol++;
	$matrix->force_attribute("ncol", $ncol);
      }
    }
    close $in if ($file);

    ## Initialize prior as equiprobable alphabet
    &InitializeEquiPriors(@matrices);

    return (@matrices);
}


=pod

=item _readFromUniprobeFile($file)

Read a matrix from a file in Uniprobe format. This method is called by
the method C<readFromFile($file, "uniprobe")>.

=cut
sub _readFromUniprobeFile {
    my ($file, %args) = @_;
    &RSAT::message::Info(join("\t", "Reading matrix from Uniprobe file\t",$file)) if ($main::verbose >= 3);


    ## open input stream
    my ($in, $dir) = &main::OpenInputFile($file);

    ## Initialize the matrix list
    my @matrices = ();
    my $matrix;
    my $current_matrix_nb = 0;
    my $l = 0;
    my $ncol = 0;
    my $multiply = 100; ## Uniprobe matrices are provided in relative frequencies

    my $in_matrix = 0;
#    my @alphabet = qw(a c g t); ## The alphabet will be read from the matrix

    while ($line = <$in>) {
	$l++;
#      next unless ($line =~ /\S/); ## Skip empty lines
	chomp($line); ## Suppress newline character
	$line =~ s/\r//; ## Suppress Windows-specific carriage return
	$line =~ s/\s+/\t/g; ## Replace spaces by tabulation

	## Read frequencies for one residue
	if ($line =~ /^(\S)\:\s*/) {
	    if ($in_matrix) {
		my $residue = lc($1);
		my @fields = split(/\s+/, $line);
		shift(@fields);
		$matrix->addIndexedRow($residue, @fields);
		&RSAT::message::Debug("Added row", $residue, join ":", @fields), if ($main::verbose >= 5);
		next;
	    } else {
		&RSAT::error::FatalError("This file does not seem to be in Uniprobe format");
	    }

	    ## New matrix identifier
	} elsif  (($line =~ /Protein:\s+(\S+)/) ||
		  ((!$in_matrix) && ($line =~ /^(\S+)/))) {
	    
	    ## Stop reading if the number of matrices to read has been restricted
	    if (defined($args{top})) {
		if (scalar(@matrices) >= $args{top}) {
		    &RSAT::message::Warning("Stop reading after",
					    $current_matrix_nb, scalar(@matrices),
					    "top matrices (-top $args{top}), but the file seems to contain more.") if ($main::verbose >= 2);
		    last;
		}
	    }
	    
	    ## Create a new matrix
	    $current_matrix_nb++;
	    $in_matrix = 1;
	    $matrix = new RSAT::matrix();
	    $matrix->set_parameter("program", "uniprobe");
	    $ncol = 0;
	    
	    ## Protein name becomes matrix name + ID
	    my $name = $1;
	    if ($line =~ /Protein:\s+(\S+)/i) {
		$name = $1;
	    }
	    $matrix->force_attribute("id", $name);
	    $matrix->set_attribute("name", $name);
	    $matrix->set_attribute("AC", $name);
	    $matrix->set_attribute("accession", $name);
	    
	    ## Seed k-mer
	    my $seed_kmer = "";
	    if ($line =~ /k-mer:\s+(\S+)/i) {
		$seed_kmer = $1;
		$matrix->set_attribute("seed_kmer", $seed_kmer);
	    }
	    ## Enrichment score
	    my $enrich_score = "";
	    if ($line =~ /Score:\s+(\S+)/i) {
		$enrich_score = $1;
		$matrix->set_attribute("enrich_score", $enrich_score);
	    }
	    
	    &RSAT::message::Info("New Uniprobe matrix", $name, $seed_kmer, $enrich_score) if ($main::verbose >= 3);
#	$matrix->setAlphabet_lc(@alphabet);
#	$matrix->force_attribute("nrow", 4);
	    push @matrices, $matrix;
	    &RSAT::message::Info("line ".$l, "matrix ".$current_matrix_nb, $name) if ($main::verbose >= 4);
	    next;

	    ## Blank rows separate multiple matrices
	} elsif ($line !~ /\S/) {
	    &RSAT::message::Debug("Separator line between matrices") if ($main::verbose >= 5);
#	print($matrix->to_tab());
	    $in_matrix = 0;
	}

    }
    close $in if ($file);


    ## Multiply the matrix frequencies to obtain numbers > 1
#    &RSAT::message::Info("Multiplying Uniprobe frequencies by ".$multiply) if ($main::verbose >= 10);
#    foreach my $marix (@matrices) {
#      $matrix->multiply($multiply);
#    }

    ## Initialize prior as equiprobable alphabet
    &InitializeEquiPriors(@matrices);

    return (@matrices);
}


=pod

=item _readFromJasparFile($file)

Read a matrix from a file in JASPAR format (files with extension
.jaspar). JASPAR is a public database of transcription factor binding
sites and motifs (http://jaspar.cgb.ki.se/).

This method is called by the method C<readFromFile($file, "jaspar")>.

=cut
sub _readFromJasparFile {
    my ($file, %args) = @_;
    &RSAT::message::Info(join("\t", "Reading matrix from JASPAR file\t",$file)) if ($main::verbose >= 3);


    ## open input stream
    my ($in, $dir) = &main::OpenInputFile($file);
#    if ($file) {
#	open INPUT, $file;
#	$in = INPUT;
#    }

    ## Special treatment for the alphabet: sometimes indicated in the first column, sometimes not.
    my @temp_alphabet = qw(a c g t);

    ## Initialize the matrix LIST
    local @matrices = (); ## updated in subroutine
    local $matrix; ## updated in subroutine
    local $current_matrix_nb = 1; ## updated in subroutine
    my $l = 0;
    my $ncol = 0;
    while ($line = <$in>) {
      $l++;
      next unless ($line =~ /\S/); ## Skip empty lines
      chomp($line); ## Suppress newline
      $line =~ s/\r//; ## Suppress carriage return
      $line =~ s/\s+/\t/g; ## Replace spaces by tabulation
      next if ($line =~ /^;/) ; # skip comment lines
      #	&RSAT::message::Debug("line", $l, $line) if ($main::verbose >= 10);
      ## Create a new matrix if required
      if  ($line =~ /^\>(\S+)/) {

	## Parse ID and name from JASPAR header line
	my $id = &clean_id($1);
	my $postmatch = $'; #'
	my $name = $id;
	if ($postmatch =~ /\S+/) {
	  $name = &RSAT::util::trim($postmatch);
	  $name =~ s/\s+/_/g;
	}

	## Reset the alphabet for the new matrix (in the Nov 2013 release, JASPAR matrix format apparently changed)
	@temp_alphabet = qw(a c g t);

	## Instantiate new matrix
	$matrix = &NewJasparMatrix($id, $name);
	next;

      } elsif ($line =~ /^\s*(\S+)\s+/) {

	## JASPAR matrices are supposed to start with a fasta-like
	## header, as in the distribution file
	## http://jaspar.genereg.net/html/DOWNLOAD/jaspar_CORE/redundant/all_species/matrix_only/matrix_only.txt
	##
	## However, the Web interface provides matrices without header
	## e.g. http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?ID=MA0001.1&rm=present&collection=CORE
	##
	## If no header is provided, 
	unless ($matrix) {
	  if ($line =~ /\[/) {
	    my $id = "m".$current_matrix_nb;
	    my $name = $id;
	    $matrix = &NewJasparMatrix($id, $name);
	    &RSAT::message::Warning("Matrix does not contain the expected JASPAR header (row starting with \">\").") if ($main::verbose >= 1);
	  } else {
	    &RSAT::message::Warning("Skipped line $l", "does not conform to JASPAR format", $line) if ($main::verbose >= 2);
	  }
	}

	$line = &main::trim($line);
	$line =~ s/\[//;
	$line =~ s/\]//;
	$line =~ s/\s+/\t/;
	my @fields = split /\t/, $line;

	## The residue associated to the row
	my $residue = "";
	my $first_value = lc($fields[0]);
	if ($first_value =~ /\d+/) {
	  ## In MSCAN format, the first field of each row contains
	  ## the first cell of the matrix, instead of the residue (cf
	  ## JASPAR format).
	  ## http://www.cisreg.ca/cgi-bin/mscan/MSCAN
	  $residue = shift (@temp_alphabet);
	} else {
	  ## In Jaspar format, the residue is explicitly written at
	  ## the beginning of each row.
	  ## See: http://jaspar.cgb.ki.se/
	  $residue = lc(shift @fields);
	}
	$matrix->addIndexedRow($residue, @fields);
	&RSAT::message::Debug($line, $first_value,  join(";", @fields)) if ($main::verbose >= 5);
      }
    }
    close $in if ($file);

    foreach my $matrix (@matrices) {
#      $matrix->setAlphabet_lc(@temp_alphabet);
      $matrix->set_alphabet_for_type();
    }

    &InitializeEquiPriors(@matrices);
    return (@matrices);
} ## END readFromJasparFile

=pod

Instantiate a new matrix with a JASPAR record.

=cut

sub NewJasparMatrix {
  my ($id, $name) = @_;
  $matrix = new RSAT::matrix();
  &RSAT::message::Debug("_readFromJasparFile", $id, $name, $matrix) if ($main::verbose >= 5);
  $matrix->set_parameter("program", "jaspar");
  $ncol = 0;
  #@temp_alphabet = qw(a c g t);

  ## JASPAR header line contains an ID and a name
  $matrix->set_attribute("id", $id);
  $matrix->set_attribute("name", $name);
  $matrix->set_attribute("accession", $id); ## For compatibility with TRANSFAC format
  $matrix->set_attribute("description", join("", $id, " ", $name, "; from JASPAR"));
#  $matrix->setAlphabet_lc(@temp_alphabet);
  $matrix->set_alphabet_for_type();
  push @matrices, $matrix;
  $current_matrix_nb++;
  &RSAT::message::Info("line", $l, "new matrix", $current_matrix_nb, $name) if ($main::verbose >= 5);

  return ($matrix);
}

=pod

=item _readFromHomerFile($file)

Read a matrix from a file in HOMER format (files with extension
.homer). HOMER is a software suite for Software for motif discovery
and next-gen sequencing analysis (http://homer.salk.edu/homer/motif/).

This method is called by the method C<readFromFile($file, "homer")>.

=cut
sub _readFromHomerFile {
    my ($file, %args) = @_;
    &RSAT::message::Info(join("\t", "Reading matrix from HOMER file\t",$file)) if ($main::verbose >= 3);


    ## open input stream
    my ($in, $dir) = &main::OpenInputFile($file);

    ## Special treatment for the alphabet: sometimes indicated in the first column, sometimes not.
#    my @temp_alphabet = qw(a c g t);

    ## Initialize the matrix LIST
    local @matrices = (); ## updated in subroutine
    local $matrix; ## updated in subroutine
    local $current_matrix_nb = 1; ## updated in subroutine
    local $description = "Homer matrix"; 
    my $l = 0;
    my $ncol = 0;
    while ($line = <$in>) {
      $l++;
      next unless ($line =~ /\S/); ## Skip empty lines
      chomp($line); ## Suppress newline
      $line =~ s/\r//; ## Suppress carriage return
      $line =~ s/\s+/\t/g; ## Replace spaces by tabulation
      next if ($line =~ /^;/) ; # skip comment lines
      #	&RSAT::message::Debug("line", $l, $line) if ($main::verbose >= 10);
      ## Create a new matrix if required
      if  ($line =~ /^\>(\S+)/) {

	## Parse ID and name from HOMER header line
	my $name = &clean_id($1); ## In Homer format, first word is matrix name
	my $matrix_nb = scalar(@matrices) + 1; ## Matrix number is used in case automatic ID assignement would be needed
	my $id = "homer_".$matrix_nb.".".$name; ## If no ID is found in the header line, use name as default id
	my $postmatch = $'; #'

	## Instantiate new matrix
	$matrix = new RSAT::matrix();
	$matrix->set_parameter("program", "Homer");
	$ncol = 0;
	$matrix->force_attribute("id", $id);
	$matrix->set_attribute("AC", $id);
	$matrix->set_attribute("accession", $id);
	$matrix->set_attribute("name", $id);

	## Parse matrix description
	if ($postmatch =~ /(\S+)/) {
#	  $matrix->set_attribute("info", $postmatch);
	  $matrix->force_attribute("description", $name);
#	  while ($postmatch =~ /^([^:]+):(\S+)/) {
#	    $postmatch = $'; #'
	    my @tuples = split(/\,/, $postmatch);
	    foreach my $tuple (@tuples) {
	      if ($tuple =~/(\S+):(\S+)/) {
		my $key = $1;
		my $value = $2;
		&RSAT::message::Debug("&RSAT::MatrixReader::_readFromHomerFile()", $key, $value) if ($main::verbose >= 5);

		## Parse specific fields
		if (($key eq "T") && ($value =~ /(\d+\.\d+)\((\d+\.\d+)\%\)/)) {
		  $matrix->set_parameter("matching_sequences", &RSAT::util::round($1));
		  $matrix->set_parameter("sequence_coverage", $2);
		} elsif (($key eq "Multiplicity") && ($value =~ /(\d+\.\d+)/)) {
		  $matrix->set_parameter("multiplicity", $1);
		  if (&RSAT::util::IsReal($matrix->get_attribute("matching_sequences"))) {
		    my $nb_sites = &RSAT::util::round($matrix->get_attribute("matching_sequences") * $matrix->get_attribute("multiplicity"));
		    $matrix->set_parameter("nb_sites", $nb_sites);
		  }
		} else {
		  $matrix->set_parameter($key, $value);
		}
	      }
#	    }
	  }
	}

	# my @alphabet = qw(a c g t);
	# $matrix->setAlphabet_lc(@alphabet);
	$matrix->set_alphabet_for_type();
#	$matrix->force_attribute("nrow", 4);
	push @matrices, $matrix;
	$current_matrix_nb++;
	&RSAT::message::Info("line", $l, "new matrix", $current_matrix_nb, $name) if ($main::verbose >= 5);
	next;

      } elsif ($line =~ /^\s*(\S+)\s+/) {
	my $nb_sites = 100;
	if (&RSAT::util::IsNatural($matrix->get_attribute("nb_sites"))) {
	  $nb_sites = $matrix->get_attribute("nb_sites");
	}
	$line = &main::trim($line);
	$line =~ s/\s+/\t/;
	my @fields = split /\t/, $line;
#	&RSAT::message::Info("line", $l, "adding column", join(";", @fields)) if ($main::verbose >= 10);
	my $first_value = $fields[0];
	next unless (&RSAT::util::IsReal($first_value));
	my @values = ();
	foreach my $field (@fields) {
	  if (&RSAT::util::IsReal($field)) {
	    my $occurrences = sprintf("%.2f", $field*$nb_sites);
	    push(@values, $occurrences); ## Multiply by an arbitrary number to convert frequencies to counts
	  }
	}
	$matrix->addColumn(@values);
	$ncol++;
	$matrix->force_attribute("ncol", $ncol);
      }
    }

    close $in if ($file);

    foreach my $matrix (@matrices) {
#      $matrix->setAlphabet_lc(@temp_alphabet);
      $matrix->set_alphabet_for_type();
    }

    &InitializeEquiPriors(@matrices);
    return (@matrices);
} ## End readFromHomerFile Homer



=pod

=item _readFromMEMEFile($file)

Read a matrix from a MEME file. This method is called by the
method C<readFromFile($file, "MEME")>.

=cut


sub _readFromMEMEFile {
  my ($file) = @_;
  &RSAT::message::Info("Reading matrix from MEME BLOCK file (obsolete)\t", $file) if ($main::verbose >= 3);

  ## open input stream
  my ($in, $dir) = &main::OpenInputFile($file);
  my @matrices = ();
  my $current_matrix_nb = 0;
  my $matrix;
  my $current_col = 0;
  my $in_proba_matrix = 0;
  my $in_blocks = 0;
  my $width_to_parse = 0;
  my %alphabet = ();
  my %residue_frequencies = ();
  my @alphabet = ();
  my @frequencies = ();
#  my $parsed_width = 0;
  my $l = 0;
  my $meme_command = "";
  my $prefix = "matrix";
  if ($file) {
      $prefix = &RSAT::util::ShortFileName($file);
  }
  while (<$in>) {
    $l++;
    next unless (/\S/);
    s/\r//;
    chomp();
    $_ = &main::trim($_);
    if (/MOTIF\s+(\d+)\s+width =\s+(\d+)\s+sites =\s+(\d+)\s+llr =\s+(\d+)\s+E-value =\s+(\S+)/) {
      &RSAT::message::Debug("line", $l, "Parsing matrix parameters") if ($main::verbose >= 5);

      $current_matrix_nb = $1;
      $width_to_parse = $2;
      $matrix = new RSAT::matrix();
      $matrix->init();
      $matrix->set_parameter("program", "meme");
      $matrix->set_parameter("matrix.nb", $current_matrix_nb);

#      my $id = "meme_motif_".$current_matrix_nb;

      $matrix->set_parameter("id", $prefix.$current_matrix_nb);

      $matrix->set_parameter("id", $id);
      $matrix->set_parameter("ac", $id); ## For TRANSFAC compatibility
      $matrix->set_parameter("name", $id); ## For readability of logos
      $matrix->set_attribute("ncol", $2);

      $matrix->set_parameter("command", $meme_command);
      $matrix->set_parameter("sites", $3);
      $matrix->set_parameter("meme.llr", $4);
      $matrix->set_parameter("meme.E-value", $5);
      $matrix->setPrior(%residue_frequencies);
#      &RSAT::message::Debug("line", $l, "Read letter frequencies", %residue_frequencies) if ($main::verbose >= 10);
      # $matrix->setAlphabet_lc(@alphabet);
      $matrix->setAlphabet(@alphabet); # For MEME, we take the alphabet from the file itself
      # $matrix->set_alphabet_for_type()
#      $matrix->force_attribute("nrow", scalar(@alphabet)); ## Specify the number of rows of the matrix
      push @matrices, $matrix;

      ## Meme command
    } elsif (/^command: /) {
      $meme_command = $'; #'

    } elsif (/Background letter frequencies/) {
      my $alphabet = <$in>;
#      $alphabet = lc($alphabet);
      $alphabet = &main::trim($alphabet);
      %residue_frequencies = split /\s+/, $alphabet;
      @alphabet = sort (keys %residue_frequencies);
#      &RSAT::message::Debug("line", $l, "Read letter frequencies", %residue_frequencies) if ($main::verbose >= 10);

      ## Index the alphabet
      foreach my $l (0..$#alphabet) {
	$alphabet{$alphabet[$l]} = $l;
      }
      $matrix->setAlphabet(@alphabet);

      ## Parse BLOCKS format
    } elsif (/Motif (\d+) in BLOCKS format/) {
      $current_matrix_nb = $1;
      $in_blocks = 1;
      &RSAT::message::Debug("line", $l, "Starting to parse BLOCKS format") if ($main::verbose >= 10);

    } elsif ($in_blocks) {
      if (/(\S+)\s+\(\s*\d+\)\s+(\S+)/) {
	my $seq_id = $1;
	my $seq = lc($2);
	my $seq_len =  length($seq);
	&RSAT::message::Debug("Reading sequences for MEME matrix", "SeqId", $seq_id, "\n\t", $seq) if ($main::verbose >= 10);
	if ($seq_len > 0) {
#	  $parsed_width = &main::max($parsed_width, $seq_len);
	  $matrix->add_site(lc($seq), score=>1, id=>$seq_id, max_score=>0);
	}

      } elsif (/\/\//) {
	&RSAT::message::Debug("line", $l, "BLOCKS format parsed") if ($main::verbose >= 5);
	$in_blocks = 0;

      }
    }
  }
  close $in if ($file);

  return @matrices;
#  return $matrices[0];
}

=pod

=item _readFromMEMEFile_2015($file)

Read a matrix from a MEME file, as described in the format on the MEME page. This method is called by the
method C<readFromFile($file, "MEME")>.

MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.303 C 0.183 G 0.209 T 0.306 

MOTIF crp alternative name
letter-probability matrix: alength= 4 w= 19 nsites= 17 E= 4.1e-009 
 0.000000  0.176471  0.000000  0.823529 
 0.000000  0.058824  0.647059  0.294118 
 0.000000  0.058824  0.000000  0.941176 
 0.176471  0.000000  0.764706  0.058824 
 0.823529  0.058824  0.000000  0.117647 
 0.294118  0.176471  0.176471  0.352941 
 0.294118  0.352941  0.235294  0.117647 
 0.117647  0.235294  0.352941  0.294118 
 0.529412  0.000000  0.176471  0.294118 
 0.058824  0.235294  0.588235  0.117647 
 0.176471  0.235294  0.294118  0.294118 
 0.000000  0.058824  0.117647  0.823529 
 0.058824  0.882353  0.000000  0.058824 
 0.764706  0.000000  0.176471  0.058824 
 0.058824  0.882353  0.000000  0.058824 
 0.823529  0.058824  0.058824  0.058824 
 0.176471  0.411765  0.058824  0.352941 
 0.411765  0.000000  0.000000  0.588235 
 0.352941  0.058824  0.000000  0.588235 

=cut


sub _readFromMEMEFile_2015 {
  my ($file, %args) = @_;
  &RSAT::message::Info("Reading matrix from meme file version 2015\t", $file) if ($main::verbose >= 3);

  ## open input stream
  my ($in, $dir) = &main::OpenInputFile($file);
  my @matrices = ();
  my $current_matrix_nb = 0;
  my $matrix;
  my $nb_sites = 0;
  my $ncol = 0;
  my $in_blocks = 0;
  my %alphabet = ();
  my @alphabet = ();
  my %residue_frequencies = ();
  my @frequencies = ();
  my $l = 0;
  my $prefix = "matrix";
  if ($file) {
      $prefix = &RSAT::util::ShortFileName($file);
  }
  while (<$in>) {
    $l++;
    s/\r//;
    chomp();
    $_ = &main::trim($_);

    if (/Background letter frequencies/) {
      ## Read alphabet from the background residue frequencies
      my $bg_freq = <$in>;
#      $bg_freq = lc($bg_freq);
      $bg_freq = &main::trim($bg_freq);
      my @residue_frequencies = split /\s+/, $bg_freq;
      while (@residue_frequencies) {
	my $residue = shift(@residue_frequencies);
	my $frequency = shift(@residue_frequencies);
	push @alphabet, $residue;
	$residue_frequencies{$residue} = $frequency;
	push @frequencies, $frequency;
	#&RSAT::message::Debug("MEME background frequencies", $residue, $frequency) if ($main::verbose >= 10);
      }

#      die ("HELLO");
#      %residue_frequencies = split /\s+/, $bg_freq;
#      @alphabet = sort (keys %residue_frequencies);     
#      &RSAT::message::Debug("line", $l, "Read letter frequencies", %residue_frequencies) if ($main::verbose >= 10);

      ## Index the alphabet
      foreach my $l (0..$#alphabet) {
	$alphabet{$alphabet[$l]} = $l;
      }

    # MOTIF crp alternative name
    } elsif (/MOTIF\s+(\S+)\s*(\S*)/) {
      $current_matrix_nb += 1;
      $ncol = 0;
      $matrix = new RSAT::matrix();
      $matrix->init();
      $matrix->set_parameter("program", "meme");
      $matrix->set_parameter("matrix.nb", $current_matrix_nb);

      my $id = $1;
      my $ac = $2 || $id;
      
      $matrix->set_parameter("id", $id);
      $matrix->set_parameter("ac", $ac); ## For TRANSFAC compatibility
      $matrix->set_parameter("name", $id); ## For readability of logos
      if (defined($args{"type"})) {
	$matrix->set_parameter("type", $args{type});
#	&RSAT::message::Debug("matrix type", $matrix->get_attribute("type")); die("HELLO");
      }
      $matrix->setAlphabet(@alphabet);
      $matrix->setPrior(%residue_frequencies);
      
      &RSAT::message::Debug("motif found line", $l, "id=$id ac=$ac name =$id") if ($main::verbose >= 5);
    
    } elsif (/letter-probability matrix:.*w=\s*(\d+)\s+nsites=\s*(\d+)\s+E=\s*(\S+)/) {
      ## letter-probability matrix: alength= alphabet length w= motif length nsites= source sites E= source E-value
      $nb_sites = $2;
      $matrix->set_parameter("sites", $nb_sites);
      $matrix->set_parameter("meme.E-value", $3);
#      $matrix->setPrior(%residue_frequencies);
      #$matrix->setAlphabet_lc(@alphabet);
 #     $matrix->setAlphabet(@alphabet);
      push @matrices, $matrix;
#       &RSAT::message::Debug("MEME matrix read", 
# 			    $matrix->get_attribute("id"),
# 			    $matrix->get_attribute("sites"),
# #			    join(" ", $matrix->get_attribute("prior"))
# 	  ) if ($main::verbose >= #0);
      
      &RSAT::message::Debug("line", $l, "ncol=$1,nbsites=$nb_sites") if ($main::verbose >= 5);

      ## next lines has the probabilities infos
      $in_blocks = 1;
      &RSAT::message::Debug("line", $l, "Starting to parse the matrix cells") if ($main::verbose >= 10);

    } elsif ($in_blocks) {
      
      if ($_ !~ /\d+/) {
	&RSAT::message::Debug("line", $l, "matrix cells parsed") if ($main::verbose >= 5);
	$in_blocks = 0;

# 	## Converting the frequency matrix to a count matrix
# 	my @matrix = $matrix->getMatrix();
# 	my $ncol = $matrix->ncol();
# 	my $nrow = $matrix->nrow();
# 	my $nb_sites = $matrix->get_attribute("sites");
# 	&RSAT::message::Debug("Parsing matrix", $matrix->get_attribute("id"), 
# 			      "ncol=".$ncol,
# 			      "nrow=".$nrow,
# 			      "sites=".$nb_sites) if($main::verbose >= 10);
# 	for my $c (0..($ncol-1)) {
# 	  for my $r (0..($nrow-1)) {
# #	    $matrix[$c][$r] *= 100; #multiply by 100 to obtain count over 100
# #	    $matrix[$c][$r] = $matrix[$c][$r] * $nb_sites / 100 ; ## put back to the number of contributing sites
# 	    $matrix[$c][$r] = $matrix[$c][$r] * $nb_sites; ## put back to the number of contributing sites	    
# 	    $matrix[$c][$r] = int($matrix[$c][$r] + 0.5); ## round to integer
# #	    $matrix[$c][$r]  = 0;
# 	  }
# 	}
# 	$matrix->setMatrix(@matrix);
			
      } else {
    	## Count column of the matrix file (row like in TRANSFAC format)
	my $values = $_;
	$values = &RSAT::util::trim($values);
	my @fields = split /\s+/, $values;
	&RSAT::message::Debug("line ".$l, "adding column", join (":", @fields)) if ($main::verbose >= 5);
	$matrix->addColumn(@fields);
	$ncol++;
	$matrix->force_attribute("ncol", $ncol);
	next;
      }
    }
    }
  close $in if ($file);

  ## Multiply frequencies by the number of sites to obtain counts
  foreach my $matrix (@matrices) {
    ## Converting the frequency matrix to a count matrix
    my @matrix = $matrix->getMatrix();
    my $ncol = $matrix->ncol();
    my $nrow = $matrix->nrow();
    my $nb_sites = $matrix->get_attribute("sites");
    # &RSAT::message::Debug("Parsing matrix", $matrix->get_attribute("id"), 
    # 			  "ncol=".$ncol,
    # 			  "nrow=".$nrow,
    # 			  "sites=".$nb_sites) if($main::verbose >= 10);
    for my $c (0..($ncol-1)) {
      for my $r (0..($nrow-1)) {
#	    $matrix[$c][$r] *= 100; #multiply by 100 to obtain count over 100
#	    $matrix[$c][$r] = $matrix[$c][$r] * $nb_sites / 100 ; ## put back to the number of contributing sites
	$matrix[$c][$r] = $matrix[$c][$r] * $nb_sites; ## put back to the number of contributing sites	    
	$matrix[$c][$r] = int($matrix[$c][$r] + 0.5); ## round to integer
#	    $matrix[$c][$r]  = 0;
      }
    }
    $matrix->setMatrix(@matrix);
  }

  return @matrices;

}


=pod

=item _readFromSeq($file)

Read a matrix from a sequence file containing the pre-aligned sites.
The method just reads the sequences and counts the residue frequencies
at each position.

=cut

sub _readFromSeq {
  my ($file, %args) = @_;
  &RSAT::message::Info("Reading matrix from feature file\t", $file) if ($main::verbose >= 3);
  my $seq_format = $args{seq_format} || "fasta";

  ## open input stream
  my ($in, $dir) = &main::OpenInputFile($file);
  my @alphabet =  ("A", "C", "G", "T");

  my $matrix = new RSAT::matrix();
  $matrix->init();
  $matrix->set_parameter("matrix.nb", 1);
  $matrix->set_attribute("name", $matrix_name);
#      $matrix->set_attribute("id", $matrix_id);
#  $matrix->setAlphabet_lc(@alphabet);
  $matrix->set_alphabet_for_type();
#      $matrix->force_attribute("nrow", scalar(@alphabet)); ## Specify the number of rows of the matrix
#  $matrices{$matrix_name} = $matrix;

  my $site_nb = 0;
  while ((($site_seq, $site_id, @comments) = &main::ReadNextSequence($in, $seq_format, $dir)) &&
	 (($site_seq) || ($site_id))) {
    $site_nb++;
#    &RSAT::message::Debug("Adding sequence to matrix", $site_nb, $site_seq, $site_id) if ($main::verbose >= 10);
    if ($site_nb == 1) {
      $matrix->set_attribute("ncol", length($site_seq));
    }
    $matrix->add_site(lc($site_seq), id=>$site_id, score=>1);
  }
  close $in if ($file);
  return ($matrix);
}


=pod

=item _readFromFeatureFile($file)

Read a matrix from a feature file (he input of feature-map). 

This method is called by the method C<readFromFile($file, "feature")>.

The main usage is to retrieve a collection of sites resulting from
matrix-scan, in order to build a new collection of matrices from these
sites. 

The third column of the feature file (containing the feature name) is
used as matrix name. If several feature names are present in the
feature file, several matrices are returned accordingly. The 7th
column, which contains the sequence of the feature, is used to build
the matrix (or matrices).

=cut
sub _readFromFeatureFile {
  my ($file) = @_;
  &RSAT::message::Info("Reading matrix from feature file\t", $file) if ($main::verbose >= 3); 

  ## open input stream
  my ($in, $dir) = &main::OpenInputFile($file);
  my @matrices = (); 
  my %matrices = (); ## Matrices are indexed by name
  my @alphabet =  ("A", "C", "G", "T");
  my $current_matrix_nb = 0;
  my $matrix;
  my $l = 0;

  while (my $line = <$in>) {
    $l++;
    next if ($line =~ /^;/);
    next if ($line =~ /^--/);
    next if ($line =~ /^#/);
    next unless ($line =~ /\S/);
    $line =~ s/\r//;
    chomp($line);
    ## Instantiate an object to store feature information
    my $feature = new RSAT::feature();
    $feature->parse_from_row($line, "ft");
    my $matrix_name = $feature->get_attribute("feature_name") || "matrix_".$current_matrix_nb;
    my $matrix_id = $feature->get_attribute("feature_name") || "matrix_".$current_matrix_nb;
    my $site_sequence = $feature->get_attribute("description");
    my $site_id = join ("_", 
			$feature->get_attribute("seq_name"),
			$feature->get_attribute("feature_name"),
			$feature->get_attribute("strand"),
			$feature->get_attribute("start"),
			$feature->get_attribute("end"),
		       );
#    &RSAT::message::Debug("&RSAT::MatrixReader", $matrix_name,"feature parsed", $l, $site_sequence, $site_id) if ($main::verbose >= 5);
    if (defined($matrices{$matrix_name})) {
      $matrix = $matrices{$matrix_name};
    } else {
      $current_matrix_nb++;
      $matrix = new RSAT::matrix();
      $matrix->init();
      $matrix->set_parameter("program", "feature");
      $matrix->set_parameter("matrix.nb", $current_matrix_nb);
      $matrix->set_attribute("name", $matrix_name);
      $matrix->set_attribute("id", $matrix_id);
      $matrix->set_attribute("ncol", length($site_sequence));
#      $matrix->set_parameter("sites", $3);
#      $matrix->setPrior(%residue_frequencies);
#      &RSAT::message::Debug("line", $l, "Read letter frequencies", %residue_frequencies) if ($main::verbose >= 10);
      #$matrix->setAlphabet_lc(@alphabet);
      #$matrix->setAlphabet(@alphabet);
      $matrix->set_alphabet_for_type();
#      $matrix->force_attribute("nrow", scalar(@alphabet)); ## Specify the number of rows of the matrix
      $matrices{$matrix_name} = $matrix;
      push @matrices, $matrix;
    }
    $matrix->add_site(lc($site_sequence), id=>$site_id,max_score=>0,
		      "score"=>1, ## Here we don't want to add up the scores, because we want to count the residue occurrences
		    );
  }
  close $in if ($file);
  return @matrices;
}


=pod

=item _readFromMotifSamplerFile($file)

Read a matrix from a I<MotifSampler> B<site> file (i.e. the file
specified with the option -o of MotifSampler).

MotifSampler is part of the software suite INCLUSive
(http://homes.esat.kuleuven.be/~thijs/download.html), developed by
Gert Thijs.

MotifSampler export two files: a description of the sites (option -o)
and of the matrices (option -m). We prefer to parse the file with
sites, because it is more informatiive for the following reasons: (1)
it contains the site sequences; (2) the matrix descriptin is in
relative frequencies rather than residue occurrences, which biases the
pseudo; (3) the site file contains additional statistics (IC, CS).

This method is called by the method C<readFromFile($file,
"MotifSampler")>.

=cut
sub _readFromMotifSamplerFile {
  my ($file, %args) = @_;
  &RSAT::message::Info(join("\t", "Reading matrix from MotifSampler output file\t",$file)) 
    if ($main::verbose >= 3);

  ## open input stream
  my ($in, $dir) = &main::OpenInputFile($file);

  ## Initialize the matrix list
  my @matrices = ();
  my @alphabet = qw(a c g t);
  my %prior;
  foreach my $letter (@alphabet) {
    $prior{lc($letter)} = 1/scalar(@alphabet);
  }
  my $ncol = 0;
  my $matrix;			## the amtrix object
  while (<$in>) {
    next unless /\S/;		## Skip empty lines
    next if /^#*$/;		## Skip empty lines
    if (/^#id:\s*(.*)/i) {

      ## The ID row also contains the matrix parameters
      my $motif_desc = $1;
      chomp($motif_desc);
      my @fields = split(/\s+/, $motif_desc);
      my $id = shift (@fields); $id = &clean_id($id);
      
      $matrix = new RSAT::matrix();
      $matrix->set_parameter("program", "MotifSampler");
      $matrix->set_parameter("id", $id);
      while (my $field = shift @fields) {
	$field =~ s/:$//;
	$value = shift @fields;
	if ($field eq "instances") {
	  $matrix->set_parameter("sites", $value);
	} else {
	  $matrix->set_parameter("MS.".$field, $value);
	}
      }
      #$matrix->setAlphabet_lc(@alphabet);
      $matrix->set_alphabet_for_type();
      $matrix->setPrior(%prior);
      $matrix->force_attribute("nrow", 4);

      push @matrices, $matrix;
    } elsif (/^#/) {
      ## Skip comment line
      next;
    } else {
      my ($seq_id, $program, $site_type, $start, $end, $score, $strand, $frame, @the_rest) = split (/\s+/, $_);
      my $desc = join(" ", @the_rest);
      $strand =~ s/\+/D/;
      $strand =~ s/\-/R/;
      my $site_seq;
      my $site_id = join ("_", $seq_id, $start, $end, $strand, $score);
      if ($desc =~ /site\s+"(\S+)"/i) {
#	die;
	$site_seq = lc($1);
	$matrix->add_site($site_seq, id=>$site_id, score=>1, max_score=>0);
#	$ncol = &RSAT::stats::max($ncol, length($site_seq));
#	$matrix->force_attribute("ncol", $ncol);
#	$matrix->treat_null_values();
      }
    }
  }
  return (@matrices);
}


=pod

=item _readFromMotifSamplerMatrixFile($file)

Read a matrix from a I<MotifSampler> B<matrix> file (i.e. the file
specified with the option -o of MotifSampler).

MotifSampler is part of the software suite INCLUSive
(http://homes.esat.kuleuven.be/~thijs/download.html), developed by
Gert Thijs.

MotifSampler export two files: a description of the sites (option -o)
and of the matrices (option -m). We prefer to parse the file with
sites, because it is more informatiive for the following reasons: (1)
it contains the site sequences; (2) the matrix descriptin is in
relative frequencies rather than residue occurrences, which biases the
pseudo; (3) the site file contains additional statistics (IC, CS).

This method is called by the method C<readFromFile($file,
"MotifSampler")>.

=cut

sub _readFromMotifSamplerMatrixFile {
    my ($file, %args) = @_;
    &RSAT::message::Info(join("\t", "Reading matrix from MotifSampler matrix file\t",$file)) 
      if ($main::verbose >= 3);

    ## open input stream
    my ($in, $dir) = &main::OpenInputFile($file);

    ## Initialize the matrix list
    my @matrices = ();
    my @alphabet = qw(a c g t);
    my $ncol = 0;
    my $matrix; ## the amtrix object
    while (<$in>) {
      next unless /\S/; ## Skip empty lines
      next if /^#*$/; ## Skip empty lines

      if(/^#ID\s*=\s*(\S+)/) {
	my $id = &clean_id($1);
	$matrix = new RSAT::matrix();
	$matrix->set_parameter("program", "MotifSampler");
	$matrix->set_attribute("AC", $id);
	$matrix->set_attribute("id", $id);
	$matrix->force_attribute("nrow", 4);
	# $matrix->setAlphabet_lc("a","c","g","t");
	$matrix->set_alphabet_for_type();
	push @matrices, $matrix;
      } elsif (/^#Score = (\S+)/i) {
	$matrix->set_parameter("MS.score", $1);
      } elsif (/^#Consensus = (\S+)/i) {
	$matrix->set_parameter("MS.consensus", $1);
	for my $i (1..$ncol) {
	  my $line = (<$in>);
	  my @values = split (/\s+/, $line);
	  $matrix->addColumn(@values);
	  $matrix->force_attribute('ncol', $i);
	}
#	$matrix->force_attribute("ncol", $ncol);
      } elsif (/^#W = (\S+)/i) {
	$ncol = $1;
      }
    }
    return (@matrices);
}


=pod

=item _readFromClustalFile($file)

Read a matrix from a multiple alignment in clustal format (extension
 .aln).  This method is called by the method C<readFromFile($file,
 "clustal")>.

=cut
sub _readFromClustalFile {
    my ($file) = @_;

    my @matrices = ();
    my $matrix = new RSAT::matrix();
    $matrix->set_parameter("program", "clustal");
    push @matrices, $matrix;

    ## open input stream
#    my $in = STDIN;
    my ($in, $dir) = &main::OpenInputFile($file);
#    if ($file) {
#	open INPUT, $file;
#	$in = INPUT;
#    }

    ## Check the header
    my $header = <$in>;
    unless ($header =~ /clustal/i) {
	&main::Warning("This file does not contain the clustal header");
    }

    ## Read the sequences
    my %sequences = ();
    &RSAT::message::Debug("&RSAT::MatrixRedaer::_readFromClustalFile()", "Reading sequences") if ($main::verbose >= 5);
    while (<$in>) {
	next unless (/\S/);
	s/\r//;
	chomp();
	if (/^\s*(\S+)\s+(.+)$/) {
	    my $seq_id = $1;
	    next if ($seq_id eq "*"); ## asterisks are used to mark conservation
	    my $new_seq = $2;
	    
	    ## index the new sequence
	    $sequences{$seq_id} .= $new_seq;
	    &RSAT::message::Debug("\tSequence", $seq_id, 
				  length($new_seq), length($sequences{$seq_id}))
		if ($main::verbose >= 5);
	}
    }
    
    ## Calculate count matrix
    my %matrix = ();
    my %prior = ();
    my $ncol = 0;
    my $nrow = 0;
    &RSAT::message::Info("Calculating profile matrix from sequences") if ($main::verbose >= 5);
    foreach my $seq_id (sort keys %sequences) {
	my $sequence = $sequences{$seq_id};
	$sequence =~ s/\s+//g;

	################################################################
	## Distinguish between insertions and leading/trailing gaps
	$terminal_gap_char = ".";

	## Substitute leading gaps
	if ($sequence =~ /^(\-+)/) {
	    $leading_gap_len = length($1);
	    my $leading_gap = ${terminal_gap_char}x$leading_gap_len;
	    $sequence =~ s|^(\-+)|${leading_gap}|;
	}
	## Substitute trailing gaps
	if ($sequence =~ /(\-+)$/) {
	    $trailing_gap_len = length($1);
	    my $trailing_gap = ${terminal_gap_char}x$trailing_gap_len;
	    $sequence =~ s|(\-+)$|${trailing_gap}|;
	}
	&RSAT::message::Debug($seq_id,$sequence) if ($main::verbose >= 5);
	    
	$ncol = &main::max($ncol, length($sequence));
	&RSAT::message::Debug("Sequence", $seq_id, length($sequence)) if ($main::verbose >= 5);
	my @sequence = split '|', $sequence;
	foreach my $i (0..$#sequence) {
	    my $res = lc($sequence[$i]);
	    next if ($res eq "N"); ## BEWARE: THIS IS FOR DNA ONLY
#	    next if ($res eq "-");
	    next if ($res eq "."); ## leading and trailing gaps
	    next if ($res eq "*");
	    $prior{$res}++;
	    $matrix{$res}->[$i] += 1;
	}
    }
    $matrix->set_attribute("ncol", $ncol);

    ## Define prior probabilities, alphabet, and matrix size
    my @alphabet = sort keys %prior;
    my $alpha_sum = 0;
    foreach my $res (@alphabet) {
	$alpha_sum += $prior{$res};
    }
    foreach my $res (@alphabet) {
	if ($alpha_sum > 0) {
	    $prior{$res} /= $alpha_sum;
	} else {
	    $prior{$res} = 0;
	}
    }
    $matrix->setPrior(%prior);

    ## Store the matrix
    my @matrix = ();
    foreach my $r (0..$#alphabet) {
	my $res = $alphabet[$r];
	my @row = @{$matrix{$res}};
	$nrow++;
	foreach $i (0..($ncol-1)) {
	    $row[$i] = 0 unless (defined($row[$i]));
	}
	$matrix->addRow(@row);
	&RSAT::message::Debug("Adding row", $r, $res, join ":", @row), if ($main::verbose >= 5);
    }
    #$matrix->setAlphabet_lc(@alphabet);
    $matrix->set_alphabet_for_type();
    $matrix->force_attribute("ncol", $ncol);
#    $matrix->force_attribute("nrow", $nrow);

    &RSAT::message::Debug ("Matrix size",
			   $nrow,
			   $ncol,
			   $matrix->nrow(),
			   $matrix->ncol())
			     if ($main::verbose >= 5);
    close $in if ($file);

    return (@matrices);
}


=pod

=item B<SortMatrices>

Sort matrices according to the value of some parameter.

Usage:
  my @sorted_matrices = &RSAT::MatrixReader::SortMatrices ($sort_key, $sort_order, @matrices);

Parameters
  $sort_key   a parameter whose value will determine the sorting
  $sort_order desc (descending), asc (ascending) or alpha (alphabetical)
  @matrices    a list of matrices (objects belonging to the class RSAT::matrix)

=cut
sub SortMatrices {
  my ($sort_key, $sort_order, @matrices) = @_;
  my $nb_matrices = scalar(@matrices);

  ## Check if there is at least one matrix
  if ($nb_matrices > 0) {

    ## Check that all matrices have the sort key as attribute
    my %key_value = ();
    my $attr_not_found = 0;
    my $not_real = 0;
    foreach my $matrix (@matrices) {
      if (my $value = $matrix->get_attribute($sort_key)) {
	$key_value{$matrix} = $value;
	unless (&RSAT::util::IsReal($value)) {
	  &RSAT::message::Debug("&RSAT::MatirxReader::SortMatrices()", "matrix", $matrix->get_attribute("id"), 
				"Non-real attribute", $sort_key, "value", $value) if ($main::verbose >= 2);
	  $not_real++;
	}
      } else {
	$attr_not_found++;
      }
    }

    if ($attr_not_found > 0) {
      &RSAT::message::Warning("Cannot sort matrices by", $sort_key,
			      "because this attribute is missing in", $attr_not_found."/".$nb_matrices, "matrices")
	if ($main::verbose >= 1);
    } elsif ($not_real > 0) {
      &RSAT::message::Warning("Cannot sort matrices by", $sort_key,
			      "because this attribute has non real values in", $not_real."/".$nb_matrices, "matrices")
	if ($main::verbose >= 1);
    } else {
      #    ## Check if the first matrix has the sort key as attribute
      #    my $first_matrix = $matrices[0];
      #    if ($first_matrix->get_attribute($sort_key)) {
      &RSAT::message::Info("Sorting", $nb_matrices, "matrices by", $sort_order, $sort_key) if ($main::verbose >= 2);

      ## Sort matrices
      if ($sort_order eq "desc") {
	&RSAT::message::Warning("Sorting matrices by descending values of", $sort_key) if ($main::verbose >= 2);
	@matrices = sort {$b->{$sort_key} <=> $a->{$sort_key}} @matrices;
      } elsif ($sort_order eq "asc") {
	&RSAT::message::Warning("Sorting matrices by ascending values of", $sort_key) if ($main::verbose >= 2);
	@matrices = sort {$a->{$sort_key} <=> $b->{$sort_key}} @matrices;
      } elsif ($sort_order eq "alpha") {
	&RSAT::message::Warning("Sorting matrices by alphabetic values of", $sort_key) if ($main::verbose >= 2);
	@matrices = sort {lc($a->{$sort_key}) cmp lc($b->{$sort_key})} @matrices;
      } else {
	&RSAT::error::FatalError($sort_order, "is not a valid sorting order. Supported: desc,asc,alpha.");
      }
    }
  }

  ## Check sorting (debugging)
  if ($main::verbose >= 5) {
    &RSAT::message::Info("Sorted", $nb_matrices, "matrices by", $sort_order, $sort_key);
    my $m = 0;
    foreach my $matrix (@matrices) {
      $m++;
      &RSAT::message::Debug("sorted matrix", $m, $matrix->get_attribute("id"), $sort_key, $matrix->get_attribute($sort_key));
    }
  }

  return @matrices;

}

=pod

=item B<SetMatrixName(matrix, matrix_number_in_input_file, matrix_file, matrix_format)>

=cut

################################################################
## Choose the name of a matrix
sub SetMatrixName {
  my ($matrix, $m, $matrix_file, $input_format) = @_;

  ## Make sure that the matrix has the mandatory attribute "id"
  my $matrix_id = $matrix->get_attribute("id");
  unless ($matrix_id) {
    if ($matrix_file) {
      ($matrix_id) =  &RSAT::util::ShortFileName($matrix_file);
      $matrix_id =~ s/\.${input_format}$//; ## suppress the extension from the file name if it corresponds to the matrix format
      $matrix_id =~ s/\.txt$//; ## suppress .txt extension
      $matrix_id .= "_m".$m;
      &RSAT::message::Debug("Matrix", $m."/".scalar(@matrices), "name", $matrix_id) if ($main::verbose >= 5);
    } else {
      $matrix_id = "matrix";
      $matrix_id .= "_m".$m;
    }
    $matrix->force_attribute("id", $matrix_id);
  }

  ## Get accession number (preferred name for TRANSFAC input format)
  my $matrix_ac = $matrix->get_attribute("accession") || $matrix->get_attribute("AC");
  unless ($matrix_ac) {
    $matrix_ac = $matrix_id;
    $matrix->force_attribute("AC", $matrix_ac);
    $matrix->force_attribute("accession", $matrix_ac);
  }

  ## Check that the matrix has a name
  my $matrix_name = $matrix->get_attribute("name");
  unless ($matrix_name){
    $matrix_name = $matrix_ac;
    $matrix->force_attribute("name", $matrix_name);
  }
  &RSAT::message::Info($m,
		       "name=".$matrix->get_attribute("name"),
		       "id=".$matrix->get_attribute("id"),
		       "format=".$input_format,
		       "file=".$matrix_file,
		      ) if ($main::verbose >= 5);
}


return 1;

__END__

