
#########################################################################################
############################## PATHWAY DIAGRAM MANAGEMENT ###############################
#########################################################################################

package GV::GraphUtil;
{
  sub calc_string_width {
    my $small_font_width = 6;
    my ($string) = @_;
    my $string_length = length($string); 
    my $string_width =  $string_length * $small_font_width;
    return $string_width;
  }

  sub calc_string_height {
    my $small_font_height = 8;
    ### trivial for the time being, but could be more elaborate for multi-line labels
    return $small_font_height;
  }
}


package classes::PathwayDiagram;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "pwdg_";
    @_objects = ();
    %_name_index = ();
    %_id_index = ();
    %_attribute_count = ();
    %_attribute_cardinality = (id=>"SCALAR",
			       title=>"SCALAR",
			       element_count=>"SCALAR",
			       xsize=>"SCALAR",
			       ysize=>"SCALAR",
			       pathway=>"SCALAR",
			       steps=>"HASH",
			       nodes=>"ARRAY", ### ARRAY of PathwayDiagramNode
			       id_index=>"HASH",
			       node_index=>"HASH",
			       arcs=>"ARRAY", ### ARRAY of PathwayDiagramArc
			       source=>"SCALAR"); 
    
    sub get_object {
	### return an object ref given its ID or one of its names
	### case insensitive : all querys are converted to uppsercase
	my ($diagram, $query) = @_;
	my $uc_query = uc($query);
	my $id = undef;
	my $obj = undef;
	
	### find query in id index
	if (defined(${$diagram->{_id_index}}{$uc_query})) {
	    return ${$diagram->{_id_index}}{$uc_query};
        } 

        ### fail to identify object
        return undef; 
    }

    #### return the number of nodes in the diagam
    sub size {
	my ($diagram) = @_;
	my @nodes = $diagram->get_attribute("nodes");
	return $#nodes+1;
    }

    sub get_node { 

	### return a node object given its ID or one of its names case
	### insensitive : the index is in uppercase and all queries
	### are converted to uppercases
	my ($diagram, $query) = @_;
	my $uc_query = uc($query);
	my $node_id;
	my $obj;
	my $r_ids;


	### find query in node index
	if ($node_id = ${$diagram->{_node_index}}{$uc_query}) {
            warn "node ID found\t", $node_id, "\n" if ($main::verbose >= 4); 
            if ($obj = ${$diagram->{_id_index}}{$node_id}) {
                warn "node object found\t", $obj, "\n" if ($main::verbose >= 4); 
	    } else {
		warn "WARNING: node object not found\n" if ($main::verbose >= 4);
	    } 
        } else {
            warn "node not found\n" if ($main::verbose >= 4); 
	}
	warn join ("\t", "getting node", $query, $uc_query, $node_id, $obj), "\n" if ($main::verbose >= 4);
        return ($id, $obj);
     }

     sub incr_element_count { 
          my ($diagram) = @_;
          my $next_count = $diagram->get_attribute("element_count") + 1;
          $diagram->force_attribute("element_count", $next_count);
          return $next_count;
     }


    sub add_node {
	#### usage :
	####    $diagram->add_node(id=>$id, label=>$label);
	my ($diagram, %args) = @_;
	my $element_nb = $diagram->incr_element_count();
	my $node = new classes::PathwayDiagramNode (id=>$element_nb, ### default id is element nb
						 %args);
	my $id = $node->get_attribute("id");
	my $label = $node->get_attribute("label");
	$diagram->push_attribute("nodes", $node);
	${$diagram->{_id_index}}{$id} = $node;
        ${$diagram->{_node_index}}{uc($id)} = $id;
        ${$diagram->{_node_index}}{uc($label)} = $id;
        return $node;
    }
  
  sub add_arc {
    my ($diagram, %args) = @_;
    warn "; new arc from ", $args{from}, "\tto ", $args{to}, "\n" if ($main::verbose >= 3);
    unless (defined($args{from}) &&
	    defined($args{to})) {
	die "Error: arc cannot be created withour specifying attributes 'from' and 'to'\n";
    }
    my $element_nb = $diagram->incr_element_count();
    my $arc = new classes::PathwayDiagramArc (id=>$element_nb, ### default id is element nb
					   %args);
    my $id = $arc->get_attribute("id");
    ${$diagram::_id_index{uc($id)}} = $arc;
    $diagram->push_attribute("arcs", $arc);
    return $arc;
  }
  
  sub calc_size {
    my ($diagram) = @_;
    my $xmax = 0;
    my $ymax = 0;
    foreach my $node ($diagram->get_attribute("nodes")) {
      my $xpos = $node->get_attribute("xpos");
      $max = $xpos if ($xpos > $xmax);
      my $ypos = $node->get_attribute("ypos");
      $max = $ypos if ($ypos > $ymax);
    }
    $diagram->set_attribute("xsize", $xmax);
    $diagram->set_attribute("ysize", $ymax);
    return ($xmax,$ymax);
  }

  sub calc_node_limits {
    my ($diagram, $node) = @_;
    $left = $diagram->calc_xpos($node) - $node->calc_label_width()/2;
    $right = $diagram->calc_xpos($node) + $node->calc_label_width()/2;
    $top = $diagram->calc_ypos($node) - $node->calc_label_height()/2;
    $bottom = $diagram->calc_ypos($node) + $node->calc_label_height()/2;
    return ($left, $top, $right, $bottom);
  }

  sub calc_limits {
    ### returns the rectangle encomprizing all node centers
    ### usage : 
    ###      my ($left,$top,$right,$bottom) = $diagram->calc_limits();
    my ($diagram) = @_;
    my $left = undef;
    my $top = undef;
    my $right = undef;
    my $bottom = undef;
    foreach my $node ($diagram->get_attribute("nodes")) {
      my ($l,$t,$r,$b) = $diagram->calc_node_limits($node);
      $left = $l if (!(defined($left)) || ($l < $left));
      $top = $t if (!(defined($top)) || ($t < $top));
      $right = $r if (!(defined($right)) || ($r > $right));
      $bottom = $b if (!(defined($bottom)) || ($b > $bottom));
    }
    return ($left, $top, $right, $bottom);
  }

  sub calc_xpos {
    ### calculate the absolute position of a node within a diagram,
    ### as a function of the other nodes it is related to
    my ($diagram, $node) = @_;
    my $xpos = undef;
    
    if ($node->get_attribute("xpos") =~ /(\S+):(.*)/) {
      my $ref_node_id = $1;
      my $x_offset = $2;
      my $ref_node = $diagram->get_object($ref_node_id);
      $xpos = $diagram->calc_xpos($ref_node) + $x_offset;
    } else {
      $xpos = $node->get_attribute("xpos");
    }
    return $xpos;
  }

  sub calc_ypos {
    ### calculate the absolute position of a node within a diagram,
    ### as a function of the other nodes it is related to
    my ($diagram, $node) = @_;
    my $ypos = undef;
    if ($node->get_attribute("ypos") =~ /(\S+):(.*)/) {
      my $ref_node = $diagram->get_object($1);
      my $y_offset = $2;
      $ypos = $diagram->calc_ypos($ref_node) + $y_offset;
    } else {
      $ypos = $node->get_attribute("ypos");
    }
    return $ypos;
  }

  sub print {
      ### usage $diagram->print($format,$file);
      ### if $file is left blank, print to the STDOUT
      ### supported formats : 
      ###    tdf   one line per node or arc
      ###    pml   coords on separate lines
      ### $absolute: convert all relative positions into absolute positions
      my ($diagram, $format, $file, $absolute) = @_;
      
      warn (";\n; ", &main::AlphaDate, 
	    " printing pathway diagram ", $diagram,
	    " in ", $format, " format ", 
	    " to file ", $file, 
	    "\n")
	  if ($main::verbose >= 1);
      
      ### open output stream
      if ($file) {
	  open STDOUT, "> $file" || 
	      die "ERROR: cannot write output file '$file'\n";
      } 
      
      ### print diagram elements
      
      ### tdf (textual diagram format)
      ### one row per node
      ### one row per arc
      ### tab-delimited fields
      ### can be used by graph_viewer.pl to generate a gif file
      if ($format eq "tdf") {
	  foreach $node ($diagram->get_attribute("nodes")) {
	      my $xpos = undef;
	      my $ypos = undef;
	      if ($absolute) {
		  $xpos = $diagram->calc_xpos($node);
	      } else {
		  $xpos = $node->get_attribute("xpos");
	      }	
	      if ($absolute) {
		  $ypos = $diagram->calc_ypos($node);
	      } else {
		  $ypos = $node->get_attribute("ypos");
	      }

	      print "node";
	      print "\t", $node->get_attribute("id");
	      print "\t", $xpos;
	      print "\t", $ypos;
	      print "\t", $node->get_attribute("label");
	      print "\t", $node->get_attribute("ref_class");
	      my $ref = $node->get_attribute("ref_db");
	      $ref .= "::".$node->get_attribute("ref_class");
	      $ref .= "::".$node->get_attribute("ref_id");
	      print "\t", lc($ref);
	      print "\n";
	  }
	  foreach $arc ($diagram->get_attribute("arcs")) {
	      print "arc";
	      print "\t", $arc->get_attribute("id");
	      print "\t", $arc->get_attribute("from");
	      print "\t", $arc->get_attribute("to");
	      print "\t", $arc->get_attribute("label");
	      print "\t", $arc->get_attribute("type");
	      print "\t", $arc->get_attribute("url");
	      print "\n";
	  }
	  
	  ### tdd (textual diagram description)
	  ### one row per node, without coordinate specification
	  ### one row per arc, without coordinate specification
	  ### one row for each coordinate specification
	  ### tab-delimited fields
	  ### can be used to generate a java view with Tom Sawyer library
      } elsif ($format eq "tdd") {
	  
	  print "; Diagram attributes\n";
	  print "ID\t", $diagram->get_attribute("id"), "\n";
	  print "NAME\t", $diagram->get_name(), "\n";
	  print "TYPE\t", $diagram->get_attribute("type"), "\n";
	  print "DESC\t", $diagram->get_attribute("description"), "\n";
	  print "; Node and Arc attributes\n";
	  print "; node description      NODE    id      label   type\n";
	  print "; node shift            SHIFT   nodeId  xShift  yShift\n";
	  print "; node anchoring        ANCHOR  nodeId  anchId  orient\n";
	  print "; node position         COORD   nodeId  xPos    yPos\n";
	  print "; arc description       ARC     arcId   fromId  toId\n";
	  print "; cross reference       XREF    DB      extID\n";

	  ### node definition
	  foreach $node ($diagram->get_attribute("nodes")) {
	      print "NODE";
	      print "\t", $node->get_attribute("id");
	      print "\t", $node->get_attribute("label");
	      print "\t", $node->get_attribute("type");
	      print "\n";
	  }
	  
	  ### node coordinates
	  foreach $node ($diagram->get_attribute("nodes")) {
	      my $xpos = undef;
	      my $ypos = undef;
	      if ($absolute) {
		  $xpos = $diagram->calc_xpos($node);
	      } else {
		  $xpos = $node->get_attribute("xpos");
	      }	
	      if ($absolute) {
		  $ypos = $diagram->calc_ypos($node);
	      } else {
		  $ypos = $node->get_attribute("ypos");
	      }
	      if ((defined($xpos)) &&
		  (defined($ypos))) {
		  print "SHIFT";
		  print "\t", $node->get_attribute("id");
		  printf "\t%f", $xpos;
		  printf "\t%f", $ypos;
		  print "\n";
	      }
	  }
	  
	  ### arcs
	  foreach $arc ($diagram->get_attribute("arcs")) {
	      print "ARC";
	      print "\t", $arc->get_attribute("id");
	      print "\t", $arc->get_attribute("from");
	      print "\t", $arc->get_attribute("to");
	      print "\t", $arc->get_attribute("type");
	      print "\n";
	  }
	  
	  ### reference to database objects
	  foreach $node ($diagram->get_attribute("nodes")) {
	      if ((my $db = $node->get_attribute("ref_db")) &&
#		  (my $class = $node->get_attribute("ref_class")) &&
		  (my $ref_id = $node->get_attribute("ref_id"))
		  ) { 
		  next if ($db eq $null);
		  next if ($ref_id eq $null);
		  print "XREF";
		  print "\t", $node->get_attribute("id");
		  print "\t", $db;
		  print "\t", $ref_id;
		  print "\n";
	      }
	  }
	  

      } else {
	  die "Error: format $format is not supported for pathway diagram printing\n";
      }
      
      close STDOUT if ($file);
  }
}


package classes::PathwayDiagramNode;
### the pathway diagram node must be of type :
### - BiochemicalEntity, 
### - BiochemicalActivity,
### - Pathway
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "pwdn_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",        ### internal ID for the pathway diagram
		      xpos=>"SCALAR",      ### absolute or relative ? To see later
		      ypos=>"SCALAR",
		      label=>"SCALAR",     ### string to display on the node in the diagram
		      ref_db=>"SCALAR",    ### database where the node comes from
		      ref_class=>"SCALAR",      ### this is the ref class
		      ref_id=>"SCALAR");    ### ID of the node in the source database

  sub calc_label_width {
    my ($node) = @_;
    my $label = $node->get_attribute("label");
    my $label_width = &GV::GraphUtil::calc_string_width($label);
    return $label_width;
  }

  sub calc_label_height {
    my ($node) = @_;
    my $label = $node->get_attribute("label");
    my $label_height = &GV::GraphUtil::calc_string_height($label);
    return $label_height;
  }
}

package classes::PathwayDiagramArc;
### the pathway diagram arc must be of type :
### - BiochemicalEntity, 
### - BiochemicalActivity,
### - Pathway
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "pwda_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",        ### internal ID for the pathway diagram
		      from=>"SCALAR",      ### absolute or relative ? To see later
		      to=>"SCALAR",
		      label=>"SCALAR",
		      type=>"SCALAR");     ### string to display on the arc in the diagram
}

return 1;
