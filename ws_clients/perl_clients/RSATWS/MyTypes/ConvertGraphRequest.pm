package MyTypes::ConvertGraphRequest;
use strict;
use warnings;


our $XML_ATTRIBUTE_CLASS;
undef $XML_ATTRIBUTE_CLASS;

sub __get_attr_class {
    return $XML_ATTRIBUTE_CLASS;
}

use Class::Std::Fast::Storable constructor => 'none';
use base qw(SOAP::WSDL::XSD::Typelib::ComplexType);

Class::Std::initialize();

{ # BLOCK to scope variables

my %output_of :ATTR(:get<output>);
my %informat_of :ATTR(:get<informat>);
my %ecolors_of :ATTR(:get<ecolors>);
my %outformat_of :ATTR(:get<outformat>);
my %inputgraph_of :ATTR(:get<inputgraph>);
my %wcol_of :ATTR(:get<wcol>);
my %scol_of :ATTR(:get<scol>);
my %tcol_of :ATTR(:get<tcol>);
my %eccol_of :ATTR(:get<eccol>);
my %sccol_of :ATTR(:get<sccol>);
my %tccol_of :ATTR(:get<tccol>);
my %pathcol_of :ATTR(:get<pathcol>);
my %undirected_of :ATTR(:get<undirected>);
my %distinct_path_of :ATTR(:get<distinct_path>);
my %layout_of :ATTR(:get<layout>);
my %ewidth_of :ATTR(:get<ewidth>);
my %target_xpos_col_of :ATTR(:get<target_xpos_col>);
my %target_ypos_col_of :ATTR(:get<target_ypos_col>);
my %source_xpos_col_of :ATTR(:get<source_xpos_col>);
my %source_ypos_col_of :ATTR(:get<source_ypos_col>);

__PACKAGE__->_factory(
    [ qw(        output
        informat
        ecolors
        outformat
        inputgraph
        wcol
        scol
        tcol
        eccol
        sccol
        tccol
        pathcol
        undirected
        distinct_path
        layout
        ewidth
        target_xpos_col
        target_ypos_col
        source_xpos_col
        source_ypos_col

    ) ],
    {
        'output' => \%output_of,
        'informat' => \%informat_of,
        'ecolors' => \%ecolors_of,
        'outformat' => \%outformat_of,
        'inputgraph' => \%inputgraph_of,
        'wcol' => \%wcol_of,
        'scol' => \%scol_of,
        'tcol' => \%tcol_of,
        'eccol' => \%eccol_of,
        'sccol' => \%sccol_of,
        'tccol' => \%tccol_of,
        'pathcol' => \%pathcol_of,
        'undirected' => \%undirected_of,
        'distinct_path' => \%distinct_path_of,
        'layout' => \%layout_of,
        'ewidth' => \%ewidth_of,
        'target_xpos_col' => \%target_xpos_col_of,
        'target_ypos_col' => \%target_ypos_col_of,
        'source_xpos_col' => \%source_xpos_col_of,
        'source_ypos_col' => \%source_ypos_col_of,
    },
    {
        'output' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'informat' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'ecolors' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'outformat' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'inputgraph' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'wcol' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'scol' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'tcol' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'eccol' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'sccol' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'tccol' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'pathcol' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'undirected' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'distinct_path' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'layout' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'ewidth' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'target_xpos_col' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'target_ypos_col' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'source_xpos_col' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'source_ypos_col' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
    },
    {

        'output' => 'output',
        'informat' => 'informat',
        'ecolors' => 'ecolors',
        'outformat' => 'outformat',
        'inputgraph' => 'inputgraph',
        'wcol' => 'wcol',
        'scol' => 'scol',
        'tcol' => 'tcol',
        'eccol' => 'eccol',
        'sccol' => 'sccol',
        'tccol' => 'tccol',
        'pathcol' => 'pathcol',
        'undirected' => 'undirected',
        'distinct_path' => 'distinct_path',
        'layout' => 'layout',
        'ewidth' => 'ewidth',
        'target_xpos_col' => 'target_xpos_col',
        'target_ypos_col' => 'target_ypos_col',
        'source_xpos_col' => 'source_xpos_col',
        'source_ypos_col' => 'source_ypos_col',
    }
);

} # end BLOCK







1;


=pod

=head1 NAME

MyTypes::ConvertGraphRequest

=head1 DESCRIPTION

Perl data type class for the XML Schema defined complexType
ConvertGraphRequest from the namespace urn:RSATWS.






=head2 PROPERTIES

The following properties may be accessed using get_PROPERTY / set_PROPERTY
methods:

=over

=item * output


=item * informat


=item * ecolors


=item * outformat


=item * inputgraph


=item * wcol


=item * scol


=item * tcol


=item * eccol


=item * sccol


=item * tccol


=item * pathcol


=item * undirected


=item * distinct_path


=item * layout


=item * ewidth


=item * target_xpos_col


=item * target_ypos_col


=item * source_xpos_col


=item * source_ypos_col




=back


=head1 METHODS

=head2 new

Constructor. The following data structure may be passed to new():

 { # MyTypes::ConvertGraphRequest
   output =>  $some_value, # string
   informat =>  $some_value, # string
   ecolors =>  $some_value, # string
   outformat =>  $some_value, # string
   inputgraph =>  $some_value, # string
   wcol =>  $some_value, # int
   scol =>  $some_value, # int
   tcol =>  $some_value, # int
   eccol =>  $some_value, # int
   sccol =>  $some_value, # int
   tccol =>  $some_value, # int
   pathcol =>  $some_value, # int
   undirected =>  $some_value, # int
   distinct_path =>  $some_value, # int
   layout =>  $some_value, # int
   ewidth =>  $some_value, # int
   target_xpos_col =>  $some_value, # int
   target_ypos_col =>  $some_value, # int
   source_xpos_col =>  $some_value, # int
   source_ypos_col =>  $some_value, # int
 },




=head1 AUTHOR

Generated by SOAP::WSDL

=cut

