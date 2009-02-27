package MyTypes::RandomGraphRequest;
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
my %outformat_of :ATTR(:get<outformat>);
my %inputgraph_of :ATTR(:get<inputgraph>);
my %random_type_of :ATTR(:get<random_type>);
my %wcol_of :ATTR(:get<wcol>);
my %scol_of :ATTR(:get<scol>);
my %tcol_of :ATTR(:get<tcol>);
my %edges_of :ATTR(:get<edges>);
my %degree_of :ATTR(:get<degree>);
my %nodes_of :ATTR(:get<nodes>);
my %mean_of :ATTR(:get<mean>);
my %sd_of :ATTR(:get<sd>);
my %directed_of :ATTR(:get<directed>);
my %no_single_of :ATTR(:get<no_single>);
my %duplicate_of :ATTR(:get<duplicate>);
my %col_conservation_of :ATTR(:get<col_conservation>);
my %normal_of :ATTR(:get<normal>);

__PACKAGE__->_factory(
    [ qw(        output
        informat
        outformat
        inputgraph
        random_type
        wcol
        scol
        tcol
        edges
        degree
        nodes
        mean
        sd
        directed
        no_single
        duplicate
        col_conservation
        normal

    ) ],
    {
        'output' => \%output_of,
        'informat' => \%informat_of,
        'outformat' => \%outformat_of,
        'inputgraph' => \%inputgraph_of,
        'random_type' => \%random_type_of,
        'wcol' => \%wcol_of,
        'scol' => \%scol_of,
        'tcol' => \%tcol_of,
        'edges' => \%edges_of,
        'degree' => \%degree_of,
        'nodes' => \%nodes_of,
        'mean' => \%mean_of,
        'sd' => \%sd_of,
        'directed' => \%directed_of,
        'no_single' => \%no_single_of,
        'duplicate' => \%duplicate_of,
        'col_conservation' => \%col_conservation_of,
        'normal' => \%normal_of,
    },
    {
        'output' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'informat' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'outformat' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'inputgraph' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'random_type' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'wcol' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'scol' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'tcol' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'edges' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'degree' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'nodes' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'mean' => 'SOAP::WSDL::XSD::Typelib::Builtin::float',
        'sd' => 'SOAP::WSDL::XSD::Typelib::Builtin::float',
        'directed' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'no_single' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'duplicate' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'col_conservation' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'normal' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
    },
    {

        'output' => 'output',
        'informat' => 'informat',
        'outformat' => 'outformat',
        'inputgraph' => 'inputgraph',
        'random_type' => 'random_type',
        'wcol' => 'wcol',
        'scol' => 'scol',
        'tcol' => 'tcol',
        'edges' => 'edges',
        'degree' => 'degree',
        'nodes' => 'nodes',
        'mean' => 'mean',
        'sd' => 'sd',
        'directed' => 'directed',
        'no_single' => 'no_single',
        'duplicate' => 'duplicate',
        'col_conservation' => 'col_conservation',
        'normal' => 'normal',
    }
);

} # end BLOCK







1;


=pod

=head1 NAME

MyTypes::RandomGraphRequest

=head1 DESCRIPTION

Perl data type class for the XML Schema defined complexType
RandomGraphRequest from the namespace urn:RSATWS.






=head2 PROPERTIES

The following properties may be accessed using get_PROPERTY / set_PROPERTY
methods:

=over

=item * output


=item * informat


=item * outformat


=item * inputgraph


=item * random_type


=item * wcol


=item * scol


=item * tcol


=item * edges


=item * degree


=item * nodes


=item * mean


=item * sd


=item * directed


=item * no_single


=item * duplicate


=item * col_conservation


=item * normal




=back


=head1 METHODS

=head2 new

Constructor. The following data structure may be passed to new():

 { # MyTypes::RandomGraphRequest
   output =>  $some_value, # string
   informat =>  $some_value, # string
   outformat =>  $some_value, # string
   inputgraph =>  $some_value, # string
   random_type =>  $some_value, # string
   wcol =>  $some_value, # int
   scol =>  $some_value, # int
   tcol =>  $some_value, # int
   edges =>  $some_value, # int
   degree =>  $some_value, # int
   nodes =>  $some_value, # int
   mean =>  $some_value, # float
   sd =>  $some_value, # float
   directed =>  $some_value, # int
   no_single =>  $some_value, # int
   duplicate =>  $some_value, # int
   col_conservation =>  $some_value, # int
   normal =>  $some_value, # int
 },




=head1 AUTHOR

Generated by SOAP::WSDL

=cut

