package MyTypes::DyadAnalysisRequest;
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
my %verbosity_of :ATTR(:get<verbosity>);
my %sequence_of :ATTR(:get<sequence>);
my %tmp_infile_of :ATTR(:get<tmp_infile>);
my %format_of :ATTR(:get<format>);
my %length_of :ATTR(:get<length>);
my %spacing_of :ATTR(:get<spacing>);
my %organism_of :ATTR(:get<organism>);
my %background_of :ATTR(:get<background>);
my %stats_of :ATTR(:get<stats>);
my %type_of :ATTR(:get<type>);
my %noov_of :ATTR(:get<noov>);
my %str_of :ATTR(:get<str>);
my %sort_of :ATTR(:get<sort>);
my %under_of :ATTR(:get<under>);
my %two_tails_of :ATTR(:get<two_tails>);
my %zeroocc_of :ATTR(:get<zeroocc>);
my %lth_of :ATTR(:get<lth>);
my %uth_of :ATTR(:get<uth>);

__PACKAGE__->_factory(
    [ qw(        output
        verbosity
        sequence
        tmp_infile
        format
        length
        spacing
        organism
        background
        stats
        type
        noov
        str
        sort
        under
        two_tails
        zeroocc
        lth
        uth

    ) ],
    {
        'output' => \%output_of,
        'verbosity' => \%verbosity_of,
        'sequence' => \%sequence_of,
        'tmp_infile' => \%tmp_infile_of,
        'format' => \%format_of,
        'length' => \%length_of,
        'spacing' => \%spacing_of,
        'organism' => \%organism_of,
        'background' => \%background_of,
        'stats' => \%stats_of,
        'type' => \%type_of,
        'noov' => \%noov_of,
        'str' => \%str_of,
        'sort' => \%sort_of,
        'under' => \%under_of,
        'two_tails' => \%two_tails_of,
        'zeroocc' => \%zeroocc_of,
        'lth' => \%lth_of,
        'uth' => \%uth_of,
    },
    {
        'output' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'verbosity' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'sequence' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_infile' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'format' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'length' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'spacing' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'organism' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'background' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'stats' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'type' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'noov' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'str' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'sort' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'under' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'two_tails' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'zeroocc' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'lth' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'uth' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
    },
    {

        'output' => 'output',
        'verbosity' => 'verbosity',
        'sequence' => 'sequence',
        'tmp_infile' => 'tmp_infile',
        'format' => 'format',
        'length' => 'length',
        'spacing' => 'spacing',
        'organism' => 'organism',
        'background' => 'background',
        'stats' => 'stats',
        'type' => 'type',
        'noov' => 'noov',
        'str' => 'str',
        'sort' => 'sort',
        'under' => 'under',
        'two_tails' => 'two_tails',
        'zeroocc' => 'zeroocc',
        'lth' => 'lth',
        'uth' => 'uth',
    }
);

} # end BLOCK







1;


=pod

=head1 NAME

MyTypes::DyadAnalysisRequest

=head1 DESCRIPTION

Perl data type class for the XML Schema defined complexType
DyadAnalysisRequest from the namespace urn:RSATWS.

Parameters for the operation dyad_analysis.




=head2 PROPERTIES

The following properties may be accessed using get_PROPERTY / set_PROPERTY
methods:

=over

=item * output


=item * verbosity


=item * sequence


=item * tmp_infile


=item * format


=item * length


=item * spacing


=item * organism


=item * background


=item * stats


=item * type


=item * noov


=item * str


=item * sort


=item * under


=item * two_tails


=item * zeroocc


=item * lth


=item * uth




=back


=head1 METHODS

=head2 new

Constructor. The following data structure may be passed to new():

 { # MyTypes::DyadAnalysisRequest
   output =>  $some_value, # string
   verbosity =>  $some_value, # int
   sequence =>  $some_value, # string
   tmp_infile =>  $some_value, # string
   format =>  $some_value, # string
   length =>  $some_value, # int
   spacing =>  $some_value, # string
   organism =>  $some_value, # string
   background =>  $some_value, # string
   stats =>  $some_value, # string
   type =>  $some_value, # string
   noov =>  $some_value, # int
   str =>  $some_value, # int
   sort =>  $some_value, # int
   under =>  $some_value, # int
   two_tails =>  $some_value, # int
   zeroocc =>  $some_value, # int
   lth =>  $some_value, # string
   uth =>  $some_value, # string
 },




=head1 AUTHOR

Generated by SOAP::WSDL

=cut

