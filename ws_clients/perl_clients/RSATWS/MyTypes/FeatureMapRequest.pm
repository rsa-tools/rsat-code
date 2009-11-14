package MyTypes::FeatureMapRequest;
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
my %features_of :ATTR(:get<features>);
my %tmp_infile_of :ATTR(:get<tmp_infile>);
my %sequence_of :ATTR(:get<sequence>);
my %tmp_sequence_file_of :ATTR(:get<tmp_sequence_file>);
my %sequence_format_of :ATTR(:get<sequence_format>);
my %format_of :ATTR(:get<format>);
my %from_of :ATTR(:get<from>);
my %to_of :ATTR(:get<to>);
my %title_of :ATTR(:get<title>);
my %label_of :ATTR(:get<label>);
my %symbol_of :ATTR(:get<symbol>);
my %dot_of :ATTR(:get<dot>);
my %mlen_of :ATTR(:get<mlen>);
my %mapthick_of :ATTR(:get<mapthick>);
my %mspacing_of :ATTR(:get<mspacing>);
my %origin_of :ATTR(:get<origin>);
my %legend_of :ATTR(:get<legend>);
my %scalebar_of :ATTR(:get<scalebar>);
my %scalestep_of :ATTR(:get<scalestep>);
my %scorethick_of :ATTR(:get<scorethick>);
my %maxscore_of :ATTR(:get<maxscore>);
my %minscore_of :ATTR(:get<minscore>);
my %maxfthick_of :ATTR(:get<maxfthick>);
my %minfthick_of :ATTR(:get<minfthick>);
my %htmap_of :ATTR(:get<htmap>);
my %mono_of :ATTR(:get<mono>);
my %orientation_of :ATTR(:get<orientation>);
my %select_of :ATTR(:get<select>);

__PACKAGE__->_factory(
    [ qw(        output
        features
        tmp_infile
        sequence
        tmp_sequence_file
        sequence_format
        format
        from
        to
        title
        label
        symbol
        dot
        mlen
        mapthick
        mspacing
        origin
        legend
        scalebar
        scalestep
        scorethick
        maxscore
        minscore
        maxfthick
        minfthick
        htmap
        mono
        orientation
        select

    ) ],
    {
        'output' => \%output_of,
        'features' => \%features_of,
        'tmp_infile' => \%tmp_infile_of,
        'sequence' => \%sequence_of,
        'tmp_sequence_file' => \%tmp_sequence_file_of,
        'sequence_format' => \%sequence_format_of,
        'format' => \%format_of,
        'from' => \%from_of,
        'to' => \%to_of,
        'title' => \%title_of,
        'label' => \%label_of,
        'symbol' => \%symbol_of,
        'dot' => \%dot_of,
        'mlen' => \%mlen_of,
        'mapthick' => \%mapthick_of,
        'mspacing' => \%mspacing_of,
        'origin' => \%origin_of,
        'legend' => \%legend_of,
        'scalebar' => \%scalebar_of,
        'scalestep' => \%scalestep_of,
        'scorethick' => \%scorethick_of,
        'maxscore' => \%maxscore_of,
        'minscore' => \%minscore_of,
        'maxfthick' => \%maxfthick_of,
        'minfthick' => \%minfthick_of,
        'htmap' => \%htmap_of,
        'mono' => \%mono_of,
        'orientation' => \%orientation_of,
        'select' => \%select_of,
    },
    {
        'output' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'features' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_infile' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'sequence' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_sequence_file' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'sequence_format' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'format' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'from' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'to' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'title' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'label' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'symbol' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'dot' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'mlen' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'mapthick' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'mspacing' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'origin' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'legend' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'scalebar' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'scalestep' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'scorethick' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'maxscore' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'minscore' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'maxfthick' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'minfthick' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'htmap' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'mono' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'orientation' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'select' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
    },
    {

        'output' => 'output',
        'features' => 'features',
        'tmp_infile' => 'tmp_infile',
        'sequence' => 'sequence',
        'tmp_sequence_file' => 'tmp_sequence_file',
        'sequence_format' => 'sequence_format',
        'format' => 'format',
        'from' => 'from',
        'to' => 'to',
        'title' => 'title',
        'label' => 'label',
        'symbol' => 'symbol',
        'dot' => 'dot',
        'mlen' => 'mlen',
        'mapthick' => 'mapthick',
        'mspacing' => 'mspacing',
        'origin' => 'origin',
        'legend' => 'legend',
        'scalebar' => 'scalebar',
        'scalestep' => 'scalestep',
        'scorethick' => 'scorethick',
        'maxscore' => 'maxscore',
        'minscore' => 'minscore',
        'maxfthick' => 'maxfthick',
        'minfthick' => 'minfthick',
        'htmap' => 'htmap',
        'mono' => 'mono',
        'orientation' => 'orientation',
        'select' => 'select',
    }
);

} # end BLOCK







1;


=pod

=head1 NAME

MyTypes::FeatureMapRequest

=head1 DESCRIPTION

Perl data type class for the XML Schema defined complexType
FeatureMapRequest from the namespace urn:RSATWS.

Parameters for the operation feature_map.




=head2 PROPERTIES

The following properties may be accessed using get_PROPERTY / set_PROPERTY
methods:

=over

=item * output


=item * features


=item * tmp_infile


=item * sequence


=item * tmp_sequence_file


=item * sequence_format


=item * format


=item * from


=item * to


=item * title


=item * label


=item * symbol


=item * dot


=item * mlen


=item * mapthick


=item * mspacing


=item * origin


=item * legend


=item * scalebar


=item * scalestep


=item * scorethick


=item * maxscore


=item * minscore


=item * maxfthick


=item * minfthick


=item * htmap


=item * mono


=item * orientation


=item * select




=back


=head1 METHODS

=head2 new

Constructor. The following data structure may be passed to new():

 { # MyTypes::FeatureMapRequest
   output =>  $some_value, # string
   features =>  $some_value, # string
   tmp_infile =>  $some_value, # string
   sequence =>  $some_value, # string
   tmp_sequence_file =>  $some_value, # string
   sequence_format =>  $some_value, # string
   format =>  $some_value, # string
   from =>  $some_value, # int
   to =>  $some_value, # int
   title =>  $some_value, # string
   label =>  $some_value, # string
   symbol =>  $some_value, # int
   dot =>  $some_value, # int
   mlen =>  $some_value, # int
   mapthick =>  $some_value, # int
   mspacing =>  $some_value, # int
   origin =>  $some_value, # int
   legend =>  $some_value, # int
   scalebar =>  $some_value, # int
   scalestep =>  $some_value, # int
   scorethick =>  $some_value, # int
   maxscore =>  $some_value, # int
   minscore =>  $some_value, # int
   maxfthick =>  $some_value, # int
   minfthick =>  $some_value, # int
   htmap =>  $some_value, # int
   mono =>  $some_value, # int
   orientation =>  $some_value, # string
   select =>  $some_value, # string
 },




=head1 AUTHOR

Generated by SOAP::WSDL

=cut

