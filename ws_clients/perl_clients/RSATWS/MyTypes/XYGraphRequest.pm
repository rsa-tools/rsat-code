package MyTypes::XYGraphRequest;
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
my %inputFile_of :ATTR(:get<inputFile>);
my %format_of :ATTR(:get<format>);
my %title1_of :ATTR(:get<title1>);
my %title2_of :ATTR(:get<title2>);
my %lines_of :ATTR(:get<lines>);
my %legend_of :ATTR(:get<legend>);
my %header_of :ATTR(:get<header>);
my %xleg1_of :ATTR(:get<xleg1>);
my %xleg2_of :ATTR(:get<xleg2>);
my %yleg1_of :ATTR(:get<yleg1>);
my %yleg2_of :ATTR(:get<yleg2>);
my %xmax_of :ATTR(:get<xmax>);
my %ymax_of :ATTR(:get<ymax>);
my %xmin_of :ATTR(:get<xmin>);
my %ymin_of :ATTR(:get<ymin>);
my %ylog_of :ATTR(:get<ylog>);
my %xlog_of :ATTR(:get<xlog>);
my %xcol_of :ATTR(:get<xcol>);
my %ycol_of :ATTR(:get<ycol>);

__PACKAGE__->_factory(
    [ qw(        output
        inputFile
        format
        title1
        title2
        lines
        legend
        header
        xleg1
        xleg2
        yleg1
        yleg2
        xmax
        ymax
        xmin
        ymin
        ylog
        xlog
        xcol
        ycol

    ) ],
    {
        'output' => \%output_of,
        'inputFile' => \%inputFile_of,
        'format' => \%format_of,
        'title1' => \%title1_of,
        'title2' => \%title2_of,
        'lines' => \%lines_of,
        'legend' => \%legend_of,
        'header' => \%header_of,
        'xleg1' => \%xleg1_of,
        'xleg2' => \%xleg2_of,
        'yleg1' => \%yleg1_of,
        'yleg2' => \%yleg2_of,
        'xmax' => \%xmax_of,
        'ymax' => \%ymax_of,
        'xmin' => \%xmin_of,
        'ymin' => \%ymin_of,
        'ylog' => \%ylog_of,
        'xlog' => \%xlog_of,
        'xcol' => \%xcol_of,
        'ycol' => \%ycol_of,
    },
    {
        'output' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'inputFile' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'format' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'title1' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'title2' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'lines' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'legend' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'header' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'xleg1' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'xleg2' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'yleg1' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'yleg2' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'xmax' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'ymax' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'xmin' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'ymin' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'ylog' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'xlog' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'xcol' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'ycol' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
    },
    {

        'output' => 'output',
        'inputFile' => 'inputFile',
        'format' => 'format',
        'title1' => 'title1',
        'title2' => 'title2',
        'lines' => 'lines',
        'legend' => 'legend',
        'header' => 'header',
        'xleg1' => 'xleg1',
        'xleg2' => 'xleg2',
        'yleg1' => 'yleg1',
        'yleg2' => 'yleg2',
        'xmax' => 'xmax',
        'ymax' => 'ymax',
        'xmin' => 'xmin',
        'ymin' => 'ymin',
        'ylog' => 'ylog',
        'xlog' => 'xlog',
        'xcol' => 'xcol',
        'ycol' => 'ycol',
    }
);

} # end BLOCK







1;


=pod

=head1 NAME

MyTypes::XYGraphRequest

=head1 DESCRIPTION

Perl data type class for the XML Schema defined complexType
XYGraphRequest from the namespace urn:RSATWS.






=head2 PROPERTIES

The following properties may be accessed using get_PROPERTY / set_PROPERTY
methods:

=over

=item * output


=item * inputFile


=item * format


=item * title1


=item * title2


=item * lines


=item * legend


=item * header


=item * xleg1


=item * xleg2


=item * yleg1


=item * yleg2


=item * xmax


=item * ymax


=item * xmin


=item * ymin


=item * ylog


=item * xlog


=item * xcol


=item * ycol




=back


=head1 METHODS

=head2 new

Constructor. The following data structure may be passed to new():

 { # MyTypes::XYGraphRequest
   output =>  $some_value, # string
   inputFile =>  $some_value, # string
   format =>  $some_value, # string
   title1 =>  $some_value, # string
   title2 =>  $some_value, # string
   lines =>  $some_value, # int
   legend =>  $some_value, # int
   header =>  $some_value, # int
   xleg1 =>  $some_value, # string
   xleg2 =>  $some_value, # string
   yleg1 =>  $some_value, # string
   yleg2 =>  $some_value, # string
   xmax =>  $some_value, # string
   ymax =>  $some_value, # string
   xmin =>  $some_value, # string
   ymin =>  $some_value, # string
   ylog =>  $some_value, # string
   xlog =>  $some_value, # string
   xcol =>  $some_value, # string
   ycol =>  $some_value, # string
 },




=head1 AUTHOR

Generated by SOAP::WSDL

=cut

