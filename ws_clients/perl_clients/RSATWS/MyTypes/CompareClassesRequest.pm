package MyTypes::CompareClassesRequest;
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
my %ref_classes_of :ATTR(:get<ref_classes>);
my %query_classes_of :ATTR(:get<query_classes>);
my %return_fields_of :ATTR(:get<return_fields>);
my %score_column_of :ATTR(:get<score_column>);
my %input_classes_of :ATTR(:get<input_classes>);
my %upper_threshold_field_of :ATTR(:get<upper_threshold_field>);
my %upper_threshold_value_of :ATTR(:get<upper_threshold_value>);
my %lower_threshold_field_of :ATTR(:get<lower_threshold_field>);
my %lower_threshold_value_of :ATTR(:get<lower_threshold_value>);
my %sort_of :ATTR(:get<sort>);
my %distinct_of :ATTR(:get<distinct>);
my %triangle_of :ATTR(:get<triangle>);
my %matrix_of :ATTR(:get<matrix>);

__PACKAGE__->_factory(
    [ qw(        output
        ref_classes
        query_classes
        return_fields
        score_column
        input_classes
        upper_threshold_field
        upper_threshold_value
        lower_threshold_field
        lower_threshold_value
        sort
        distinct
        triangle
        matrix

    ) ],
    {
        'output' => \%output_of,
        'ref_classes' => \%ref_classes_of,
        'query_classes' => \%query_classes_of,
        'return_fields' => \%return_fields_of,
        'score_column' => \%score_column_of,
        'input_classes' => \%input_classes_of,
        'upper_threshold_field' => \%upper_threshold_field_of,
        'upper_threshold_value' => \%upper_threshold_value_of,
        'lower_threshold_field' => \%lower_threshold_field_of,
        'lower_threshold_value' => \%lower_threshold_value_of,
        'sort' => \%sort_of,
        'distinct' => \%distinct_of,
        'triangle' => \%triangle_of,
        'matrix' => \%matrix_of,
    },
    {
        'output' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'ref_classes' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'query_classes' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'return_fields' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'score_column' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'input_classes' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'upper_threshold_field' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'upper_threshold_value' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'lower_threshold_field' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'lower_threshold_value' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'sort' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'distinct' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'triangle' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'matrix' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
    },
    {

        'output' => 'output',
        'ref_classes' => 'ref_classes',
        'query_classes' => 'query_classes',
        'return_fields' => 'return_fields',
        'score_column' => 'score_column',
        'input_classes' => 'input_classes',
        'upper_threshold_field' => 'upper_threshold_field',
        'upper_threshold_value' => 'upper_threshold_value',
        'lower_threshold_field' => 'lower_threshold_field',
        'lower_threshold_value' => 'lower_threshold_value',
        'sort' => 'sort',
        'distinct' => 'distinct',
        'triangle' => 'triangle',
        'matrix' => 'matrix',
    }
);

} # end BLOCK







1;


=pod

=head1 NAME

MyTypes::CompareClassesRequest

=head1 DESCRIPTION

Perl data type class for the XML Schema defined complexType
CompareClassesRequest from the namespace urn:RSATWS.

Parameters for the operation compare_classes.




=head2 PROPERTIES

The following properties may be accessed using get_PROPERTY / set_PROPERTY
methods:

=over

=item * output


=item * ref_classes


=item * query_classes


=item * return_fields


=item * score_column


=item * input_classes


=item * upper_threshold_field


=item * upper_threshold_value


=item * lower_threshold_field


=item * lower_threshold_value


=item * sort


=item * distinct


=item * triangle


=item * matrix




=back


=head1 METHODS

=head2 new

Constructor. The following data structure may be passed to new():

 { # MyTypes::CompareClassesRequest
   output =>  $some_value, # string
   ref_classes =>  $some_value, # string
   query_classes =>  $some_value, # string
   return_fields =>  $some_value, # string
   score_column =>  $some_value, # int
   input_classes =>  $some_value, # string
   upper_threshold_field =>  $some_value, # string
   upper_threshold_value =>  $some_value, # string
   lower_threshold_field =>  $some_value, # string
   lower_threshold_value =>  $some_value, # string
   sort =>  $some_value, # string
   distinct =>  $some_value, # int
   triangle =>  $some_value, # int
   matrix =>  $some_value, # string
 },




=head1 AUTHOR

Generated by SOAP::WSDL

=cut

