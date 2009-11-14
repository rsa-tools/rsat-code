package MyTypes::RetrieveSequenceMultigenomeRequest;
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
my %input_of :ATTR(:get<input>);
my %tmp_input_file_of :ATTR(:get<tmp_input_file>);
my %all_of :ATTR(:get<all>);
my %noorf_of :ATTR(:get<noorf>);
my %from_of :ATTR(:get<from>);
my %to_of :ATTR(:get<to>);
my %feattype_of :ATTR(:get<feattype>);
my %type_of :ATTR(:get<type>);
my %format_of :ATTR(:get<format>);
my %lw_of :ATTR(:get<lw>);
my %label_of :ATTR(:get<label>);
my %label_sep_of :ATTR(:get<label_sep>);
my %nocom_of :ATTR(:get<nocom>);
my %repeat_of :ATTR(:get<repeat>);
my %imp_pos_of :ATTR(:get<imp_pos>);
my %gene_col_of :ATTR(:get<gene_col>);
my %org_col_of :ATTR(:get<org_col>);

__PACKAGE__->_factory(
    [ qw(        output
        input
        tmp_input_file
        all
        noorf
        from
        to
        feattype
        type
        format
        lw
        label
        label_sep
        nocom
        repeat
        imp_pos
        gene_col
        org_col

    ) ],
    {
        'output' => \%output_of,
        'input' => \%input_of,
        'tmp_input_file' => \%tmp_input_file_of,
        'all' => \%all_of,
        'noorf' => \%noorf_of,
        'from' => \%from_of,
        'to' => \%to_of,
        'feattype' => \%feattype_of,
        'type' => \%type_of,
        'format' => \%format_of,
        'lw' => \%lw_of,
        'label' => \%label_of,
        'label_sep' => \%label_sep_of,
        'nocom' => \%nocom_of,
        'repeat' => \%repeat_of,
        'imp_pos' => \%imp_pos_of,
        'gene_col' => \%gene_col_of,
        'org_col' => \%org_col_of,
    },
    {
        'output' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'input' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_input_file' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'all' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'noorf' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'from' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'to' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'feattype' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'type' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'format' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'lw' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'label' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'label_sep' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'nocom' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'repeat' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'imp_pos' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'gene_col' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'org_col' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
    },
    {

        'output' => 'output',
        'input' => 'input',
        'tmp_input_file' => 'tmp_input_file',
        'all' => 'all',
        'noorf' => 'noorf',
        'from' => 'from',
        'to' => 'to',
        'feattype' => 'feattype',
        'type' => 'type',
        'format' => 'format',
        'lw' => 'lw',
        'label' => 'label',
        'label_sep' => 'label_sep',
        'nocom' => 'nocom',
        'repeat' => 'repeat',
        'imp_pos' => 'imp_pos',
        'gene_col' => 'gene_col',
        'org_col' => 'org_col',
    }
);

} # end BLOCK







1;


=pod

=head1 NAME

MyTypes::RetrieveSequenceMultigenomeRequest

=head1 DESCRIPTION

Perl data type class for the XML Schema defined complexType
RetrieveSequenceMultigenomeRequest from the namespace urn:RSATWS.

Parameters for the operation retrieve_seq_multigenome.




=head2 PROPERTIES

The following properties may be accessed using get_PROPERTY / set_PROPERTY
methods:

=over

=item * output


=item * input


=item * tmp_input_file


=item * all


=item * noorf


=item * from


=item * to


=item * feattype


=item * type


=item * format


=item * lw


=item * label


=item * label_sep


=item * nocom


=item * repeat


=item * imp_pos


=item * gene_col


=item * org_col




=back


=head1 METHODS

=head2 new

Constructor. The following data structure may be passed to new():

 { # MyTypes::RetrieveSequenceMultigenomeRequest
   output =>  $some_value, # string
   input =>  $some_value, # string
   tmp_input_file =>  $some_value, # string
   all =>  $some_value, # int
   noorf =>  $some_value, # int
   from =>  $some_value, # int
   to =>  $some_value, # int
   feattype =>  $some_value, # string
   type =>  $some_value, # string
   format =>  $some_value, # string
   lw =>  $some_value, # int
   label =>  $some_value, # string
   label_sep =>  $some_value, # string
   nocom =>  $some_value, # int
   repeat =>  $some_value, # int
   imp_pos =>  $some_value, # int
   gene_col =>  $some_value, # int
   org_col =>  $some_value, # int
 },




=head1 AUTHOR

Generated by SOAP::WSDL

=cut

