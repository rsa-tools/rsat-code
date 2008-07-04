package MyTypes::RetrieveEnsemblSequenceRequest;
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
my %organism_of :ATTR(:get<organism>);
my %ensembl_host_of :ATTR(:get<ensembl_host>);
my %db_name_of :ATTR(:get<db_name>);
my %query_of :ATTR(:get<query>);
my %all_of :ATTR(:get<all>);
my %noorf_of :ATTR(:get<noorf>);
my %nogene_of :ATTR(:get<nogene>);
my %from_of :ATTR(:get<from>);
my %to_of :ATTR(:get<to>);
my %feattype_of :ATTR(:get<feattype>);
my %type_of :ATTR(:get<type>);
my %chromosome_of :ATTR(:get<chromosome>);
my %left_of :ATTR(:get<left>);
my %right_of :ATTR(:get<right>);
my %strand_of :ATTR(:get<strand>);
my %features_of :ATTR(:get<features>);
my %feat_format_of :ATTR(:get<feat_format>);
my %mask_coding_of :ATTR(:get<mask_coding>);
my %repeat_of :ATTR(:get<repeat>);
my %all_transcripts_of :ATTR(:get<all_transcripts>);
my %first_intron_of :ATTR(:get<first_intron>);
my %non_coding_of :ATTR(:get<non_coding>);

__PACKAGE__->_factory(
    [ qw(        output
        organism
        ensembl_host
        db_name
        query
        all
        noorf
        nogene
        from
        to
        feattype
        type
        chromosome
        left
        right
        strand
        features
        feat_format
        mask_coding
        repeat
        all_transcripts
        first_intron
        non_coding

    ) ],
    {
        'output' => \%output_of,
        'organism' => \%organism_of,
        'ensembl_host' => \%ensembl_host_of,
        'db_name' => \%db_name_of,
        'query' => \%query_of,
        'all' => \%all_of,
        'noorf' => \%noorf_of,
        'nogene' => \%nogene_of,
        'from' => \%from_of,
        'to' => \%to_of,
        'feattype' => \%feattype_of,
        'type' => \%type_of,
        'chromosome' => \%chromosome_of,
        'left' => \%left_of,
        'right' => \%right_of,
        'strand' => \%strand_of,
        'features' => \%features_of,
        'feat_format' => \%feat_format_of,
        'mask_coding' => \%mask_coding_of,
        'repeat' => \%repeat_of,
        'all_transcripts' => \%all_transcripts_of,
        'first_intron' => \%first_intron_of,
        'non_coding' => \%non_coding_of,
    },
    {
        'output' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'organism' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'ensembl_host' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'db_name' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'query' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'all' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'noorf' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'nogene' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'from' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'to' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'feattype' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'type' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'chromosome' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'left' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'right' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'strand' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'features' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'feat_format' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'mask_coding' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'repeat' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'all_transcripts' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'first_intron' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'non_coding' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
    },
    {

        'output' => 'output',
        'organism' => 'organism',
        'ensembl_host' => 'ensembl_host',
        'db_name' => 'db_name',
        'query' => 'query',
        'all' => 'all',
        'noorf' => 'noorf',
        'nogene' => 'nogene',
        'from' => 'from',
        'to' => 'to',
        'feattype' => 'feattype',
        'type' => 'type',
        'chromosome' => 'chromosome',
        'left' => 'left',
        'right' => 'right',
        'strand' => 'strand',
        'features' => 'features',
        'feat_format' => 'feat_format',
        'mask_coding' => 'mask_coding',
        'repeat' => 'repeat',
        'all_transcripts' => 'all_transcripts',
        'first_intron' => 'first_intron',
        'non_coding' => 'non_coding',
    }
);

} # end BLOCK







1;


=pod

=head1 NAME

MyTypes::RetrieveEnsemblSequenceRequest

=head1 DESCRIPTION

Perl data type class for the XML Schema defined complexType
RetrieveEnsemblSequenceRequest from the namespace urn:RSATWS.

Parameters for the operation retrieve_ensembl_seq.




=head2 PROPERTIES

The following properties may be accessed using get_PROPERTY / set_PROPERTY
methods:

=over

=item * output


=item * organism


=item * ensembl_host


=item * db_name


=item * query


=item * all


=item * noorf


=item * nogene


=item * from


=item * to


=item * feattype


=item * type


=item * chromosome


=item * left


=item * right


=item * strand


=item * features


=item * feat_format


=item * mask_coding


=item * repeat


=item * all_transcripts


=item * first_intron


=item * non_coding




=back


=head1 METHODS

=head2 new

Constructor. The following data structure may be passed to new():

 { # MyTypes::RetrieveEnsemblSequenceRequest
   output =>  $some_value, # string
   organism =>  $some_value, # string
   ensembl_host =>  $some_value, # string
   db_name =>  $some_value, # string
   query =>  $some_value, # string
   all =>  $some_value, # int
   noorf =>  $some_value, # int
   nogene =>  $some_value, # int
   from =>  $some_value, # int
   to =>  $some_value, # int
   feattype =>  $some_value, # string
   type =>  $some_value, # string
   chromosome =>  $some_value, # string
   left =>  $some_value, # int
   right =>  $some_value, # int
   strand =>  $some_value, # int
   features =>  $some_value, # string
   feat_format =>  $some_value, # string
   mask_coding =>  $some_value, # int
   repeat =>  $some_value, # int
   all_transcripts =>  $some_value, # int
   first_intron =>  $some_value, # int
   non_coding =>  $some_value, # int
 },




=head1 AUTHOR

Generated by SOAP::WSDL

=cut

