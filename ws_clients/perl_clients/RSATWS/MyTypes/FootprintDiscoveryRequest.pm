package MyTypes::FootprintDiscoveryRequest;
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
my %genes_of :ATTR(:get<genes>);
my %tmp_infile_of :ATTR(:get<tmp_infile>);
my %all_genes_of :ATTR(:get<all_genes>);
my %max_genes_of :ATTR(:get<max_genes>);
my %output_prefix_of :ATTR(:get<output_prefix>);
my %query_of :ATTR(:get<query>);
my %sep_genes_of :ATTR(:get<sep_genes>);
my %organism_of :ATTR(:get<organism>);
my %taxon_of :ATTR(:get<taxon>);
my %index_of :ATTR(:get<index>);
my %lth_of :ATTR(:get<lth>);
my %uth_of :ATTR(:get<uth>);
my %return_of :ATTR(:get<return>);
my %to_matrix_of :ATTR(:get<to_matrix>);
my %bg_model_of :ATTR(:get<bg_model>);
my %no_filter_of :ATTR(:get<no_filter>);
my %infer_operons_of :ATTR(:get<infer_operons>);
my %dist_thr_of :ATTR(:get<dist_thr>);

__PACKAGE__->_factory(
    [ qw(        output
        verbosity
        genes
        tmp_infile
        all_genes
        max_genes
        output_prefix
        query
        sep_genes
        organism
        taxon
        index
        lth
        uth
        return
        to_matrix
        bg_model
        no_filter
        infer_operons
        dist_thr

    ) ],
    {
        'output' => \%output_of,
        'verbosity' => \%verbosity_of,
        'genes' => \%genes_of,
        'tmp_infile' => \%tmp_infile_of,
        'all_genes' => \%all_genes_of,
        'max_genes' => \%max_genes_of,
        'output_prefix' => \%output_prefix_of,
        'query' => \%query_of,
        'sep_genes' => \%sep_genes_of,
        'organism' => \%organism_of,
        'taxon' => \%taxon_of,
        'index' => \%index_of,
        'lth' => \%lth_of,
        'uth' => \%uth_of,
        'return' => \%return_of,
        'to_matrix' => \%to_matrix_of,
        'bg_model' => \%bg_model_of,
        'no_filter' => \%no_filter_of,
        'infer_operons' => \%infer_operons_of,
        'dist_thr' => \%dist_thr_of,
    },
    {
        'output' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'verbosity' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'genes' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_infile' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'all_genes' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'max_genes' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'output_prefix' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'query' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'sep_genes' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'organism' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'taxon' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'index' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'lth' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'uth' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'return' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'to_matrix' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'bg_model' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'no_filter' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'infer_operons' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'dist_thr' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
    },
    {

        'output' => 'output',
        'verbosity' => 'verbosity',
        'genes' => 'genes',
        'tmp_infile' => 'tmp_infile',
        'all_genes' => 'all_genes',
        'max_genes' => 'max_genes',
        'output_prefix' => 'output_prefix',
        'query' => 'query',
        'sep_genes' => 'sep_genes',
        'organism' => 'organism',
        'taxon' => 'taxon',
        'index' => 'index',
        'lth' => 'lth',
        'uth' => 'uth',
        'return' => 'return',
        'to_matrix' => 'to_matrix',
        'bg_model' => 'bg_model',
        'no_filter' => 'no_filter',
        'infer_operons' => 'infer_operons',
        'dist_thr' => 'dist_thr',
    }
);

} # end BLOCK







1;


=pod

=head1 NAME

MyTypes::FootprintDiscoveryRequest

=head1 DESCRIPTION

Perl data type class for the XML Schema defined complexType
FootprintDiscoveryRequest from the namespace urn:RSATWS.

Parameters for the operation footprint_discovery.




=head2 PROPERTIES

The following properties may be accessed using get_PROPERTY / set_PROPERTY
methods:

=over

=item * output


=item * verbosity


=item * genes


=item * tmp_infile


=item * all_genes


=item * max_genes


=item * output_prefix


=item * query


=item * sep_genes


=item * organism


=item * taxon


=item * index


=item * lth


=item * uth


=item * return


=item * to_matrix


=item * bg_model


=item * no_filter


=item * infer_operons


=item * dist_thr




=back


=head1 METHODS

=head2 new

Constructor. The following data structure may be passed to new():

 { # MyTypes::FootprintDiscoveryRequest
   output =>  $some_value, # string
   verbosity =>  $some_value, # int
   genes =>  $some_value, # string
   tmp_infile =>  $some_value, # string
   all_genes =>  $some_value, # int
   max_genes =>  $some_value, # int
   output_prefix =>  $some_value, # string
   query =>  $some_value, # string
   sep_genes =>  $some_value, # int
   organism =>  $some_value, # string
   taxon =>  $some_value, # string
   index =>  $some_value, # int
   lth =>  $some_value, # string
   uth =>  $some_value, # string
   return =>  $some_value, # string
   to_matrix =>  $some_value, # int
   bg_model =>  $some_value, # string
   no_filter =>  $some_value, # int
   infer_operons =>  $some_value, # int
   dist_thr =>  $some_value, # int
 },




=head1 AUTHOR

Generated by SOAP::WSDL

=cut

