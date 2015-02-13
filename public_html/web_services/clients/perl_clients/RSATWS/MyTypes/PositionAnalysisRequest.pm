package MyTypes::PositionAnalysisRequest;
use strict;
use warnings;


__PACKAGE__->_set_element_form_qualified(0);

sub get_xmlns { 'urn:RSATWS' };

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
my %seq_type_of :ATTR(:get<seq_type>);
my %last_of :ATTR(:get<last>);
my %mask_of :ATTR(:get<mask>);
my %noov_of :ATTR(:get<noov>);
my %str_of :ATTR(:get<str>);
my %class_int_of :ATTR(:get<class_int>);
my %origin_of :ATTR(:get<origin>);
my %offset_of :ATTR(:get<offset>);
my %group_rc_of :ATTR(:get<group_rc>);
my %sort_of :ATTR(:get<sort>);
my %return_of :ATTR(:get<return>);
my %lth_chi_of :ATTR(:get<lth_chi>);
my %lth_sig_of :ATTR(:get<lth_sig>);
my %lth_occ_of :ATTR(:get<lth_occ>);
my %uth_rank_of :ATTR(:get<uth_rank>);
my %max_graphs_of :ATTR(:get<max_graphs>);
my %pattern_of :ATTR(:get<pattern>);
my %tmp_pattern_infile_of :ATTR(:get<tmp_pattern_infile>);
my %score_column_of :ATTR(:get<score_column>);
my %min_pos_of :ATTR(:get<min_pos>);
my %max_pos_of :ATTR(:get<max_pos>);
my %no_check_of :ATTR(:get<no_check>);
my %no_filter_of :ATTR(:get<no_filter>);
my %image_format_of :ATTR(:get<image_format>);
my %title_of :ATTR(:get<title>);

__PACKAGE__->_factory(
    [ qw(        output
        verbosity
        sequence
        tmp_infile
        format
        length
        seq_type
        last
        mask
        noov
        str
        class_int
        origin
        offset
        group_rc
        sort
        return
        lth_chi
        lth_sig
        lth_occ
        uth_rank
        max_graphs
        pattern
        tmp_pattern_infile
        score_column
        min_pos
        max_pos
        no_check
        no_filter
        image_format
        title

    ) ],
    {
        'output' => \%output_of,
        'verbosity' => \%verbosity_of,
        'sequence' => \%sequence_of,
        'tmp_infile' => \%tmp_infile_of,
        'format' => \%format_of,
        'length' => \%length_of,
        'seq_type' => \%seq_type_of,
        'last' => \%last_of,
        'mask' => \%mask_of,
        'noov' => \%noov_of,
        'str' => \%str_of,
        'class_int' => \%class_int_of,
        'origin' => \%origin_of,
        'offset' => \%offset_of,
        'group_rc' => \%group_rc_of,
        'sort' => \%sort_of,
        'return' => \%return_of,
        'lth_chi' => \%lth_chi_of,
        'lth_sig' => \%lth_sig_of,
        'lth_occ' => \%lth_occ_of,
        'uth_rank' => \%uth_rank_of,
        'max_graphs' => \%max_graphs_of,
        'pattern' => \%pattern_of,
        'tmp_pattern_infile' => \%tmp_pattern_infile_of,
        'score_column' => \%score_column_of,
        'min_pos' => \%min_pos_of,
        'max_pos' => \%max_pos_of,
        'no_check' => \%no_check_of,
        'no_filter' => \%no_filter_of,
        'image_format' => \%image_format_of,
        'title' => \%title_of,
    },
    {
        'output' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'verbosity' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'sequence' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_infile' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'format' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'length' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'seq_type' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'last' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'mask' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'noov' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'str' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'class_int' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'origin' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'offset' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'group_rc' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'sort' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'return' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'lth_chi' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'lth_sig' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'lth_occ' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'uth_rank' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'max_graphs' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'pattern' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_pattern_infile' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'score_column' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'min_pos' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'max_pos' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'no_check' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'no_filter' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'image_format' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'title' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
    },
    {

        'output' => 'output',
        'verbosity' => 'verbosity',
        'sequence' => 'sequence',
        'tmp_infile' => 'tmp_infile',
        'format' => 'format',
        'length' => 'length',
        'seq_type' => 'seq_type',
        'last' => 'last',
        'mask' => 'mask',
        'noov' => 'noov',
        'str' => 'str',
        'class_int' => 'class_int',
        'origin' => 'origin',
        'offset' => 'offset',
        'group_rc' => 'group_rc',
        'sort' => 'sort',
        'return' => 'return',
        'lth_chi' => 'lth_chi',
        'lth_sig' => 'lth_sig',
        'lth_occ' => 'lth_occ',
        'uth_rank' => 'uth_rank',
        'max_graphs' => 'max_graphs',
        'pattern' => 'pattern',
        'tmp_pattern_infile' => 'tmp_pattern_infile',
        'score_column' => 'score_column',
        'min_pos' => 'min_pos',
        'max_pos' => 'max_pos',
        'no_check' => 'no_check',
        'no_filter' => 'no_filter',
        'image_format' => 'image_format',
        'title' => 'title',
    }
);

} # end BLOCK







1;


=pod

=head1 NAME

MyTypes::PositionAnalysisRequest

=head1 DESCRIPTION

Perl data type class for the XML Schema defined complexType
PositionAnalysisRequest from the namespace urn:RSATWS.

Parameters for the operation position_analysis.




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


=item * seq_type


=item * last


=item * mask


=item * noov


=item * str


=item * class_int


=item * origin


=item * offset


=item * group_rc


=item * sort


=item * return


=item * lth_chi


=item * lth_sig


=item * lth_occ


=item * uth_rank


=item * max_graphs


=item * pattern


=item * tmp_pattern_infile


=item * score_column


=item * min_pos


=item * max_pos


=item * no_check


=item * no_filter


=item * image_format


=item * title




=back


=head1 METHODS

=head2 new

Constructor. The following data structure may be passed to new():

 { # MyTypes::PositionAnalysisRequest
   output =>  $some_value, # string
   verbosity =>  $some_value, # int
   sequence =>  $some_value, # string
   tmp_infile =>  $some_value, # string
   format =>  $some_value, # string
   length =>  $some_value, # int
   seq_type =>  $some_value, # string
   last =>  $some_value, # int
   mask =>  $some_value, # string
   noov =>  $some_value, # int
   str =>  $some_value, # int
   class_int =>  $some_value, # int
   origin =>  $some_value, # string
   offset =>  $some_value, # int
   group_rc =>  $some_value, # int
   sort =>  $some_value, # int
   return =>  $some_value, # string
   lth_chi =>  $some_value, # int
   lth_sig =>  $some_value, # int
   lth_occ =>  $some_value, # int
   uth_rank =>  $some_value, # int
   max_graphs =>  $some_value, # int
   pattern =>  $some_value, # string
   tmp_pattern_infile =>  $some_value, # string
   score_column =>  $some_value, # int
   min_pos =>  $some_value, # int
   max_pos =>  $some_value, # int
   no_check =>  $some_value, # int
   no_filter =>  $some_value, # int
   image_format =>  $some_value, # string
   title =>  $some_value, # string
 },




=head1 AUTHOR

Generated by SOAP::WSDL

=cut

