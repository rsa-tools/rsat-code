package MyInterfaces::RSATWebServices::RSATWSPortType;
use strict;
use warnings;
use Class::Std::Fast::Storable;
use Scalar::Util qw(blessed);
use base qw(SOAP::WSDL::Client::Base);

# only load if it hasn't been loaded before
require MyTypemaps::RSATWebServices
    if not MyTypemaps::RSATWebServices->can('get_class');

sub START {
    $_[0]->set_proxy('http://rsat.ulb.ac.be/rsat/web_services/RSATWS.cgi') if not $_[2]->{proxy};
    $_[0]->set_class_resolver('MyTypemaps::RSATWebServices')
        if not $_[2]->{class_resolver};

    $_[0]->set_prefix($_[2]->{use_prefix}) if exists $_[2]->{use_prefix};
}

sub retrieve_seq {
    my ($self, $body, $header) = @_;
    die "retrieve_seq must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'retrieve_seq',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::retrieve_seq )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub retrieve_seq_multigenome {
    my ($self, $body, $header) = @_;
    die "retrieve_seq_multigenome must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'retrieve_seq_multigenome',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::retrieve_seq_multigenome )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub retrieve_ensembl_seq {
    my ($self, $body, $header) = @_;
    die "retrieve_ensembl_seq must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'retrieve_ensembl_seq',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::retrieve_ensembl_seq )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub purge_seq {
    my ($self, $body, $header) = @_;
    die "purge_seq must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'purge_seq',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::purge_seq )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub oligo_analysis {
    my ($self, $body, $header) = @_;
    die "oligo_analysis must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'oligo_analysis',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::oligo_analysis )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub dyad_analysis {
    my ($self, $body, $header) = @_;
    die "dyad_analysis must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'dyad_analysis',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::dyad_analysis )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub pattern_assembly {
    my ($self, $body, $header) = @_;
    die "pattern_assembly must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'pattern_assembly',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::pattern_assembly )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub dna_pattern {
    my ($self, $body, $header) = @_;
    die "dna_pattern must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'dna_pattern',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::dna_pattern )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub convert_features {
    my ($self, $body, $header) = @_;
    die "convert_features must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'convert_features',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::convert_features )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub feature_map {
    my ($self, $body, $header) = @_;
    die "feature_map must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'feature_map',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::feature_map )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub footprint_discovery {
    my ($self, $body, $header) = @_;
    die "footprint_discovery must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'footprint_discovery',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::footprint_discovery )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub get_orthologs {
    my ($self, $body, $header) = @_;
    die "get_orthologs must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'get_orthologs',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::get_orthologs )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub infer_operon {
    my ($self, $body, $header) = @_;
    die "infer_operon must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'infer_operon',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::infer_operon )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub gene_info {
    my ($self, $body, $header) = @_;
    die "gene_info must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'gene_info',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::gene_info )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub supported_organisms {
    my ($self, $body, $header) = @_;
    die "supported_organisms must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'supported_organisms',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::supported_organisms )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub text_to_html {
    my ($self, $body, $header) = @_;
    die "text_to_html must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'text_to_html',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::text_to_html )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub parse_psi_xml {
    my ($self, $body, $header) = @_;
    die "parse_psi_xml must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'parse_psi_xml',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::parse_psi_xml )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub roc_stats {
    my ($self, $body, $header) = @_;
    die "roc_stats must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'roc_stats',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::roc_stats )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub classfreq {
    my ($self, $body, $header) = @_;
    die "classfreq must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'classfreq',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::classfreq )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub xygraph {
    my ($self, $body, $header) = @_;
    die "xygraph must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'xygraph',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::xygraph )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub convert_seq {
    my ($self, $body, $header) = @_;
    die "convert_seq must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'convert_seq',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::convert_seq )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub compare_classes {
    my ($self, $body, $header) = @_;
    die "compare_classes must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'compare_classes',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::compare_classes )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub convert_classes {
    my ($self, $body, $header) = @_;
    die "convert_classes must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'convert_classes',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::convert_classes )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub contingency_stats {
    my ($self, $body, $header) = @_;
    die "contingency_stats must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'contingency_stats',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::contingency_stats )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub contingency_table {
    my ($self, $body, $header) = @_;
    die "contingency_table must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'contingency_table',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::contingency_table )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub matrix_scan {
    my ($self, $body, $header) = @_;
    die "matrix_scan must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'matrix_scan',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::matrix_scan )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub convert_matrix {
    my ($self, $body, $header) = @_;
    die "convert_matrix must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'convert_matrix',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::convert_matrix )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub matrix_distrib {
    my ($self, $body, $header) = @_;
    die "matrix_distrib must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'matrix_distrib',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::matrix_distrib )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub random_seq {
    my ($self, $body, $header) = @_;
    die "random_seq must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'random_seq',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::random_seq )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub convert_graph {
    my ($self, $body, $header) = @_;
    die "convert_graph must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'convert_graph',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::convert_graph )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub alter_graph {
    my ($self, $body, $header) = @_;
    die "alter_graph must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'alter_graph',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::alter_graph )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub graph_cliques {
    my ($self, $body, $header) = @_;
    die "graph_cliques must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'graph_cliques',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::graph_cliques )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub display_graph {
    my ($self, $body, $header) = @_;
    die "display_graph must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'display_graph',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::display_graph )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub draw_heatmap {
    my ($self, $body, $header) = @_;
    die "draw_heatmap must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'draw_heatmap',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::draw_heatmap )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub compare_graphs {
    my ($self, $body, $header) = @_;
    die "compare_graphs must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'compare_graphs',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::compare_graphs )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub graph_neighbours {
    my ($self, $body, $header) = @_;
    die "graph_neighbours must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'graph_neighbours',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::graph_neighbours )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub mcl {
    my ($self, $body, $header) = @_;
    die "mcl must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'mcl',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::mcl )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub rnsc {
    my ($self, $body, $header) = @_;
    die "rnsc must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'rnsc',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::rnsc )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub graph_node_degree {
    my ($self, $body, $header) = @_;
    die "graph_node_degree must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'graph_node_degree',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::graph_node_degree )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub graph_topology {
    my ($self, $body, $header) = @_;
    die "graph_topology must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'graph_topology',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::graph_topology )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub graph_cluster_membership {
    my ($self, $body, $header) = @_;
    die "graph_cluster_membership must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'graph_cluster_membership',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::graph_cluster_membership )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub graph_get_clusters {
    my ($self, $body, $header) = @_;
    die "graph_get_clusters must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'graph_get_clusters',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::graph_get_clusters )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub random_graph {
    my ($self, $body, $header) = @_;
    die "random_graph must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'random_graph',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::random_graph )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub monitor {
    my ($self, $body, $header) = @_;
    die "monitor must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'monitor',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::monitor )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}


sub get_result {
    my ($self, $body, $header) = @_;
    die "get_result must be called as object method (\$self is <$self>)" if not blessed($self);
    return $self->SUPER::call({
        operation => 'get_result',
        soap_action => '',
        style => 'document',
        body => {
            

           'use'            => 'literal',
            namespace       => 'http://schemas.xmlsoap.org/wsdl/soap/',
            encodingStyle   => '',
            parts           =>  [qw( MyElements::get_result )],
        },
        header => {
            
        },
        headerfault => {
            
        }
    }, $body, $header);
}




1;



__END__

=pod

=head1 NAME

MyInterfaces::RSATWebServices::RSATWSPortType - SOAP Interface for the RSATWebServices Web Service

=head1 SYNOPSIS

 use MyInterfaces::RSATWebServices::RSATWSPortType;
 my $interface = MyInterfaces::RSATWebServices::RSATWSPortType->new();

 my $response;
 $response = $interface->retrieve_seq();
 $response = $interface->retrieve_seq_multigenome();
 $response = $interface->retrieve_ensembl_seq();
 $response = $interface->purge_seq();
 $response = $interface->oligo_analysis();
 $response = $interface->dyad_analysis();
 $response = $interface->pattern_assembly();
 $response = $interface->dna_pattern();
 $response = $interface->convert_features();
 $response = $interface->feature_map();
 $response = $interface->footprint_discovery();
 $response = $interface->get_orthologs();
 $response = $interface->infer_operon();
 $response = $interface->gene_info();
 $response = $interface->supported_organisms();
 $response = $interface->text_to_html();
 $response = $interface->parse_psi_xml();
 $response = $interface->roc_stats();
 $response = $interface->classfreq();
 $response = $interface->xygraph();
 $response = $interface->convert_seq();
 $response = $interface->compare_classes();
 $response = $interface->convert_classes();
 $response = $interface->contingency_stats();
 $response = $interface->contingency_table();
 $response = $interface->matrix_scan();
 $response = $interface->convert_matrix();
 $response = $interface->matrix_distrib();
 $response = $interface->random_seq();
 $response = $interface->convert_graph();
 $response = $interface->alter_graph();
 $response = $interface->graph_cliques();
 $response = $interface->display_graph();
 $response = $interface->draw_heatmap();
 $response = $interface->compare_graphs();
 $response = $interface->graph_neighbours();
 $response = $interface->mcl();
 $response = $interface->rnsc();
 $response = $interface->graph_node_degree();
 $response = $interface->graph_topology();
 $response = $interface->graph_cluster_membership();
 $response = $interface->graph_get_clusters();
 $response = $interface->random_graph();
 $response = $interface->monitor();
 $response = $interface->get_result();



=head1 DESCRIPTION

SOAP Interface for the RSATWebServices web service
located at http://rsat.ulb.ac.be/rsat/web_services/RSATWS.cgi.

=head1 SERVICE RSATWebServices

Web services for the Regulatory Sequence Analysis Tools (RSAT). Tools developed by Jacques van Helden (jvanheld@bigre.ulb.ac.be), SOAP/WSDL interface developed by Olivier Sand (oly@bigre.ulb.ac.be).

=head2 Port RSATWSPortType



=head1 METHODS

=head2 General methods

=head3 new

Constructor.

All arguments are forwarded to L<SOAP::WSDL::Client|SOAP::WSDL::Client>.

=head2 SOAP Service methods

Method synopsis is displayed with hash refs as parameters.

The commented class names in the method's parameters denote that objects
of the corresponding class can be passed instead of the marked hash ref.

You may pass any combination of objects, hash and list refs to these
methods, as long as you meet the structure.

List items (i.e. multiple occurences) are not displayed in the synopsis.
You may generally pass a list ref of hash refs (or objects) instead of a hash
ref - this may result in invalid XML if used improperly, though. Note that
SOAP::WSDL always expects list references at maximum depth position.

XML attributes are not displayed in this synopsis and cannot be set using
hash refs. See the respective class' documentation for additional information.



=head3 retrieve_seq

Returns upstream, downstream or coding DNA sequences for list of query genes.

 $interface->retrieve_seq( {
    request =>  { # MyTypes::RetrieveSequenceRequest
      output =>  $some_value, # string
      organism =>  $some_value, # string
      query =>  $some_value, # string
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
    },
  },,
 );

=head3 retrieve_seq_multigenome

Returns upstream, downstream or coding DNA sequencesfor list of query genes and organisms.

 $interface->retrieve_seq_multigenome( {
    request =>  { # MyTypes::RetrieveSequenceMultigenomeRequest
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
  },,
 );

=head3 retrieve_ensembl_seq

Returns upstream, downstream or coding DNA sequences for list of query genes (in EnsEMBL database).

 $interface->retrieve_ensembl_seq( {
    request =>  { # MyTypes::RetrieveEnsemblSequenceRequest
      output =>  $some_value, # string
      organism =>  $some_value, # string
      ensembl_host =>  $some_value, # string
      db_name =>  $some_value, # string
      query =>  $some_value, # string
      tmp_infile =>  $some_value, # string
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
      unique_sequences =>  $some_value, # int
      first_intron =>  $some_value, # int
      non_coding =>  $some_value, # int
      utr =>  $some_value, # string
      line_width =>  $some_value, # int
      ortho =>  $some_value, # int
      taxon =>  $some_value, # string
      homology_type =>  $some_value, # string
      header_organism =>  $some_value, # string
    },
  },,
 );

=head3 purge_seq

Mask repeated fragments of an input sequence.

 $interface->purge_seq( {
    request =>  { # MyTypes::PurgeSequenceRequest
      output =>  $some_value, # string
      sequence =>  $some_value, # string
      tmp_infile =>  $some_value, # string
      format =>  $some_value, # string
      match_length =>  $some_value, # int
      mismatch =>  $some_value, # int
      str =>  $some_value, # int
      delete =>  $some_value, # int
      mask_short =>  $some_value, # int
    },
  },,
 );

=head3 oligo_analysis

Analysis of the statistical significance of all the oligomers of a given size in a sequence. Commonly used to detect over-represented oligonucleotides in a set of promoter sequences.

 $interface->oligo_analysis( {
    request =>  { # MyTypes::OligoAnalysisRequest
      output =>  $some_value, # string
      verbosity =>  $some_value, # int
      sequence =>  $some_value, # string
      tmp_infile =>  $some_value, # string
      format =>  $some_value, # string
      length =>  $some_value, # int
      organism =>  $some_value, # string
      background =>  $some_value, # string
      stats =>  $some_value, # string
      noov =>  $some_value, # int
      str =>  $some_value, # int
      sort =>  $some_value, # int
      lth =>  $some_value, # string
      uth =>  $some_value, # string
      pseudo =>  $some_value, # string
    },
  },,
 );

=head3 dyad_analysis

Analysis of the statistical significance of all the spaced dyads of a given size in a sequence. Commonly used to detect over-represented spaced dyads in a set of promoter sequences.

 $interface->dyad_analysis( {
    request =>  { # MyTypes::DyadAnalysisRequest
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
  },,
 );

=head3 pattern_assembly

Assemble a set of oligonucleotides or dyads into groups of overlapping patterns (assemblies).

 $interface->pattern_assembly( {
    request =>  { # MyTypes::PatternAssemblyRequest
      output =>  $some_value, # string
      input =>  $some_value, # string
      tmp_infile =>  $some_value, # string
      verbosity =>  $some_value, # int
      score_col =>  $some_value, # int
      str =>  $some_value, # int
      maxfl =>  $some_value, # int
      subst =>  $some_value, # int
      max_asmb_nb =>  $some_value, # int
      max_asmb_size =>  $some_value, # int
      maxpat =>  $some_value, # int
      toppat =>  $some_value, # int
    },
  },,
 );

=head3 dna_pattern

Searches all occurrences of a pattern within DNA sequences.

 $interface->dna_pattern( {
    request =>  { # MyTypes::DnaPatternRequest
      output =>  $some_value, # string
      sequence =>  $some_value, # string
      tmp_infile =>  $some_value, # string
      format =>  $some_value, # string
      subst =>  $some_value, # int
      pattern =>  $some_value, # string
      pattern_file =>  $some_value, # string
      tmp_pattern_file =>  $some_value, # string
      id =>  $some_value, # string
      origin =>  $some_value, # string
      noov =>  $some_value, # int
      score =>  $some_value, # int
      str =>  $some_value, # int
      sort =>  $some_value, # int
      th =>  $some_value, # int
      return =>  $some_value, # string
    },
  },,
 );

=head3 convert_features

Interconversions between various formats of feature description.

 $interface->convert_features( {
    request =>  { # MyTypes::ConvertFeaturesRequest
      output =>  $some_value, # string
      input =>  $some_value, # string
      tmp_infile =>  $some_value, # string
      from =>  $some_value, # string
      to =>  $some_value, # string
    },
  },,
 );

=head3 feature_map

Draws a graphical map of features (e.g. results of pattern matching) in a set of sequences.

 $interface->feature_map( {
    request =>  { # MyTypes::FeatureMapRequest
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
  },,
 );

=head3 footprint_discovery

Detect phylogenetic footprints by applying dyad-analysis in promoters of a set of orthologous genes.

 $interface->footprint_discovery( {
    request =>  { # MyTypes::FootprintDiscoveryRequest
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
  },,
 );

=head3 get_orthologs

Get orthologuous genes.

 $interface->get_orthologs( {
    request =>  { # MyTypes::GetOrthologsRequest
      output =>  $some_value, # string
      organism =>  $some_value, # string
      taxon =>  $some_value, # string
      query =>  $some_value, # string
      all =>  $some_value, # int
      nogrep =>  $some_value, # int
      return =>  $some_value, # string
      lth =>  $some_value, # string
      uth =>  $some_value, # string
    },
  },,
 );

=head3 infer_operon

Infer operon.

 $interface->infer_operon( {
    request =>  { # MyTypes::InferOperonRequest
      output =>  $some_value, # string
      organism =>  $some_value, # string
      query =>  $some_value, # string
      tmp_infile =>  $some_value, # string
      all =>  $some_value, # int
      distance =>  $some_value, # int
      return =>  $some_value, # string
    },
  },,
 );

=head3 gene_info

Get information about genes.

 $interface->gene_info( {
    request =>  { # MyTypes::GeneInfoRequest
      output =>  $some_value, # string
      organism =>  $some_value, # string
      query =>  $some_value, # string
      full =>  $some_value, # int
      noquery =>  $some_value, # int
      descr =>  $some_value, # int
      feattype =>  $some_value, # string
    },
  },,
 );

=head3 supported_organisms

List RSAT suppported organisms.

 $interface->supported_organisms( {
    request =>  { # MyTypes::SupportedOrganismsRequest
      output =>  $some_value, # string
      format =>  $some_value, # string
      taxon =>  $some_value, # string
    },
  },,
 );

=head3 text_to_html

Converts a tab-delimited file into a HTML table

 $interface->text_to_html( {
    request =>  { # MyTypes::TextToHtmlRequest
      output =>  $some_value, # string
      inputfile =>  $some_value, # string
      chunk =>  $some_value, # int
      no_sort =>  $some_value, # int
      font =>  $some_value, # string
    },
  },,
 );

=head3 parse_psi_xml

Converts a psi xml file in a tab delimited file

 $interface->parse_psi_xml( {
    request =>  { # MyTypes::parsepsixmlRequest
      output =>  $some_value, # string
      inputfile =>  $some_value, # string
      channels =>  $some_value, # string
      interactor_type =>  $some_value, # string
      uth =>  $some_value, # float
      lth =>  $some_value, # float
    },
  },,
 );

=head3 roc_stats

Computes, from a set of scored results associated with validation labels, the derived statistics (Sn, PPV, FPR), which can be further used to draw a ROC curve.

 $interface->roc_stats( {
    request =>  { # MyTypes::RocStatsRequest
      output =>  $some_value, # string
      inputfile =>  $some_value, # string
      scol =>  $some_value, # int
      lcol =>  $some_value, # int
      status =>  $some_value, # string
      total =>  $some_value, # int
    },
  },,
 );

=head3 classfreq

This script takes a group of numbers (real or integers) and outputs their distribution among classes.

 $interface->classfreq( {
    request =>  { # MyTypes::ClassFreqRequest
      output =>  $some_value, # string
      inputFile =>  $some_value, # string
      classinterval =>  $some_value, # string
      col =>  $some_value, # string
      min =>  $some_value, # string
      max =>  $some_value, # string
      from =>  $some_value, # string
      to =>  $some_value, # string
    },
  },,
 );

=head3 xygraph

Plot a graph and export it.

 $interface->xygraph( {
    request =>  { # MyTypes::XYGraphRequest
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
  },,
 );

=head3 convert_seq

Converts a sequence between two formats (e.g. fasta -> raw).

 $interface->convert_seq( {
    request =>  { # MyTypes::ConvertSeqRequest
      output =>  $some_value, # string
      sequence =>  $some_value, # string
      tmp_infile =>  $some_value, # string
      from =>  $some_value, # string
      to =>  $some_value, # string
    },
  },,
 );

=head3 compare_classes

Compare two class files(the query file and the reference file). Each class of the query file is compared to each class of the reference file. The number of common elements is reported, as well as the probability to observe at least this number of common elements by chance alone.

 $interface->compare_classes( {
    request =>  { # MyTypes::CompareClassesRequest
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
  },,
 );

=head3 convert_classes

Interconversions between different formats of cluster files.

 $interface->convert_classes( {
    request =>  { # MyTypes::ConvertClassesRequest
      output =>  $some_value, # string
      informat =>  $some_value, # string
      outformat =>  $some_value, # string
      member_col =>  $some_value, # string
      class_col =>  $some_value, # string
      score_col =>  $some_value, # int
      min_score =>  $some_value, # string
      inputclasses =>  $some_value, # string
      names =>  $some_value, # string
    },
  },,
 );

=head3 contingency_stats

This programs takes as input a contingency table, and calculates various matching statistics between the rows and columns. The description of these statistics can be found in Brohee and van Helden (2006).

 $interface->contingency_stats( {
    request =>  { # MyTypes::ContingencyStatsRequest
      output =>  $some_value, # string
      inputfile =>  $some_value, # string
      decimals =>  $some_value, # int
      return =>  $some_value, # string
      rsizes =>  $some_value, # string
      csizes =>  $some_value, # string
    },
  },,
 );

=head3 contingency_table

 Create a contingency table from a two-column file.

 $interface->contingency_table( {
    request =>  { # MyTypes::ContingencyTableRequest
      output =>  $some_value, # string
      inputfile =>  $some_value, # string
      col1 =>  $some_value, # int
      col2 =>  $some_value, # int
      margin =>  $some_value, # int
      null =>  $some_value, # int
    },
  },,
 );

=head3 matrix_scan

Scan sequences with one or several position-specific scoring matrices (PSSM) to identify instances of the corresponding motifs(putative sites). This program supports a variety of background models (Bernoulli, Markov chains of any order).

 $interface->matrix_scan( {
    request =>  { # MyTypes::MatrixScanRequest
      output =>  $some_value, # string
      sequence =>  $some_value, # string
      tmp_sequence_infile =>  $some_value, # string
      matrix =>  $some_value, # string
      tmp_matrix_infile =>  $some_value, # string
      sequence_format =>  $some_value, # string
      matrix_format =>  $some_value, # string
      n_treatment =>  $some_value, # string
      consensus_name =>  $some_value, # string
      pseudo =>  $some_value, # int
      equi_pseudo =>  $some_value, # int
      top_matrices =>  $some_value, # int
      background_model =>  $some_value, # string
      tmp_background_infile =>  $some_value, # string
      organism =>  $some_value, # string
      background =>  $some_value, # string
      background_input =>  $some_value, # int
      background_window =>  $some_value, # int
      markov =>  $some_value, # int
      background_pseudo =>  $some_value, # float
      return_fields =>  $some_value, # string
      sort_distrib =>  $some_value, # int
      lth =>  $some_value, # string
      uth =>  $some_value, # string
      str =>  $some_value, # int
      verbosity =>  $some_value, # int
      origin =>  $some_value, # string
      decimals =>  $some_value, # int
      crer_ids =>  $some_value, # int
    },
  },,
 );

=head3 convert_matrix

Performs inter-conversions between various formats of position-specific scoring matrices (PSSM). The program also performs a statistical analysis of the original matrix to provide different position-specific scores (weight, frequencies, information contents), general statistics (E-value, total information content), and synthetic descriptions (consensus).

 $interface->convert_matrix( {
    request =>  { # MyTypes::ConvertMatrixRequest
      output =>  $some_value, # string
      matrix =>  $some_value, # string
      background_format =>  $some_value, # string
      background_pseudo =>  $some_value, # float
      from =>  $some_value, # string
      to =>  $some_value, # string
      return =>  $some_value, # string
      sort =>  $some_value, # string
      top =>  $some_value, # int
      pseudo =>  $some_value, # float
      equi_pseudo =>  $some_value, # int
      base =>  $some_value, # string
      decimals =>  $some_value, # int
      perm =>  $some_value, # int
      max_profile =>  $some_value, # int
      rc =>  $some_value, # int
    },
  },,
 );

=head3 matrix_distrib

Returns the theoretical distribution of matrix weight within the defined background model.

 $interface->matrix_distrib( {
    request =>  { # MyTypes::MatrixDistribRequest
      output =>  $some_value, # string
      matrix_file =>  $some_value, # string
      tmp_matrix_file =>  $some_value, # string
      matrix_format =>  $some_value, # string
      matrix_pseudo =>  $some_value, # int
      background =>  $some_value, # string
      background_pseudo =>  $some_value, # float
      decimals =>  $some_value, # int
      background_format =>  $some_value, # string
    },
  },,
 );

=head3 random_seq

Generates random sequences.

 $interface->random_seq( {
    request =>  { # MyTypes::RandomSequenceRequest
      output =>  $some_value, # string
      sequence_length =>  $some_value, # int
      repetition =>  $some_value, # int
      format =>  $some_value, # string
      line_width =>  $some_value, # int
      type =>  $some_value, # string
      seed =>  $some_value, # int
      alphabet =>  $some_value, # string
      expfreq =>  $some_value, # string
      tmp_expfreq_file =>  $some_value, # string
      bg_model =>  $some_value, # string
      organism =>  $some_value, # string
      oligo_length =>  $some_value, # int
      length_file =>  $some_value, # string
    },
  },,
 );

=head3 convert_graph

Convert graphs between different formats

 $interface->convert_graph( {
    request =>  { # MyTypes::ConvertGraphRequest
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
  },,
 );

=head3 alter_graph

Alter a graph either by adding or removing edges or nodes

 $interface->alter_graph( {
    request =>  { # MyTypes::AlterGraphRequest
      output =>  $some_value, # string
      informat =>  $some_value, # string
      outformat =>  $some_value, # string
      inputgraph =>  $some_value, # string
      wcol =>  $some_value, # int
      scol =>  $some_value, # int
      tcol =>  $some_value, # int
      directed =>  $some_value, # int
      duplicate =>  $some_value, # int
      self =>  $some_value, # int
      target =>  $some_value, # string
      add_nodes =>  $some_value, # string
      rm_nodes =>  $some_value, # string
      add_edges =>  $some_value, # string
      rm_edges =>  $some_value, # string
    },
  },,
 );

=head3 graph_cliques

Find all cliques in a graph

 $interface->graph_cliques( {
    request =>  { # MyTypes::GraphCliquesRequest
      output =>  $some_value, # string
      informat =>  $some_value, # string
      inputgraph =>  $some_value, # string
      scol =>  $some_value, # int
      tcol =>  $some_value, # int
      min_size =>  $some_value, # int
    },
  },,
 );

=head3 display_graph

Produces the figure of a graph

 $interface->display_graph( {
    request =>  { # MyTypes::DisplayGraphRequest
      output =>  $some_value, # string
      informat =>  $some_value, # string
      outformat =>  $some_value, # string
      ewidth =>  $some_value, # int
      inputgraph =>  $some_value, # string
      wcol =>  $some_value, # int
      scol =>  $some_value, # int
      tcol =>  $some_value, # int
      eccol =>  $some_value, # int
      sccol =>  $some_value, # int
      tccol =>  $some_value, # int
      layout =>  $some_value, # int
    },
  },,
 );

=head3 draw_heatmap

Produces the figure of a heatmap

 $interface->draw_heatmap( {
    request =>  { # MyTypes::DrawHeatmapRequest
      output =>  $some_value, # string
      outformat =>  $some_value, # string
      html =>  $some_value, # int
      inputfile =>  $some_value, # string
      row_names =>  $some_value, # int
      no_text =>  $some_value, # int
      col_width =>  $some_value, # int
      row_height =>  $some_value, # int
      min =>  $some_value, # int
      max =>  $some_value, # int
      gradient =>  $some_value, # string
    },
  },,
 );

=head3 compare_graphs

Computes the union / difference or intersection of two graphs

 $interface->compare_graphs( {
    request =>  { # MyTypes::CompareGraphsRequest
      output =>  $some_value, # string
      Qinformat =>  $some_value, # string
      Rinformat =>  $some_value, # string
      outformat =>  $some_value, # string
      outweight =>  $some_value, # string
      Rinputgraph =>  $some_value, # string
      Qinputgraph =>  $some_value, # string
      Qwcol =>  $some_value, # int
      Qscol =>  $some_value, # int
      Qtcol =>  $some_value, # int
      Rwcol =>  $some_value, # int
      Rscol =>  $some_value, # int
      Rtcol =>  $some_value, # int
      return =>  $some_value, # string
      directed =>  $some_value, # int
      self =>  $some_value, # int
    },
  },,
 );

=head3 graph_neighbours

Find the neihbours up to a certain distance of a collection of seed nodes

 $interface->graph_neighbours( {
    request =>  { # MyTypes::GraphNeighboursRequest
      output =>  $some_value, # string
      informat =>  $some_value, # string
      direction =>  $some_value, # string
      all =>  $some_value, # int
      stats =>  $some_value, # int
      self =>  $some_value, # int
      inputgraph =>  $some_value, # string
      seedfile =>  $some_value, # string
      wcol =>  $some_value, # int
      scol =>  $some_value, # int
      tcol =>  $some_value, # int
      steps =>  $some_value, # int
    },
  },,
 );

=head3 mcl

Clustering via Stijn van Dongen MCL algorithm

 $interface->mcl( {
    request =>  { # MyTypes::MCLRequest
      output =>  $some_value, # string
      inputgraph =>  $some_value, # string
      inflation =>  $some_value, # float
    },
  },,
 );

=head3 rnsc

Clustering via Andrew King RNSC algorithm

 $interface->rnsc( {
    request =>  { # MyTypes::RNSCRequest
      output =>  $some_value, # string
      inputgraph =>  $some_value, # string
      max_clust =>  $some_value, # int
      tabulength =>  $some_value, # int
      tabulist =>  $some_value, # int
      naive_stop =>  $some_value, # int
      scale_stop =>  $some_value, # int
      exp_nb =>  $some_value, # int
      div_freq =>  $some_value, # int
      shf_div_len =>  $some_value, # int
    },
  },,
 );

=head3 graph_node_degree

Calculates the in / out / global degree for a selection of seed nodes

 $interface->graph_node_degree( {
    request =>  { # MyTypes::GraphNodeDegreeRequest
      output =>  $some_value, # string
      informat =>  $some_value, # string
      all =>  $some_value, # int
      inputgraph =>  $some_value, # string
      nodefile =>  $some_value, # string
      wcol =>  $some_value, # int
      scol =>  $some_value, # int
      tcol =>  $some_value, # int
    },
  },,
 );

=head3 graph_topology

Calculate the node degree, the closeness and the betweenness of each node and specifies if this node is a seed or a target node.

 $interface->graph_topology( {
    request =>  { # MyTypes::GraphTopologyRequest
      output =>  $some_value, # string
      informat =>  $some_value, # string
      all =>  $some_value, # int
      return =>  $some_value, # string
      inputgraph =>  $some_value, # string
      nodefile =>  $some_value, # string
      directed =>  $some_value, # int
      wcol =>  $some_value, # int
      scol =>  $some_value, # int
      tcol =>  $some_value, # int
    },
  },,
 );

=head3 graph_cluster_membership

Map a clustering result onto a graph, and compute the membership degree between each node and each cluster, on the basis of egdes linking this node to the cluster.

 $interface->graph_cluster_membership( {
    request =>  { # MyTypes::GraphClusterMembershipRequest
      output =>  $some_value, # string
      informat =>  $some_value, # string
      inputgraph =>  $some_value, # string
      clusters =>  $some_value, # string
      stat =>  $some_value, # string
      decimals =>  $some_value, # int
      wcol =>  $some_value, # int
      scol =>  $some_value, # int
      tcol =>  $some_value, # int
    },
  },,
 );

=head3 graph_get_clusters

Compares a graph with a classification/clustering file.

 $interface->graph_get_clusters( {
    request =>  { # MyTypes::GraphGetClustersRequest
      output =>  $some_value, # string
      informat =>  $some_value, # string
      return =>  $some_value, # string
      outformat =>  $some_value, # string
      inputgraph =>  $some_value, # string
      clusters =>  $some_value, # string
      wcol =>  $some_value, # int
      scol =>  $some_value, # int
      tcol =>  $some_value, # int
      distinct =>  $some_value, # int
      induced =>  $some_value, # int
    },
  },,
 );

=head3 random_graph

Generate random graphs either from scratch of from an existing graph using different randomization models

 $interface->random_graph( {
    request =>  { # MyTypes::RandomGraphRequest
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
  },,
 );

=head3 monitor

Monitoring the status of a job

 $interface->monitor( {
    request =>  { # MyTypes::MonitorRequest
      ticket =>  $some_value, # string
    },
  },,
 );

=head3 get_result

Get result of a job

 $interface->get_result( {
    request =>  { # MyTypes::GetResultRequest
      ticket =>  $some_value, # string
    },
  },,
 );



=head1 AUTHOR

Generated by SOAP::WSDL on Tue Feb 24 09:35:47 2009

=cut
