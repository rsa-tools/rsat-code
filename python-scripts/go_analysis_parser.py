import argparse
import sys
import os
### Start parser ###


def go_parser():
    parser = argparse.ArgumentParser(description="This program will do things with GO & stuff")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    subparsers = parser.add_subparsers(dest='task')

    #Task 1 : download_go
    parser_dl_go = subparsers.add_parser("download_go",
                                         description="Download the GO term definitions from the gene ontoloy web site "
                                              "http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo")
    parser_dl_go.add_argument('-o', '--outputFile',
                              metavar='<file>',
                              type=str,
                              help='Output file.',
                              dest='outFile',
                              default='GO_hierarchy.obo')
    
    #Task 2 : parse_go
    parser_parse_go = subparsers.add_parser("parse_go",
                                            description="Parse the obo-formatted file from GO and exportits content in tab-delimited files.")
    parser_parse_go.add_argument('-f', '--file',
                                 metavar='<file>',
                                 type=str,
                                 help="Input file path for the GO hierarchy. If not specified, the obo file will automatically be download, and removed after parsing.",
                                 default=None)

    parser_parse_go.add_argument('--output_dir',
                                 metavar='<path>',
                                 type=str,
                                 dest='outdir',
                                 help="Path where the created files will be stored.",
                                 default=os.getcwd())
    #Task 3 : get_annotations
    parser_get_annotations = subparsers.add_parser("get_annotations",
                                                   description="Return primary table Gene annotation - GO ID")

    parser_get_annotations.add_argument('-org',
                                        type=str,
                                        metavar='<org_name>',
                                        dest='org_name',
                                        default=None)
    parser_get_annotations.add_argument('--output_dir',
                                        type=str,
                                        metavar='<path>',
                                        dest='outdir',
                                        help='Path where the created files will be stored.',
                                        default=os.getcwd())

    # Task 4 : get_ancestors
    parser_ancestor = subparsers.add_parser("ancestor",
                                            description='Search for go term ancestors. Using a set of '
                                                 'gene to term relations, the command will collect '
                                                 'additional relations by retrieving ancestors of each GO terms.')

    parser_ancestor.add_argument('--inputFile', '-i',
                                 metavar='File',
                                 type=argparse.FileType('r'),
                                 required=True,
                                 help='A tabulated flat file with gene ids as a first column and '
                                      'go term as second column.',
                                 dest='inFile',
                                 default=None)

    parser_ancestor.add_argument('--goFile', '-g',
                                 metavar='File',
                                 type=argparse.FileType('r'),
                                 required=True,
                                 help='A tabulated flat file with gene ids as a first column and '
                                      'go term as second column.',
                                 dest='goFile',
                                 default=None)

    parser_ancestor.add_argument('--outputFile', '-o',
                                 metavar='File',
                                 type=argparse.FileType('w'),
                                 help='Output file.',
                                 dest='outFile',
                                 default=None)
    
    # Task 5
    
    parser_expand = subparsers.add_parser("expand",
                                          description='expand program')
    parser_expand.add_argument('-org', '--organism',
                               type=str,
                               metavar='<org_name>',
                               dest='org_name',
                               default=None)
    parser_expand.add_argument('-d', '--description_file',
                               type=str,
                               metavar='<description_file>',
                               dest='description_file',
                               default=None)
    parser_expand.add_argument('-r', '--relation_file',
                               type=str,
                               metavar='<relation_file>',
                               dest='relation_file',
                               default=None)
    parser_expand.add_argument('-a', '--annotation_file',
                               type=str,
                               metavar='<annotation_file>',
                               dest='annotation_file',
                               default=None)


    parser_enrichment = subparsers.add_parser("enrichment",
                                              description="enrichment help message")

    args = argparse.Namespace()
    args.task = None

    # If there are any arguments..
    if sys.argv[1:]:

        args = parser.parse_args()

    go_error_message(args)

    return args
### End parser ###


### Start command line error gestion ###
def go_error_message(args):
    if args.task is None:
        print("\nUsage : python go_analysis.py <task> [<options>]"
              "\n\nCommand 'task' expects one argument, received zero.\nList of accepted arguments:\n"
              "\tdownload_go\n"
              "\tparse_go\n"
              "\tget_annotations\n"
              "\tancestor\n"
              "\tenrichment")
        sys.exit(-1)

    elif args.task == "get_annotations" and args.org_name is None:
        print("\nUsage : python go_analysis.py get_annotations -org <org_name>"
              "\n\nOption '-org' expects one argument, received zero."
              "\nExemple of argument :"
              "\npython go_analysis.py get_annotations -org mycoplasma_genitalium_g37\n")
        sys.exit(-1)

    elif args.task == "expand" and args.org_name is None:
        if ((args.annotation_file is None) or (args.description_file is None) or (args.relation_file is None)):
            print("\nUsage : python go_analysis.py expand -org <org_name>"
                  "\nYou should either specify an organism (-org), or the three required input files (-d, -r, -a)."
                  "\nExamples of arguments :"
                  "\npython go_analysis.py expand [-org mycoplasma_genitalium_g37 | "
                  "-a annotations_table_campylobacter_jejuni.tab -r GO_relations.tab -d GO_descriptions.tab]\n")
            sys.exit(-1)

### End command line error gestion ###
