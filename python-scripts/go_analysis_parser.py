import argparse
import sys


### Start parser ###
def go_parser():
    parser = argparse.ArgumentParser(description="This program will do things with GO & stuff")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    subparsers = parser.add_subparsers(dest='task')

    #Task 1 : download_go
    parser_dl_go = subparsers.add_parser("download_go",
                                         help="Use it to DL the GO reference at "
                                              "http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo")
    parser_dl_go.add_argument('--outputFile',
                               '-o',
                               metavar='File',
                               type=argparse.FileType('w'),
                               help='Output file.',
                               dest='outFile',
                               default=None)
    
    #Task 2 : parse_go
    parser_parse_go = subparsers.add_parser("parse_go",
                                            help="Will create tables with infos in it")

    #Task 3 : get_annotations
    parser_get_annotations = subparsers.add_parser("get_annotations",
                                                   help="Return primary table Gene annotation - GO ID")

    parser_get_annotations.add_argument('-org',
                                        type=str,
                                        dest='org_name',
                                        default=None)

    # Task 4 : get_ancestors
    parser_ancestor = subparsers.add_parser("ancestor",
                                          help='Search for go term ancestors. Using a set of '
                                          'gene to term relations, the command will collect '
                                          'additional relations by retrieving ancestors of each GO terms.')

    parser_ancestor.add_argument('--inputFile',
                               '-i',
                               metavar='File',
                               type=argparse.FileType('r'),
                               required=True,
                               help='A tabulated flat file with gene ids as a first column and '
                               'go term as second column.',
                               dest='inFile',
                               default=None)

    parser_ancestor.add_argument('--goFile',
                               '-g',
                               metavar='File',
                               type=argparse.FileType('r'),
                               required=True,
                               help='A tabulated flat file with gene ids as a first column and '
                               'go term as second column.',
                               dest='goFile',
                               default=None)

    parser_ancestor.add_argument('--outputFile',
                               '-o',
                               metavar='File',
                               type=argparse.FileType('w'),
                               help='Output file.',
                               dest='outFile',
                               default=None)
    
    # Task 5
    
    parser_expand = subparsers.add_parser("expand",
                                          help='expand program')

    parser_expand.add_argument('-org',
                               type=str,
                               dest='org_name',
                               default=None)

    parser_enrichment = subparsers.add_parser("enrichment",
                                              help="enrichment help message")

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
        sys.exit()

    elif args.task == "get_annotations" and args.org_name is None:
        print("\nUsage : python go_analysis.py get_annotations -org <org_name>"
              "\n\nOption '-org' expects one argument, received zero."
              "\nExemple of argument :"
              "\npython go_analysis.py get_annotations -org mycoplasma_genitalium_g37\n")
        sys.exit()

    elif args.task == "expand" and args.org_name is None:
        print("\nUsage : python go_analysis.py expand -org <org_name>"
              "\n\nOption '-org' expects one argument, received zero."
              "\nExemple of argument :"
              "\npython go_analysis.py expand -org mycoplasma_genitalium_g37\n")
        sys.exit()
        


### End command line error gestion ###