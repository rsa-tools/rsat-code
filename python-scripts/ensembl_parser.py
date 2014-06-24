# !/usr/bin/python
#-*-coding: utf-8-*-

import argparse
import sys
import subprocess as subp
import platform
import os


class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        # this is the RawTextHelpFormatter._split_lines
        if text.startswith('R|'):
            return text[2:].splitlines()
        return argparse.HelpFormatter._split_lines(self, text, width)


def help_message():
    if platform.system() == 'Windows':
        subp.call('cls', shell=True)
    elif platform.system() == 'Linux':
        subp.call('clear', shell=True)

    print("\t\t\t-------------------------")
    print("\t\t\t     ENSEMBLGENOMES      ")
    print("\t\t\t-------------------------")
    print("\nDescription:\n\tThis program fetches species-related data from ensemblGenomes\n"
          "\tDataBase, such as the list of supported species, and their\n"
          "\tannotations (genes, transcripts, proteins).")

    print("\n\nUsage: ensembl_genomes <command> [<args>]")
    print("\n\nCommand list:")
    print("\tretrieve_species\tCollect informations on all species of EnsemblGenomes")
    print("\tretrieve_features\tCollect genes informations of one or several species")

    print("\nGetting help:")
    print("\t-h\t\tDisplay overall help message")
    print("\t<command> -h\tDisplay help message for the command")

    print("\nretrieve_species args:")
    print("\t-o --output\tChoose the output file name containing the informations about all species")
    print("\t-v --verbosity\tChoose the verbosity level\n\t\t\t\tdefault: None \n\t\t\t\t1: few informations\n\t\t\t\t2: detailed informations\n\t\t\t\t3: complete")
    print("\t-f --file\t-f <name> create a file with only the species name, one per line")

    print("\nretrieve_features args:")
    print("\t-org\t\t-org <name> to collect genes infos related to <name>")
    print("\t-f\t\t-f <file> to collect genes infos related to species into <file>")
    print("\t-d --depth\tChoose the parsing depth level (each include the previous ones) :\n"
          "\t\t\t\t1: genes\n"
          "\t\t\t\t2: transcripts\n"
          "\t\t\t\t3: translations\n"
          "\t\t\t\t4: cross-refs of proteins\n"
          "\t\t\t\t5: linkage types")
    print("\t-v --verbosity\tChoose the verbosity level")
    print("\t--outputdir\tCreate a folder tree where one folder stands for one specie.\n"
          "\t\t\tWorks well when a specie list is given\n")


def ensembl_parser():

    ##### ---------- Parsing Arguments ---------- #####
    parser = argparse.ArgumentParser(description="This program fetches species-related data from ensemblGenomes DataBase, "
                                                 "such as the list of supported species, and their annotations "
                                                 "(genes, transcripts, proteins).")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('--dump',
                        action='store_true',
                        help="debugging tool, prints the json output")

    subparsers = parser.add_subparsers(dest='sp')

    #Sub parser 1 :
    parser_retrieve_species = subparsers.add_parser("retrieve_species",
                                                    help="Use this command to obtain a file with informations"
                                                         " about the whole species of ensemblGenomes\n ",
                                                    formatter_class=SmartFormatter)
    parser_retrieve_species.required = False
    #Args of the subparser 1
    parser_retrieve_species.add_argument('-o', '--output',
                                         type=str,
                                         default='species_info.txt',
                                         help="Choose the output file name containing the informations about all species\n ")

    parser_retrieve_species.add_argument('-v', '--verbosity',
                                         type=int,
                                         default=0,
                                         choices=[1, 2, 3],
                                         help="R|Choose the verbosity level\ndefault: None \n1: few informations\n2: detailed informations\n3: complete\n ")

    parser_retrieve_species.add_argument('-f', '--file',
                                         type=str,
                                         default=None,
                                         help="If enabled, this option will ensure the creation of a separate "
                                              "file containing only species names (one species per line).\n ")

    parser_retrieve_species.add_argument('-n', '--null',
                        		 type=str,
		                         default="<NULL>",
                                         help="Value to display for NULL field values.")

    parser_retrieve_species.add_argument('-d', '--division',
                                         type=str,
		                         default=None,
                                         help="Ensembl division (e.g. EnsemblBacteria, EnsemblFungi, EnsemblPlants, ...).")

    #Sub parser 2 :
    parser_retrieve_features = subparsers.add_parser("retrieve_features",
                                                     help="Use this command to obtain the genes information related to one "
                                                          "or several input species",
                                                     formatter_class=SmartFormatter)
    parser_retrieve_features.required = False
    #Args of the subparser 2
    exclusiv_group = parser_retrieve_features.add_mutually_exclusive_group(required=True)
    exclusiv_group.add_argument('-org',
                                type=str,
                                help="If you want the genes related to only one organism, enter its name here.\n ")

    exclusiv_group.add_argument('-f',
                                type=str,
                                help="Name of a species list file (one line per species).\n ")

    parser_retrieve_features.add_argument('-d', '--depth',
                                          type=int,
                                          dest='depth',
                                          default=5,
                                          choices=[1, 2, 3, 4, 5],
                                          help="\n".join(["R|Choose the parsing depth, i.e. the number of feature levels to "
                                                          "extract (each level includes all the previous ones)",
                                                          "1: genes",
                                                          "2: transcripts",
                                                          "3: translations",
                                                          "4: cross-refs of proteins",
                                                          "5: linkage types.\n "]))

    parser_retrieve_features.add_argument('-v', '--verbosity',
                                          type=int,
                                          default=0,
                                          choices=[1, 2, 3],
                                          help="R|Choose the verbosity level\ndefault: None \n1: few informations\n2: detailed informations\n3: complete\n ")

    parser_retrieve_features.add_argument('--outputdir',
                                          type=str,
                                          default=None,
                                          help="Create a folder tree where one folder stands for one specie.\n")
    parser_retrieve_features.add_argument('-n', '--null',
                        		  type=str,
		                          default="<NULL>",
                                          help="Value to display for null field values")

    # A way to ensure python 2.7 compatibility :
    # Used to make python believe that you enter a task when you didn't
    # Creating a namespace object
    args = argparse.Namespace()
    # Add a fake subparser
    args.sp = None
    # If there is a subparser chosen || if args.sp is not None
    if sys.argv[1:]:
        # Parsing the user's arguments
        args = parser.parse_args()
    return args
