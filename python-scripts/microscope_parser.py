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
    print("\t\t\t     MicroScope_get      ")
    print("\t\t\t-------------------------")
    print("\nDescription:\n\tThis program fetches species-related data from MicroScope\n"
          "\tDataBase, such as the list of supported species,\n"
          "\ttheir annotations (genes, transcripts, proteins),\n"
          "\tand metabolic netwoks (reactions, pathways.)")

    print("\n\nUsage: microscope_get <command> [<args>]")
    print("\n\nCommand list:")
    print("\torganisms\tCollect informations on all species of MicroScope")
    print("\tgpr\t\tCollect genes informations for one specie")
    print("\treactionInfo\tCollect a reaction information")
    print("\treactionList\tCollect reaction-gene list for one specie")
    print("\treactions\tCollect all reactions informations (for all species)")   
    
    print("\nGetting help:")
    print("\t-h\t\tDisplay overall help message")
    print("\t<command> -h\tDisplay help message for the command")

    print("\norganisms args:")
    print("\t-o --output\tChoose the output file name containing the informations about all species")
#    print("\t-v --verbosity\tChoose the verbosity level\n\t\t\t\tdefault: None \n\t\t\t\t1: few informations\n\t\t\t\t2: detailed informations\n\t\t\t\t3: complete")
#    print("\t-f --file\t-f <name> create a file with only the species name, one per line")

    print("\ngpr args:")
    print("\t-org\t\t-org <name> to collect gene informations related to <name> (ex: 'Escherichia_coli_K-12_DH10B')")
    print("\t-o --output\tChoose the output file name containing the informations about gene information for one specie")
#    print("\t-f\t\t-f <file> to collect genes infos related to species into <file>")
#    print("\t-d --depth\tChoose the parsing depth level (each include the previous ones) :\n"
#          "\t\t\t\t1: genes\n"
#          "\t\t\t\t2: transcripts\n"
#          "\t\t\t\t3: translations\n"
#          "\t\t\t\t4: cross-refs of proteins\n"
#          "\t\t\t\t5: linkage types")
#    print("\t-v --verbosity\tChoose the verbosity level")
#    print("\t--outputdir\tCreate a folder tree where one folder stands for one specie.\n"
#          "\t\t\tWorks well when a specie list is given\n")
    print("\nreactionInfo args:")
    print("\t-mrId\t\t-mrId <name> to collect reaction information related to <name> (ex: '6-ACETYLGLUCOSE-DEACETYLASE-RXN')")
    print("\t-o --output\tChoose the output file name containing the informations about a reaction")

    print("\nreactionList args:")
    print("\t-org\t\t-org <name> to collect reactions list related to <name> (ex: 'Escherichia_coli_K-12_DH10B')")
    print("\t-o --output\tChoose the output file name containing the informations about all reactions for a specie")
    
    print("\nreactions args:")
    print("\t-o --output\tChoose the output file name containing the reaction informations about all reactions")
#    print("\t-v --verbosity\tChoose the verbosity level\n\t\t\t\tdefault: None \n\t\t\t\t1: few informations\n\t\t\t\t2: detailed informations\n\t\t\t\t3: complete")
#    print("\t-f --file\t-f <name> create a file with only the species name, one per line")

    
def microscope_parser():

    ##### ---------- Parsing Arguments ---------- #####
    parser = argparse.ArgumentParser(description="This program fetches species-related data from MicroScope DataBase, "
                                                 "such as the list of supported species, and their annotations "
                                                 "(genes, transcripts, proteins), metabolic networks.")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    subparsers = parser.add_subparsers(dest='sp')

    #Sub parser 1 :
    parser_organisms = subparsers.add_parser("organisms",
                                                    help="Use this command to obtain a file with informations"
                                                         " about the whole species of MicroScope\n ",
                                                    formatter_class=SmartFormatter)
    parser_organisms.required = False
    #Args of the subparser 1
    parser_organisms.add_argument('-o', '--output',
                                         type=str,
                                         default=None,
                                         help="Choose the output file name containing the informations about all species\n ")

#    parser_organisms.add_argument('-v', '--verbosity',
#                                         type=int,
#                                         default=0,
#                                         choices=[1, 2, 3],
#                                         help="R|Choose the verbosity level\ndefault: None \n1: few informations\n2: detailed informations\n3: complete\n ")

#    parser_organisms.add_argument('-f', '--file',
#                                         type=str,
#                                         default=None,
#                                         help="If enabled, this option will ensure the creation of a separate "
#                                              "file containing only species names (one species per line).\n ")

    #Sub parser 2 :
    parser_gpr = subparsers.add_parser("gpr",
                                                     help="Use this command to obtain the gene informations related to one "
                                                          "input species",
                                                     formatter_class=SmartFormatter)
    parser_gpr.required = False
    #Args of the subparser 2
    exclusiv_group = parser_gpr.add_mutually_exclusive_group(required=True)
    exclusiv_group.add_argument('-org',
                                type=str,
                                help="If you want the genes related to only one organism, enter its name here (ex: 'Escherichia_coli_K-12_DH10B').\n ")
    parser_gpr.add_argument('-o', '--output', 
                                 type=str, 
                                 default=None, 
                                 help="Choose the output file name containing the informations about all gpr\n ")

#    exclusiv_group.add_argument('-f',
#                                type=str,
#                                help="Name of a species list file (one line per species).\n ")

#    parser_retrieve_features.add_argument('-d', '--depth',
#                                          type=int,
#                                          dest='depth',
#                                          default=5,
#                                          choices=[1, 2, 3, 4, 5],
#                                          help="\n".join(["R|Choose the parsing depth, i.e. the number of feature levels to "
#                                                          "extract (each level includes all the previous ones)",
#                                                          "1: genes",
#                                                          "2: transcripts",
#                                                          "3: translations",
#                                                          "4: cross-refs of proteins",
#                                                          "5: linkage types.\n "]))

#    parser_retrieve_features.add_argument('-v', '--verbosity',
#                                          type=int,
#                                          default=0,
#                                          choices=[1, 2, 3],
#                                          help="R|Choose the verbosity level\ndefault: None \n1: few informations\n2: detailed informations\n3: complete\n ")

#    parser_retrieve_features.add_argument('--outputdir',
#                                          type=str,
#                                          default=None,
#                                          help="Create a folder tree where one folder stands for one specie.\n")


    #Sub parser 3 :
    parser_reactionInfo = subparsers.add_parser("reactionInfo",
                                                     help="Use this command to obtain the reaction information related to one "
                                                          "input reaction",
                                                     formatter_class=SmartFormatter)
    parser_reactionInfo.required = False
    #Args of the subparser 3
    exclusiv_group = parser_reactionInfo.add_mutually_exclusive_group(required=True)
    exclusiv_group.add_argument('-mrId',
                                type=str,
                                help="If you want the genes related to only one reaction, enter its name here (ex: '6-ACETYLGLUCOSE-DEACETYLASE-RXN').\n ")
    parser_reactionInfo.add_argument('-o', '--output',
                                         type=str,
                                         default=None,
                                         help="Choose the output file name containing the informations about all species\n ")

    #Sub parser 4 :
    parser_reactionList = subparsers.add_parser("reactionList",
							  help="Use this command to obtain the reaction list related to one input specie",
                                                     formatter_class=SmartFormatter)
    parser_reactionList.required = False
    #Args of the subparser 4
    exclusiv_group = parser_reactionList.add_mutually_exclusive_group(required=True)
    exclusiv_group.add_argument('-org',
                                type=str,
                                help="If you want the genes related to only one reaction, enter its name here (ex: 'Escherichia_coli_K-12_DH10B').\n ")
    parser_reactionList.add_argument('-o', '--output',
                                         type=str,
                                         default=None,
                                         help="Choose the output file name containing the informations about all species\n ")

                                
    #Sub parser 5 :
    parser_reactions = subparsers.add_parser("reactions",
							  help="Use this command to obtain the reaction list",
                                                     formatter_class=SmartFormatter)
    parser_reactions.required = False
    #Args of the subparser 4
   # exclusiv_group = parser_reactions.add_mutually_exclusive_group(required=True)
    parser_reactions.add_argument('-o', '--output',
                                         type=str,
                                         default=None,
                                         help="Choose the output file name containing the reaction informations about all species\n ")




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
