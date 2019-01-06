#!/usr/bin/python
#-*-coding: utf-8-*-

##### ---------- imports ---------- #####

#Used to set parser and subparsers
import argparse
from ensembl_parser import ensembl_parser
from ensembl_parser import help_message
from ensembl_functions import retrieve_features
from ensembl_functions import retrieve_species
import sys


##### ---------- Main Program ---------- #####
if __name__ == '__main__':

    args = ensembl_parser()

    if args.sp is None:
        help_message()
    elif args.sp == 'retrieve_features':
        retrieve_features(args.org, args.f, args.depth, args.verbosity, args.outputdir, args.dump, args.null)
    elif args.sp == 'retrieve_species':
        retrieve_species(args.output, args.verbosity, args.file, args.dump, args.null, args.division)
