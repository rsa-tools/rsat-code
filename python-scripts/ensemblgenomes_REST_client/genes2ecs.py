#!/usr/bin/python
# -*-coding: utf-8-*-

import sys, argparse
from ensembl_rest_client_lib import *


def extract_match(ids_list, dico):
    """Extract from a GPE file the genes corresponding to the IDs from ids_list.
        If there's a dictionnary, add the gene names before the gene IDs in the output.

    :param id_list: list of gene IDs of interest
    :param dico: dictionnary with IDs for keys and the corresponding names for values
    :return: None
    """
    fh = open_output_handle(args.outfile, verbosity_level=0)

    gpe_file = open(args.gpe_file, "r")

    # Parse the given file
    while True:
        line = gpe_file.readline()
        # Stop the while at the end of the file
        if line == "":
            break
        # Ignore comment line beginning with ';'
        if line == '\n' or line.startswith(';') or line.startswith('";'):
            pass
        else:
            list_temp = line.split('\t')

            if list_temp[0] in ids_list:
                if dico:
                    fh.write(dico[list_temp[0]].strip("\n") + "\t" + line)
                else:
                    fh.write(line)
    fh.close()


if __name__ == '__main__':

    # Specify and docment the command line arguments
    parser = argparse.ArgumentParser(description="Given a set of gene IDs and a gene-protein-EC (GPE) file, " +
                                                 "return the lines of the GPE that match the query genes.",
                                     epilog="Jeanne Cheneby, Justine Long and Jacques van Helden",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-g", "--gpe_file",
                        help="File containing gene-protein-EC (GPE) relationships for a species of interest.",
                        action="store",
                        required=True)

    parser.add_argument("-q", "--query_file",
                        help="\n".join(("File indicating the query genes.",
                                        "Each gene must be specified by its ID (in future version,",
                                        "gene names should be supported).",
                                        "By default the queries are expected to be found in the first column,",
                                        "but the query column can be changed with the option -c.")),
                        action="store",
                        required=True)

    parser.add_argument("-c", "--column_nb",
                        help="\n".join(
                            ("Index of the column in which the species' id are found in the input file(query_file).",
                             "By default the information extract are from the first column.")),
                        type=int,
                        action="store",
                        default=1)

    parser.add_argument("-n", "--name_file",
                        help="\n".join(("File with names and ids of genes.",
                                        "The names are expected to be found in the first column,",
                                        "and ids in the second")),
                        action="store", )

    parser.add_argument("-o", "--outfile", help="Put the result in an outfile named by the user.", action="store")

    args = parser.parse_args()

    # Get the info of interest from the query list
    query_list = parsing_file(args.query_file, args.column_nb)

    # Initialize a dictionnary
    dic_id = {}

    # Check if it's a name file supplied
    if args.name_file:
        names_list = []
        ids_list = []

        # Create a list with each name cleaned from the query_list 
        for name in query_list:
            names_list.append(name.strip("\n"))

        # Get all the names and IDs from the name_file
        all_ids = parsing_file(args.name_file, 1)
        all_names = parsing_file(args.name_file, 2)

        for i in range(len(all_ids)):
            if all_names[i].strip("\n") in names_list:
                dic_id[all_ids[i]] = all_names[i]
                ids_list.append(all_ids[i].strip("\n"))


    else:
        ids_list = []
        for id in query_list:
            ids_list.append(id.strip("\n"))

    extract_match(ids_list, dic_id)
