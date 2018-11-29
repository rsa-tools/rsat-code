#!/usr/bin/python
# -*-coding: utf-8-*-


import sys, argparse

from ensembl_rest_client_lib import *


def get_species_for_taxon(taxon):
    """Collect all the species belonging to a given taxon.

    :param taxon: (string) indicating the taxon name or taxon ID (NCBI taxonomic database)
    :return: (list) list of integers (taxon_id)
    """

    # Beware ! The two databases (ensembl and ensemblgenomes)
    # - have different extensions to collect the species belonging to a given taxon,
    # - return results with different structures.

    # Choose the extension depending on the database
    if args.database == "ensembl":
        ext = "/taxonomy/id/"

    else:
        ext = "/info/genomes/taxonomy/"

    # Send request and collect the decoded result
    decoded = treat_one_request(args.database, ext + taxon, args.verbose)

    # Extract the species IDs from the result
    # Must be treated differently depending on the database
    species_ids = []
    if args.database == 'ensemblgenomes':

        for dico in decoded:
            species_ids.append(dico['taxonomy_id'])

    else:

        for dico in decoded['children']:
            species_ids.append(dico['id'])

    return species_ids


def print_species_info(decoded_all, species_ids=None):
    """ Collect species attributes from a json file.
    Write these attributes into the handle fh.

    :param decoded_all: dictionary of a list of species
    :param species_ids: optional list of selected species (info will be printed only for those species)
    :return:None
    """

    # Initialize variables
    groups = str()
    aliases = str()
    n = 0
    to_print = {}

    # Put every info in the database in a list. column separate by "\t", line by "\n" 
    for i in range(len(decoded_all['species'])):

        # Separate each element of groups list with ','
        for x in range(len(decoded_all['species'][i]['groups']) - 1):
            groups += decoded_all['species'][i]['groups'][x] + ', '
        groups += decoded_all['species'][i]['groups'][-1]

        # Separate each element of aliases list with ','
        for x in range(len(decoded_all['species'][i]['aliases']) - 1):
            aliases += decoded_all['species'][i]['aliases'][x] + ', '
        # check if the aliases is an empty list
        if len(decoded_all['species'][i]['aliases']) > 0:
            aliases += decoded_all['species'][i]['aliases'][-1]

        if (not args.taxon) or (decoded_all['species'][i]['taxon_id'] in species_ids):
            n += 1  # Increment species counter

            # Define the line describing the current organism
            to_print[decoded_all['species'][i]['name']] = "\t".join((
                str(decoded_all['species'][i]['accession']),
                str(decoded_all['species'][i]['name']),
                str(decoded_all['species'][i]['display_name']),
                str(decoded_all['species'][i]['release']),
                str(decoded_all['species'][i]['division']),
                str(decoded_all['species'][i]['taxon_id']),
                str(decoded_all['species'][i]['assembly']),
                str(decoded_all['species'][i]['common_name']),
                groups,
                aliases.replace('"', "'"))) + '\n'

        groups = str()
        aliases = str()

    for species in sorted(to_print.keys()):
        try:
            fh.write(to_print[species])
        except UnicodeEncodeError:
            fh.write(repr(to_print[species]))

    verbosity('; Numbers of species found: ' + str(n), 1, args.verbose)


if __name__ == '__main__':

    # Specify and docment the command line arguments
    parser = argparse.ArgumentParser(
        description="\n".join(("This program gets information about species  available"
                               "at ensemblgenomes (http://ensemblgenomes.org/)",
                               "or ensembl (http://ensembl.org/).",
                               "By default all species are returned, " +
                               "but the query can be restricted to a user-specified taxon.")),
        epilog="Jeanne Cheneby, Justine Long and Jacques van Helden",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-v", "--verbose",
                        help="\n".join((
                            "Verbosity. Supported: 0, 1, 2, 3. Default = 0.",
                            "Level 1: include the command line in the output.",
                            "Level 2: provide time warns about ongoing tasks.",
                            "Level 3: detailed information about the requests.")),
                        type=int,
                        choices=[0, 1, 2, 3], default=0)

    parser.add_argument("-d", "--database",
                        help="\n".join(("supported: ensembl | ensemblgenomes",
                                        "default = ensemblgenomes",
                                        "ensembl: collect species from http://ensembl.org/",
                                        "ensemblgenomes: collect species from  http://ensemblgenomes.org/")),
                        choices=["ensembl", "ensemblgenomes"], default="ensemblgenomes")

    parser.add_argument("-t", "--taxon",
                        help="\n".join(("Taxon of interest (name or taxon-id).",
                                        "Will return only the species belonging to the selected taxon. ",
                                        "If not specified, will return all species from the selected database.",
                                        "For better results, use the taxon-id. ",
                                        "Example: --taxon 83333")), action="store")

    parser.add_argument("-o", "--outfile", help="output file", action="store")

    args = parser.parse_args()

    # Start the job
    start_time = current_time()
    verbosity(current_time() + "\t" + "Starting job", 2, args.verbose)

    # Open output file if requested
    fh = open_output_handle(args.outfile, args.verbose)

    # Print the command as a comment in the beginning of the output file
    report_command(sys.argv, fh, args.verbose)

    # Create header
    create_headers(fh, args.verbose, ["accession", "name", "display_name", "release", "division",
                                      "taxon_ID", "assembly", "common_name", "groups", "aliases"])

    # Initialize result structures
    species_ids = None

    # Collect the list of IDs for the selected taxon (if specified),
    # which will then serve as filter for the output.
    if args.taxon:
        species_ids = get_species_for_taxon(args.taxon)

    # Collect info about all species
    decoded_all = treat_one_request(args.database, "/info/species", args.verbose)

    # Print info for the selected species
    print_species_info(decoded_all, species_ids)

    verbosity(current_time() + "\t" + "Job done", 2, args.verbose)

    processing_time(start_time, current_time(), fh, args.verbose)

    if args.outfile:
        fh.close()
