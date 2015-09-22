#!/usr/bin/python
# -*-coding: utf-8-*-

import sys, argparse, pprint, os, ntpath

from ensembl_rest_client_lib import *


def initializing_dir(species):
    """ Choose the species directory and create a dictionary of the output files.

    :param species: (string) indicating the species name
    :return: dictionary of all open files and the start time
    """
    # Specify the species directory
    species_dir = os.path.join(args.outdir, species)

    # Specify the output file names
    output_files = {}
    output_files["GPE"] = os.path.join(species_dir, "gene_protein_ec.tab")

    output_files["Genes (json)"] = dico_entries(args.export_json, species_dir, "gene.json")
    output_files["Gene attributes"] = dico_entries(args.export_genes, species_dir, "gene.tab")
    output_files["Gene names"] = dico_entries(args.export_names, species_dir, "gene_names.tab")

    # Initialize time for the current species
    start_time = current_time()

    return output_files, start_time


def dico_entries(args, species_dir, file_name):
    """Create key in the dictionary output_files.

    :param args: args of the option
    :param species_dir: directory of the species
    :param file_name: name of the file
    :return: value of the key
    """
    if args:
        value = os.path.join(species_dir, file_name)
    else:
        value = None

    return value


def get_GPE(species, output_files, start_time):
    """Collect the gene - EC - protein relationships (GEP) for a given species.
    Write into the handle fh the geneID (from NCBI), proteinID (from SWISSPROT) and ECs (from IntEnz).

    :param species: (string) indicating the species name
    :param output_files: (dictionary) containing name of the output files
    :param start_time: time of the starting job
    :return: None
    """

    # Define extension for the REST URL
    ext = "/lookup/genome/" + species + "?level=protein_feature&xrefs=1/"

    verbosity(current_time() + "\t" + "Starting job", 2, args.verbose)

    result_json = treat_one_request("ensemblgenomes", ext, args.verbose)

    # Open file to store Gene-Protein-EC relationships
    fh = open_output_handle(output_files["GPE"], args.verbose)
    # Print the command as a comment in the beginning of the output file
    report_command(sys.argv, fh, args.verbose)
    # Print headers
    create_headers(fh, args.verbose, ['gene_id', 'protein_id', 'EC'])

    # Open file to store genes attributes
    if args.export_genes:
        fa = open_output_handle(output_files["Gene attributes"], args.verbose)
        create_headers(fa, args.verbose, ['id', 'biotype', 'name', 'contig', 'start', 'end', 'strand', 'description'])

    # Export the raw JSON object in a text file if specified
    if args.export_json:
        jh = open_output_handle(output_files["Genes (json)"], args.verbose)
        pprint.pprint(result_json, stream=jh)
        jh.close()

    if args.export_names:
        nh = open_output_handle(output_files["Gene names"], args.verbose)
        create_headers(nh, args.verbose, ["id", "name"])

    # Get EC and protein ID (from UNIPROT) ONLY if EC exist.
    for j in range(len(result_json)):
        list_temp_ec = []
        indice_uniprot = None
        try:
            # Parse databases per gene
            for i in range(len(result_json[j]['transcripts'][0]['translations'][0]['xrefs'])):
                if 'IntEnz' in result_json[j]['transcripts'][0]['translations'][0]['xrefs'][i]['dbname']:
                    list_temp_ec.append(result_json[j]['transcripts'][0]['translations'][0]['xrefs'][i]['primary_id'])
                if 'Uniprot' in result_json[j]['transcripts'][0]['translations'][0]['xrefs'][i]['dbname']:
                    indice_uniprot = i
            # Put information in list only if an EC is found
            if len(list_temp_ec) != 0:
                for EC in list_temp_ec:
                    fh.write("\t".join((
                        result_json[j]['id'],  # gene ID
                        result_json[j]['transcripts'][0]['translations'][0]['xrefs'][indice_uniprot]['primary_id'],
                        # protein name
                        EC  # EC number)
                    )) + "\n")

            # Option, get attributes of the gene and/or the names
            if args.export_genes or args.export_names:
                gene_attributes, names = get_genes_attributes_and_names(result_json, j)
                if args.export_genes:
                    # write gene's attributes into file handle
                    fa.write(gene_attributes + "\n")
                if args.export_names:
                    # write names into file handle
                    nh.write(names + "\n")

        except TypeError:
            pass

    # Report execution time
    end_time = current_time()
    processing_time(start_time, end_time, fh, args.verbose)

    fh.close()

    if args.export_genes:
        fa.close()

    if args.export_names:
        nh.close()


def get_genes_attributes_and_names(result_json, index):
    """Get the genes attributes, the names and id for each gene.

    :param result_json: the JSON result
    :param index: (int) index of gene
    :return: (tupple) gene_attributes
    """

    id = check_if_none(result_json, index, "id")
    biotype = check_if_none(result_json, index, "biotype")

    # Check if there's a name. Otherwise, use the id.
    if result_json[index]["name"] != None:
        name = result_json[index]["name"]
    else:
        name = id

    if "ENA_FEATURE_GENE" in result_json[index]["xrefs"][0]["dbname"]:
        contig = result_json[index]["xrefs"][0]["primary_id"]
    else:
        contig = "NA"

    start = check_if_none(result_json, index, "start")
    end = check_if_none(result_json, index, "end")
    strand = check_if_none(result_json, index, "strand")
    description = check_if_none(result_json, index, "description")

    gene_attributes = ("\t".join((id, biotype,
                                  name,
                                  contig,
                                  start,
                                  end,
                                  strand,
                                  description)))

    names = ("\t".join((id, name)))

    return gene_attributes, names


def check_if_none(result_json, index, attribute):
    """ Check if the attribute is None.

    :param result_json: the JSON result
    :param index: (int) index of gene
    :param attribute: (str) attribute of the gene
    :return: NA if the attribute is None, the attribute otherwise
    """

    if result_json[index][attribute] != None:
        return result_json[index][attribute]

    return "NA"


def write_index(dictionary, file, species):
    """Write the content of dictionary into file

    :param dictionary: dictionary which values are the name of outfiles.
    :param file: index file
    :param species: name of the species
    :return: None
    """

    line = "| " + species + " | "

    for key in sorted(dictionary.keys()):
        if dictionary[key]:
            line += "[" + ntpath.basename(dictionary[key]) + "](" + os.path.join(os.getcwd(), dictionary[key]) + ") | "
        else:
            line += "Not exported"

    file.write(line + "\n")


if __name__ == '__main__':

    # Specify and docment the command line arguments
    parser = argparse.ArgumentParser(description="For a given species, collect all the " +
                                                 "relationships between genes and enzymatic functions (via EC number).",
                                     epilog="Jeanne Cheneby, Justine Long and Jacques van Helden",
                                     formatter_class=argparse.RawTextHelpFormatter)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-s", "--species", help="\n".join(("Species of interest.",
                                                          "example: escherichia_coli_str_k_12_substr_mg1655")))
    group.add_argument("-sf", "--species_file",
                       help="\n".join(("File containing species of interest.",
                                       "By default the species are expected to be found in the second column" +
                                       "(due to the organisation of the result file of ensemblgenomes_get_organisms.py",
                                       "but the query column can be changed with the option -c.")))

    parser.add_argument("-o", "--outdir",
                        help="Output directory." +
                             "Since the program can treat multiple organisms, and exports several files per organisms, "
                             + "the output file names are defined automatically. "
                             + "Each species is exported in a specific sub-directory named after the species name. ",
                        required=True,
                        action="store")

    parser.add_argument("-c", "--column_nb",
                        help="\n".join((
                            "Index of the column in which the species' name are found in the input file(species_file).",
                            "Is not needed if the input file was created by ensemblgenomes_get_organisms.py")),
                        type=int,
                        action="store",
                        default=2)

    parser.add_argument("-v", "--verbose",
                        help="\n".join((
                            "Verbosity. Supported: 0, 1, 2, 3. Default = 0.",
                            "Level 1: include the command line in the output.",
                            "Level 2: provide time warns about ongoing tasks.",
                            "Level 3: detailed information about the requests.")),
                        type=int,
                        choices=[0, 1, 2, 3],
                        default=0)

    parser.add_argument("-j", "--export_json", help="export JSON file",
                        action="store_true")

    parser.add_argument("-g", "--export_genes",
                        help="Create a file for each species containing extra information about genes",
                        action="store_true")

    parser.add_argument("-n", "--export_names",
                        help="Create a file for each species containing the name of the genes and their id.",
                        action="store_true")

    parser.add_argument("-a", "--export_all",
                        help="Export all the file than can be exported.",
                        action="store_true")

    args = parser.parse_args()

    # Check if all the file have to be exported
    if args.export_all:
        args.export_genes = True
        args.export_json = True
        args.export_names = True

    if args.species_file:
        species_list = parsing_file(args.species_file, args.column_nb)
    else:
        species_list = [args.species]

    index_file = open_output_handle(args.outdir + "/index.md", args.verbose)
    index_file.write("# Index \n")
    index_file.write("| Species |  GPE file  |  Gene attributes file  |  Gene names file  |  Gene (json) file  |  \n")
    index_file.write("|---------|:----------:|:----------------------:|:-----------------:|:------------------:| \n")

    for species in species_list:
        output_files, start_time = initializing_dir(species)
        get_GPE(species, output_files, start_time)
        write_index(output_files, index_file, species)

    index_file.write("|||\n")

    index_file.close()

    verbosity(current_time() + "\t" + "Job done", 2, args.verbose)
