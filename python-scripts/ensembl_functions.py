#!/usr/bin/python
#-*-coding: utf-8-*-

### imports ###
#Used to use print 3.x function
from __future__ import print_function
import os
import sys
import json
import platform
#if using python 3.x
try:
    from urllib.request import urlopen
#else
except ImportError:
    from urllib import urlopen


python_version = platform.python_version_tuple()

### functions needed for retrieve Species###


def connection_server_species(v, dump, division):

    """
    This function will send a request to a server and then return the answer
    @param v: verbosity index
    @return: A Json object. The response from the server to the request
    """

    url = "http://beta.rest.ensemblgenomes.org/info/species?"
#    url = "http://test.rest.ensemblgenomes.org/info/species?"
    url = url + "content-type=application/json"
    if (division != None):
        url = url + ";division="+division

    if v >= 2:
        print("Sending a request to the http server...")
        print(url)

    try:
        response = urlopen(url)
        data = response.read()
        decoded_data = data.decode('utf-8')
        json_dict = json.loads(decoded_data)

        if dump:
            with open('dump_species', 'w') as d:
                d.write(json.dumps(json_dict, indent=4))

        if v >= 3:
            print("Server sent a response")
    except:
        print("Can't reach server at : %s" % url)
        print(sys.exc_info()[1])
        if os.path.exists('org'):
            os.remove('org')
        sys.exit()

    return json_dict


def parsing_content_species(org_supported, v):

    """
    This function will parse a Json object into a string and then return it
    @param org_supported_json: A json object, it is the response from a server to a html request
    @param v: vrebosity index
    @return: String object. Result of the parsing of the Json object
    """

    out_columns = ['division', 'display_name', 'name', 'common_name', 'assembly', 'groups', 'taxon_id', 'release', 'aliases']

    ## Print the header
    file_content += "#"
    file_content += "\t".join(out_columns)
    file_content += '\n'

    if v >= 3:
        print('Parsing the server response (Json file)')

    #org_supported = json.loads(org_supported_json)

    file_content = ""
    if v >= 3:
        print('Retrieving data, this may take a while...')


    for value in org_supported.values():
        ## Value contained in the dictionary is a list of dictionaries containing the species informations
        list_species_dict = value


        for i in range(len(list_species_dict)):
            ## Dictionary containing informations about a given species
            current_species_dict = list_species_dict[i]

            ## The value of current_species_dict[u'groups'] is a list:
            ## the for loop allows to display each element of this list without
            ## brackets and separated by commas
            the_groups = ""
            for group in current_species_dict[u"groups"]:
                the_groups = the_groups + group + ',' + ' '
                current_species_dict[u'groups'] = the_groups

            ## The value of current_species_dict[u'aliases'] is a list:
            ## the for loop allows to display each element of this
            ## list without brackets and separated by commas
            the_aliases = ""
            for alias in current_species_dict[u'aliases']:
                the_aliases = the_aliases + alias + ',' + ' '
                current_species_dict[u'aliases'] = the_aliases

            ## Feeding of file_content variable
            values = []
            for j in range(len(out_columns)):
                if out_columns[j] in current_species_dict.keys():
                    value = str(current_species_dict[out_columns[j]])
                    values.append(value)
#                    file_content += str(current_species_dict[out_columns[j]]) + '\t'
            file_content += "\t".join(values)
            file_content += '\n'

    if v >= 2:
        print("Parsing done\n")

    return file_content


def create_file_species(ff, o, v):

    """This function will create a file filled with the informations
    about all the species collected on ensemblGenomes

    @param ff: String containing all of the informations about the
    species

    @param o: output file name

    @param v: verbosity index

    """

    ## Default directory to store species list in case the -o option
    ## is not specified
    parent_folder = 'species_results'
    file_name = o

    ## If output is not a filename but a path:
    if not os.path.dirname(o) == '':
        parent_folder = os.path.dirname(o)
        file_name = o[len(parent_folder)+1:]

    if v >= 2:
            print("Creation of a new folder named {}".format(parent_folder))
    if not os.path.exists(parent_folder):
        os.makedirs(parent_folder)

    if v >= 3:
        print("Creation of a new file named {}".format(file_name))
    with open(os.path.join(parent_folder, file_name), "w") as newfile:
        newfile.write(ff)

    if v >= 3:
        print("File created\n")


def create_name_only_file(o, v, f):

    """
    This function will create a file containing the name the species in ensemblGenomes. This file can be used directly
    in the retrieveGenes function as an input of the '-f' option.
    @param o: output file name
    @param v: verbosity index
    @param f: name of the file containing species name only
    """
    parent_folder = 'species_results'
    file_name = o

    ## If output is not a filename but a path:
    if not os.path.dirname(o) == '':
        parent_folder = os.path.dirname(o)
        file_name = o[len(parent_folder)+1:]

    i = 0
    stri = ""
    content = ""
    if v >= 2:
        print("Trying to create a file with the name of species alone")
    if v >= 3:
        print("Searching species informations at {}...".format(os.path.join(os.getcwd(), parent_folder, file_name)))

    try:
        with open(os.path.join(os.getcwd(), parent_folder, file_name), 'r') as fi:
            if v >= 3:
                print("File {} found".format(o), "Parsing the file...")
            for line in fi:
                text = line

                for char in text:
                    if char == '\t':
                        i += 1
                    if i == 2:
                        if char == '\t':
                            continue
                        stri += char
                i = 0
                content += stri + '\n'
                stri = ""
        if v >= 3:
            print("Parsing done.`\nCreating the file {}".format(f))
        with open(os.path.join(parent_folder, f), 'w') as fi:
            fi.write(content)

        if v >= 2:
            print("File created\n")
    except IOError:
        print("Error: File {} cannot be found in directory {}\n".format(o, os.path.join(os.getcwd(), parent_folder)))

### end of functions needed for retrieve Species###


def retrieve_species(o, v, f, dump, null, division):

    """
    This function will create a file containing all the species of ensemblGenomes and the informations about it. If the
    -f option is enabled, it will also create a file containing only the name of each species so that the file can be
    used as an input to the retrieveGenes function.
    @param o: output file name
    @param v: verbosity index
    @param f: name of the file containing species name only
    """

    if v >= 1:
        print("Connecting server ...")
    org_supported_json = connection_server_species(v, dump, division)

    if v >= 1:
        print("Server's response received\n\nParsing informations...")
    ff = parsing_content_species(org_supported_json, v)

    if v >= 2:
        print("Creating new file with informations...")
    create_file_species(ff, o, v)

    if f is not None:
        if v >= 2:
            print("Creating file with specie's name only..")
        create_name_only_file(o, v, f)

    if v >= 1:
        print("\nWork done, the file {} has been created".format(o))
    if f is not None and v >= 1:
        print("the file {} has been created".format(f))


### functions needed for retrieve Genes ###

def connection_server_genes(line, v):

    """
    This function connect with the ensemblGenomes server and bring back the genes for a given specie
    @param line : String of char. It is also a specie name
    @param v    : verbosity index
    @rtype      : Json object
    @return     : The server's response
    """

    urlprefixe = "http://beta.rest.ensemblgenomes.org"
    urlsuffixe = "/lookup/genome/" + line + "?content-type=application/json"
    url = urlprefixe + urlsuffixe

    try:
        if v == 3:
            print("Trying to reach server...")
        response = urlopen(url)
        data = response.read()
        decoded_data = data.decode('utf-8')
        if v == 3:
            print("Response received")
    except:
        print("Can't reach server at %s" % url)
        print(sys.exc_info()[1])
        if os.path.exists('org'):
            remove('org')
        sys.exit()

    return decoded_data


def parsing_content_genes(decoded_json, v, null):
    """
    This function will go through the decoded_json object and extract and format every information
    by parsing and separate each sens unit by a '\t' char
    @param decoded_json: Deserialized Json object turned into a Python object following this table :
    ####################    JSON    ######    Python     ###################
    #                     object      ->      dict                         #
    #                     array       ->      list                         #
    #                     string      ->      str                          #
    #                     number(int) ->      int, long                    #
    #                     number(real)->      float                        #
    #                     true        ->      True                         #
    #                     false       ->      False                        #
    #                     null        ->      None                         #
    ########################################################################
    @param v: verbosity index
    @return: a tuple (gene_content, transcript_content, protein_feature_content, xrefs_content, linkage_types_content)
    It contains all the content for the future files
    """

    ##  Declaration of the columns of each files
    columns_key_genes = ['id',
                         'seq_region_name',
                         'start',
                         'end',
                         'strand',
                         'name',
                         'biotype',
                         'description']

    columns_key_transcripts = ['id',
                               'seq_region_name',
                               'start',
                               'end',
                               'strand',
                               'name',
                               'biotype',
                               'description',
                               'id_gene']

    columns_key_protein_features = ['id_transcript',
                                    'id_gene/id_translation',
                                    'dbname',
                                    'start',
                                    'end',
                                    'name',
                                    'interpro_ac',
                                    'description']

    columns_key_protein_xrefs = ['id_transcript',
                                 'id_gene/id_translation',
                                 'dbname',
                                 'display_id',
                                 'primary_id',
                                 'linkage_types']

    columns_key_linkage_types = ['GO_id',
                                 'id_transcript',
                                 'id_gene',
                                 'evidence',
                                 'dbname',
                                 'display_id',
                                 'primary_id']

    ##  ---------------- Creating all of the empty content ---------------- ##
    gene_content = ""
    transcript_content = ""
    protein_feature_content = ""
    xrefs_content = ""
    linkage_types_content = ""

    ##  ---------------- Beginning of the main loop which parse the content of the server's response ---------------- ##

    if v == 3:
        print("Parsing begin...")
    ##  The for permit to get into the decoded_json (type(decoded_json) = <list>)
    for value in decoded_json:
        ##  Extracting the dictionnary inside the list
        current_gene_dict = value

        ##  Matching all attributes of the dictionnary with its value and writing it into the specific file
        gene_attributes = []
        for i in range(len(columns_key_genes)):
            attr_key = columns_key_genes[i]
            if attr_key in current_gene_dict.keys():
                attr_value = current_gene_dict[columns_key_genes[i]]

                ## Check if gene name is defined. If not, use gene ID instead.
                if attr_key == 'name' and attr_value is None:
                    attr_value = current_gene_dict['id']
                if attr_key == 'description' and attr_value is None:
                    attr_value = null
                gene_attributes.append(str(attr_value))

        ## Append the line for the parsed gene
        gene_content += "\t".join(gene_attributes)
        gene_content += "\n"

        transcript_list = current_gene_dict['transcripts']
        for t_value in transcript_list:
            current_transcripts_dict = t_value

            ##  Matching all attributes of the dictionnary with its value and writing it into the specific file
            transcript_attributes = []
            for j in range(len(columns_key_transcripts)):
                attr_key = columns_key_transcripts[j]
                if attr_key in current_transcripts_dict.keys():
                    attr_value = current_transcripts_dict[columns_key_transcripts[j]]

                    ## Check if transcript name is defined. If not, use transcript ID instead.
                    if (attr_key == 'name') and (attr_value is None):
                        attr_value = current_transcripts_dict['id']
                    if attr_key == 'interpro_ac' and attr_value == "":
                        attr_value = null
                    if attr_key == 'description' and attr_value is None:
                        attr_value = null

                    transcript_attributes.append(str(attr_value))
            transcript_attributes.append(str(current_gene_dict['id']))
            ## Append the line for the parsed transcript
            transcript_content += "\t".join(transcript_attributes)
            transcript_content += '\n'

            ##  Extracting a list inside the dictionnary

            if 'translations' in current_transcripts_dict.keys():
                translation_list = current_transcripts_dict['translations']

                ##  Getting into the list
                for tr_value in translation_list:
                    ## Extracting the dictionnary inside the list
                    current_translations_dict = tr_value

                    ## This attribute is not in each translation.
                    if 'protein_features' in current_translations_dict.keys():
                        protein_features_list = current_translations_dict['protein_features']

                        for prot_f in protein_features_list:
                            current_protein_dict = prot_f
                            protein_feature_content += str(current_translations_dict['id']) + '\t'
                            protein_feature_content += str(current_gene_dict['id']) + '\t'

                            ## Matching all attributes of the dictionnary with its value and writing it into the specific file
                            for k in range(len(columns_key_protein_features)):
                                if columns_key_protein_features[k] in current_protein_dict.keys():
                                    protein_feature_content += str(current_protein_dict[columns_key_protein_features[k]]) + '\t'
                            protein_feature_content += '\n'

                    xrefs_list = current_translations_dict['xrefs']
                    for x in xrefs_list:
                        ## Extracting the dictionnary inside the list
                        current_xrefs_dict = x
                        xrefs_content += str(current_gene_dict['id']) + '\t'
                        xrefs_content += str(current_translations_dict['id']) + '\t'

                        ## Matching all attributes of the dictionnary with its value and writing it into the specific file
                        for l in range(len(columns_key_protein_xrefs)):
                            if columns_key_protein_xrefs[l] in current_xrefs_dict.keys() and columns_key_protein_xrefs[l] != 'linkage_types':
                                xrefs_content += str(current_xrefs_dict[columns_key_protein_xrefs[l]]) + '\t'

                        ## This attribute is not in each xref
                        if 'linkage_types' in current_xrefs_dict.keys():
                            xrefs_content += 'YES' + '\t'

                            linkage_types_list = current_xrefs_dict['linkage_types']
                            for link in linkage_types_list:
                                current_linkage_types_dict = link
                                linkage_types_content += str(current_xrefs_dict['display_id']) + '\t'
                                linkage_types_content += str(current_transcripts_dict['id']) + '\t'
                                linkage_types_content += str(current_gene_dict['id']) + '\t'

                                ## Matching all attributes of the dictionnary with its value and writing it into the specific file
                                for m in range(len(columns_key_linkage_types)):
                                    if columns_key_linkage_types[m] in current_linkage_types_dict.keys():
                                        linkage_types_content += str(current_linkage_types_dict[columns_key_linkage_types[m]]) + '\t'

                                    current_source_dict = current_linkage_types_dict['source']
                                    if columns_key_linkage_types[m] in current_source_dict.keys():
                                        linkage_types_content += str(current_source_dict[columns_key_linkage_types[m]]) + '\t'

                                linkage_types_content += '\n'
                        else:
                            xrefs_content += 'NO' + '\t'
                        xrefs_content += '\n'
    if v == 3:
        print("Parsing done")
    return gene_content, transcript_content, protein_feature_content, xrefs_content, linkage_types_content


def make_file(specie_name, arg_list, file_name, file_content, v, outputdir):
    """
    This little function just create a file named <file_name> and fill it with the <file_content>
    and the arg_list on top of the file
    @param specie_name: The name of the parent folder
    @param arg_list: the title of each column into a list
    @param file_name: the file output name
    @param file_content: the content of the file
    @param v: verbosity index
    """
    if outputdir is not None:
        f = open(os.path.join(outputdir, specie_name, file_name), 'w')
    else:
        if not os.path.exists('features_results'):
            if v >= 2:
                print("Creating new folder named features_results")
            os.makedirs('features_results')
        f = open(os.path.join('features_results', specie_name + '_' + file_name), 'w')

    f.write('#' + arg_list[0] + '\t' + '\t'.join(arg_list[1:]) + '\n')
    f.write(file_content)
    f.close()

    if v == 3:
        print("file {} created".format(file_name))


def file_folder_manager(specie_name, depth, parsed_tuple, v, outputdir):
    """
    This function will create the root folder 'Species Informations' and the sub folders with the species name on each.
    Then call the function that create the files with the appropriate args depending the user # of file wanted and name of file chosen.
    @param specie_name: the name of the folder
    @param depth: the number of files the user wants
    @param parsed_tuple: a tuple with the content of each file in it
    @param v: verbosity index
    @param outputdir: if not None, creation of a folder tree
    """

    if outputdir is not None:
        species_dir = os.path.join(outputdir, specie_name)
    else:
        species_dir = specie_name

    columns_key_genes = ['id',
                         'seq_region_name',
                         'start',
                         'end',
                         'strand',
                         'name',
                         'biotype',
                         'description']

    columns_key_transcripts = ['id',
                               'seq_region_name',
                               'start',
                               'end',
                               'strand',
                               'name',
                               'biotype',
                               'description',
                               'id_gene']

    columns_key_protein_features = ['id_transcript',
                                    'id_gene/id_translation',
                                    'dbname',
                                    'start',
                                    'end',
                                    'name',
                                    'interpro_ac',
                                    'description']

    columns_key_protein_xrefs = ['id_transcript',
                                 'id_gene/id_translation',
                                 'dbname',
                                 'display_id',
                                 'primary_id',
                                 'linkage_types']

    columns_key_linkage_types = ['GO_id',
                                 'id_transcript',
                                 'id_gene',
                                 'evidence',
                                 'dbname',
                                 'display_id',
                                 'primary_id']

    ## If the option --outputdir is active...
    if outputdir is not None:

        ## Creation of the root folder if it does not exist
        if not os.path.exists(outputdir):
            if v >= 2:
                print("Creating new folder named {}".format(outputdir))
            os.makedirs(outputdir)

        ## Creation of the sub folder if it does not exist
        if not os.path.exists(species_dir):
            os.makedirs(species_dir)

    ## Initialisation of default name list
    default_name_output_list = ['genes.tab', 'transcripts.tab', 'translations.tab', 'proteins_xrefs.tab', 'linkage_types.tab']
    ## The column are placed into a list to ease the use of the for loop below
    columns_keys = [columns_key_genes, columns_key_transcripts, columns_key_protein_features, columns_key_protein_xrefs, columns_key_linkage_types]

    ## Creating all files with default names
    for i in range(depth):
        make_file(specie_name, columns_keys[i], default_name_output_list[i], parsed_tuple[i], v, outputdir)


def retrieve_features_file(file_name, depth, v, outputdir, dump, null):

    """
    This is the main function. file_name is a file wheter the user used -f or -org. That simplify the code by far.
    This function will call all of the others above.
    @param file_name: the file containing all species' name
    @param depth: the # of file the user wants
    @param v: verbosity index
    @param outputdir: if outputdir != None, a tree will be created with a folder for each specie
    """

    with open(file_name, 'r') as f:
        ## Count the number of lines in f
        c = 0
        for lines in f:
            if (lines[0] != '#') and (lines[0] != "\n"):
                c += 1
        lines_number = c

    i = 0
    with open(file_name, 'r') as f:
        for line in f:
            if line[0] != '#' and line[0] != '\n':
                percent = float(i)/lines_number * 100

                if v >= 1:
                    print("\n\tWorking on : {}.\n\tDone: {}/{} ".format(line[:-1], i, lines_number), end="")
                    print("({0:.2f}%)".format(percent), end="\n\n")

                i += 1
                if v >= 2:
                    print('Connection...')

                json_ = connection_server_genes(line[:-1], v)
                if json_ == "":
                    decoded_json = ""
                else:
                    try:
                        decoded_json = json.loads(json_)
                        if dump:
                            with open('dump_features', 'w') as d:
                                d.write(json.dumps(decoded_json, indent=4))

                    except ValueError:
                        print("The Json file received is corrupt or invalid")
                        print(sys.exc_info()[1])
                        os.remove('org')
                        sys.exit()

                if v >= 2:
                    print('\nParsing:')
                parsed_content_tuple = parsing_content_genes(decoded_json, v, null)

                if outputdir is not None and v >= 2:
                    print('\nCreating folder and files...')
                elif v >= 2:
                    print('\nCreating files...')

                file_folder_manager(line[:-1], depth, parsed_content_tuple, v, outputdir)

                if v >= 2:
                    print('Done')


def retrieve_features(org, f, depth,  v, outputdir, dump, null):

    """
    This function will call retrieve_features_file. If -org option was chosen, the specie's name is put into a blank file
    and then the retrieve_features_file function is called.
    @param org: a specie's name or None if -f option was chosen
    @param f: a file containing several species' name or None if -org option was chosen
    @param depth: # of file the user wants
    @param verbose: verbosity index
    @param outputdir: if outputdir != None, a tree will be created with a folder for each specie
    """

    if org:
        if v == 3:
            print("\nSingle organism mode selected\n")
        with open('org', 'w') as org_file:
            org_file.write(org + '\n')
        retrieve_features_file('org', depth, v, outputdir, dump, null)
        os.remove('org')
    elif f:
        if v == 3:
            print("\nMultiple organisms mode selected\n")
        assert os.path.exists(f)
        retrieve_features_file(f, depth, v, outputdir, dump, null)
