from __future__ import print_function
from itertools import chain
from pip._vendor.distlib.compat import raw_input
from datetime import datetime
import requests
import pprint
import os
import sys
import json
import math
import time
try:
    from urllib.request import urlopen
except ImportError:
    from urllib import urlopen

from datetime import datetime

def warning(message, display_level=0, verbosity_level=0):
    """Display a warning message.

    :param message: (str) information about the current steps
    :param display_level: (int) level above which the message has to be displayed.
    :param verbosity_level: verbosity level
    :return: None
    """
    if verbosity_level >= display_level:
        verbosity("Warning\t" + message + "\n", display_level, verbosity_level)

def current_time():
    """Return the current time.

    :return: (datetime) the current datetime
    """
    time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    return time

def download_go(outfile):
    """Download the term definition of the Gene Ontology from the official
    Web site, and store it in a local file.

    :param outfile:output file (default: GO_hierarchy.obo)
    :return: no return value. The download result is printed to the output file.

    """
    #url = "http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo"
    url = "http://current.geneontology.org/ontology/go.obo"

    if not os.path.isfile(outfile):
        print("Downloading gene ontology terms from \n\t{} ...".format(url))
        try:
            response = urlopen(url)
            data = response.read()
            str_data = data.decode('utf-8')
            print("\tDone\nGene ontology file in obo format:\n\t{}".format(outfile))
        except:
            print("Can't reach server at : %s" % url)
            print(sys.exc_info()[1])
            sys.exit(-1)
    else:
        print("File {} already exists.".format(outfile))
        sys.exit(0)

    try:
        with open(outfile, "w") as f:
            f.write(str_data)
        print("\tFile successfully created\n")
    except:
        print("Can't open file {}".format(outfile))
        print(sys.exc_info()[1])


def after_colon(line):
    # macro for getting anything after the :
    return line.split(":", 1)[1].strip()


def run_parse_go(filename):

    print("Start parsing the file {}".format(filename))

    list_rec = []
    rec = GOTerm()
    with open(filename, 'r') as f:

        for line in f:
            if line.startswith("[Term]"):
                if rec.id != "" and rec.is_obsolete is False:
                    list_rec.append(rec)
                rec = GOTerm()
            if line.startswith("[Typedef]"):
                break
            if line.find("-namespace: ") != -1:
                continue
            if line.startswith("id:"):
                rec.id = after_colon(line)
            if line.startswith("alt_id:"):
                rec.alt_ids.append(after_colon(line))
            elif line.startswith("name:"):
                rec.name = after_colon(line)
            elif line.startswith("namespace:"):
                rec.namespace = after_colon(line)
            elif line.startswith("xref"):
                rec.xrefs.append(after_colon(line).split()[0])
            elif line.startswith("is_a:"):
                rec.parents.append(after_colon(line).split()[0])
            elif (line.startswith("is_obsolete:") and
                  after_colon(line) == "true"):
                rec.is_obsolete = True
    #last term
    list_rec.append(rec)
    print("\tParsing done\n")
    return list_rec


def file_creation(outputdir, list_rec):
    #create files
    print("At {} :".format(outputdir))
    print("\tCreation of file GO_xrefs.tab")
    print("\tCreation of file GO_relations.tab")
    print("\tCreation of file GO_description.tab\n")
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)
    with open(os.path.join(outputdir, "GO_description.tab"), 'w') as f1:
        with open(os.path.join(outputdir, 'GO_relations.tab'), 'w') as f2:
            with open(os.path.join(outputdir, 'GO_xrefs.tab'), 'w') as f3:
                f1.write("#GO_ID" + '\t' + "GO Term" + '\t' + "Ontology Type" + '\t' + "GO_ALT_ID" + '\n')
                f2.write("GO_Child_ID" + '\t' + "GO_Parent_ID" + '\n')
                f3.write('#GO_ID' + '\t' + 'Source_GO_Xref' + '\t' + 'ID_source_GO_Xref' + '\n')

                for term in list_rec:
                    f1.write(term.id + '\t' + term.name + '\t' + term.namespace + '\t' + ', '.join(term.alt_ids) +'\n')

                    if term.xrefs:
                        for xref in term.xrefs:
                            separate_xref = xref.split(":")
                            f3.write(term.id + '\t' + separate_xref[0] + '\t' + separate_xref[1] + '\n')
#f3.write(term.id + '\t' + ', '.join(term.xrefs) + '\n')
                    if len(term.parents) > 0:
                        for i in range(len(term.parents)):
                            f2.write(term.id + '\t' + term.parents[i] + '\n')
                    else:
                        f2.write(term.id + '\t' + 'None' + '\n')

    return


def parse_go(f, outputdir):

    print("\n=== Extracting informations from an OBO file ===\n")

    if f:
        print("Start parsing file\n\t".format(f))
        if os.path.isfile(f):
            file_creation(outputdir, run_parse_go(f))

        else:
            print("Specified file path is incorrect, no such file found")
            if raw_input("Would you like to download it ? (y/n) ") == "y":
                download_go(f)
                file_creation(outputdir, run_parse_go(f))
            else:
                sys.exit(0)
    else:
        if os.path.isfile('GO_hierarchy'):
            print("Found a file 'GO_hierarchy' in current folder.")
            if raw_input("Would you like to parse it ? (y/n) ") == "y":
                file_creation(outputdir, run_parse_go('GO_hierarchy'))
            else:
                if raw_input("Would you like to download the required file ? (y/n) ") == "y":
                    filename = "__GO_hierarchy__"
                    while os.path.isfile(filename):
                        filename += "_"
                    download_go(filename)
                    file_creation(outputdir, run_parse_go(filename))
                else:
                    sys.exit(0)
        else:
            filename = "GO_hierarchy"
            download_go(filename)
            file_creation(outputdir, run_parse_go(filename))

def open_output_handle(file_name=None, verbosity_level=2):
    """Open a write handle either on a file (if the argument file_name is specified) or on the STDOUT device.

    :param file_name: name of the output file. If not specified, a handle on the STDOUT is returned.
    :param verbosity_level: level of verbosity from which the optional message should be displayed.
    :return: a writable file handle
    """
    # check if outfile name and directory already exist
    if file_name:
        # Check the directory and create it if required
        file_dir = os.path.dirname(file_name)
        if file_dir != '':
            # Check if directory exists
            if os.path.exists(file_dir):
                # If a file exists in the place of the directory, die !
                if not os.path.isdir(file_dir):
                    raise Exception("Cannot create directory " + file_dir + " file exists with same name")

            # Create directory if doesn't exist
            else:
                os.makedirs(file_dir)

        # Check if output file already exists
        if os.path.exists(file_name):
            warning("Overwriting existing file " + file_name, 2, verbosity_level)

        actual_time = current_time()
        verbosity(actual_time + "\t" + "Opening file " + file_name, 2, verbosity_level)

        try:
            fh = open(file_name, 'w')
        except PermissionError:
            sys.stderr.write("Permission error: can't create: " + file_name +
                             "\n File may be open in an other program.")
            sys.exit()

    else:
        fh = sys.stdout

    return fh

def verbosity(message, display_level, verbosity_level=0):
    """Display the message and the time if verbosity level is superior to 1.

    :param message: (str) information about the current steps
    :param display_level: (int) level above which the message has to be displayed.
    :param verbosity_level: verbosity level
    :return: None
    """
    if verbosity_level >= display_level:
        sys.stderr.write(message + "\n")

def treat_one_request(database, species, outputdir):
    """Treat one REST request on a given database ("ensemblgenomes" or "ensembl"),
    and return the result as a list or dictionary produced by applying decode on the JSON result.

    :param database:  (str) database name ("ensemblgenomes" or "ensembl").
    :param verbosity_level: verbosity level
    :return: a list or dictionary produced by applying the function json() on the result
    """

    if database == "ensembl":
        server = "http://rest.ensembl.org"

    else:
        server = "http://rest.ensemblgenomes.org"

    # Define extension for the REST URL
    ext = "/lookup/genome/" + species + "?level=protein_feature&xrefs=1/"

    # Send the request to the REST server
    rest_url = server + ext + "?"

    verbosity(current_time() + "\t" + "Getting result from " + server +
              " REST URL (" + rest_url + "content-type=application/json" + ")" + "\n", 2, 2)

    try:
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    # Treat connection error
    except requests.exceptions.ConnectionError:
        sys.stderr.write("ConnectionError: Can't connect to server: " + server + "\n")

    # Check the validity of the result
    if not r.ok:
        if r.status_code == 400:
            sys.stderr.write("ArgumentError: URL invalid" + "\n")
            sys.exit()

    verbosity(current_time() + "\t" + "Got result from " + server +
              " REST URL (" + rest_url + "content-type=application/json" + ")" + "\n", 2, 2)

    jh = open_output_handle(outputdir + "/annotations.json")
    pprint.pprint(r.json(), stream=jh)
    jh.close()

    parsing_results_list = parsing_json(r.json())

    create_file_gene_go(parsing_results_list, species, outputdir)

def connection_server_gene(org_name):

    """Collect information abot all the genes of the specified organism
    via the REST interface of Ensemblgenomes.

    @param org_name: organims name. Must be an organism supported at Ensemblgenomes.
    @param v    : verbosity index
    @rtype      : Json object
    @return     : The server's response

    """

#    url_prefix = "http://beta.rest.ensemblgenomes.org"
    url_prefix = "http://rest.ensemblgenomes.org"
    url_suffix = "/lookup/genome/" + org_name + "?content-type=application/json"
    url = url_prefix + url_suffix

    try:
        print("Trying to reach ensembl server at \n\t%s ..." % url)
        response = urlopen(url)
        data = response.read()
        decoded_data = data.decode('utf-8')
    except:
        print("Can't reach server at %s" % url)
        print(sys.exc_info()[1])
        sys.exit(-1)

    print("\tRespond received\n")
    return decoded_data


def parsing_json(j_obj):
    #Goal of this function : return a list with list[0] = a list of gene id and list[1] = a list of GO di associated

    print("Start parsing of JSON file")

    x_list = []
    gene_id_list = []

    for value in j_obj:
        current_gene = value
        current_gene_id = current_gene['id']

        transcripts_list = current_gene['transcripts']
        for t_value in transcripts_list:
            current_transcript = t_value

            translations_list = current_transcript['translations']
            if translations_list is not None:
                for tr_value in translations_list:
                    current_translation = tr_value

                    xrefs_list = current_translation['xrefs']
                    x_current_list = []
                    for x in xrefs_list:
                        current_xref = x
                        if "GO:" in current_xref['display_id']:
                            x_current_list.append(current_xref['display_id'])
                    if x_current_list:
                        gene_id_list.append(current_gene_id)
                        x_list.append(x_current_list)

    print("\tParsing done")
    return [gene_id_list, x_list]


def create_file_gene_go(l, org_name, outdir):

    filename = "annotations_table_" + org_name + ".tab"
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    with open(os.path.join(outdir, filename), 'w') as f:
        f.write("Gene_ID" + '\t' + "GO_ID" + '\n')

        for i in range(len(l[0])):
            counter = 0
            while counter < len(l[1][i]):
                f.write(l[0][i] + '\t' + l[1][i][counter] + '\n')
                counter += 1

    print("At {} :".format(outdir))
    print("\tCreation of file {}".format(filename))


def get_annotations(org_name, outputdir):

    json_response = connection_server_gene(org_name)

    try:
        decoded_json = json.loads(json_response)
    except ValueError:
        print("The Json file received is corrupted or invalid")
        print(sys.exc_info()[1])
        sys.exit(-1)

    parsing_results_list = parsing_json(decoded_json)

    create_file_gene_go(parsing_results_list, org_name, outputdir)



class GOTerm:
    """Each instance of this class corresponds to a GO term.

    The GO term contains information about its direct relatives:
    children and parents. It is equipped with ..
    """

    def __init__(self):
        self.id = ""                    # GO:xxxxxx
        self.name = ""                  # description
        self.namespace = ""             # BP, CC, MF
        self.parents = []               # parent records
        self.children = []              # children records
        self.level = -1                 # distance from root node
        self.is_obsolete = False        # is_obsolete
        self.alt_ids = []               # alternative identifiers
        self.xrefs = []                 # cross references
        self.genes = set()              # Gene name
        self.expanded_genes = set()     # All genes with child's one
        self.visited_flag = False       # Set to True when node already explored


class DAG:
    """
   DOC TO BE ADDED
    """

    def __init__(self):
        """
    DOC TO BE ADDED
        """
        
        self.Roots = []
        self.Dict = dict()
        self.level = 0

    def load_from_file(self, filename):

        print("\nLoading GO relations from file " + filename + "\n")
        with open(filename, "r") as f:
            l = 0
            rel = 0
            for line in f:
                l += 1

             ## Read child and parent
                (child, parent) = line.rstrip().split("\t")[0:2]

                if child in self.Dict and not (line.startswith(';') or line.startswith('#')):
                    rel += 1
                    if parent == "None":
                        sys.stderr.write("\t".join(["Root term", child, self.Dict[child].name]) + "\n")
                        self.Roots.append(child)
                    elif parent in self.Dict:
                        self.Dict[child].parents.append(self.Dict[parent])
                        self.Dict[parent].children.append(self.Dict[child])
                    else:
                        sys.stderr.write("\t".join(["line", str(l), "skipping invalid parent GO term :", "'" + parent + "'"]) + "\n")
                else:
                    sys.stderr.write("\t".join(["line", str(l), "skipping invalid child GO term :", "'" + child + "'"]) + "\n")
        print("\nFinished loading\t {} GO relations found\n".format(rel))

    def expand_annot(self, annot_file):
        current_gene_dict = dict()
        gene = ""
        old_gene = -1
        entries = 0

        print("\nBegin expansion of file \n\t" + annot_file + "\n")

        with open(annot_file, 'r') as a:
            with open('expanded_' + annot_file, 'w') as exp_a:
                for line in a:
                    if gene:
                        old_gene = gene
                    (gene, go) = line.rstrip().split("\t")[0:2]

                    if old_gene != -1 and old_gene != gene:
                        for go_id in current_gene_dict:
                            entries += 1
                            exp_a.write(old_gene + '\t' + go_id + '\n')
                        current_gene_dict = dict()

                    #if it is an exploitable line:
                    if go in self.Dict and not (line.startswith(';') or line.startswith('#')):
                        current_gene_dict = dict(chain(parent_recursion([self.Dict[go]], current_gene_dict, gene).items(), current_gene_dict.items()))
                #witring the last dictionary
                for go_id in current_gene_dict:
                    entries += 1
                    exp_a.write(old_gene + '\t' + go_id + '\n')
        print("Expand done\t%d entries written\n" % entries)

## End of DAG class definition
################################################################


def parent_recursion(list_child, current_gene_dict, gene):
    if list_child:
        next_list = []
        for child in list_child:
            current_gene_dict[child.id] = child
            for parent in child.parents:
                if not gene in parent.genes:
                    parent.genes = parent.genes.union(gene)
                    next_list.append(parent)
        return parent_recursion(next_list, current_gene_dict, gene)
    else:
        return current_gene_dict


def create_go_dictionary(desc_file):
    d = dict()
    try:
        with open(desc_file, "r") as f:
            f.readline()
            for line in f:
                first_tab = line.find('\t')
                second_tab = line.find('\t', first_tab + 1)
                third_tab = line.find('\t', second_tab + 1)
                GO = GOTerm()
                GO.id = line[:first_tab]
                GO.alt_ids = line[first_tab + 1:second_tab]
                GO.name = line[second_tab + 1:third_tab]
                GO.namespace = line[line.find('\t', 11) + 1:]
                d[GO.id] = GO
            return d
    except:
        sys.stderr.write(str(sys.exc_info()[1]))
        sys.exit(2)


def expand_parsing_check(o, d, r, a):
    if o:
        if d or r or a:
            sys.stderr.write("Choose between [-org <org>] OR [-d <desc_file> -r <rel_file> -a <annot_file>]\n")
            sys.exit(0)
        #Check if annotation file exists for this organism
        if not (os.path.exists("GO_relations.tab") and os.path.exists("GO_description.tab")):
            if os.path.exists('GO_hierarchy'):
                file_creation(os.getcwd(), run_parse_go('GO_hierarchy'))
            else:
                parse_go(None, outputdir=os.getcwd())

        if not os.path.exists("annotations_table_" + o + ".tab"):
            get_annotations(o, os.getcwd())

    else:
        if not os.path.exists(d):
            sys.stderr.write("file {} does not exists at specified location.\n".format(d))
            return -1
        if not os.path.exists(r):
            sys.stderr.write("file {} does not exists at specified location.\n".format(r))
            return -1
        if not os.path.exists(a):
            sys.stderr.write("file {} does not exists at specified location.\n".format(a))
            return -1

    return 0


def expand(org, desc_file, rel_file, annot_file):

    err = expand_parsing_check(org, desc_file, rel_file, annot_file)

    if err != 0:
        sys.exit(-1)

    if org:
        desc_file = 'GO_description.tab'
        rel_file = 'GO_relations.tab'
        annot_file = 'annotations_table_' + org + '.tab'

    ## DAG construction
    diag = DAG()
    dict_go = create_go_dictionary(desc_file)
    diag.Dict = dict_go
    diag.load_from_file(rel_file)

    #Annotations table go-through
    diag.expand_annot(annot_file)


def enrichment():
    print("function enrichment")

