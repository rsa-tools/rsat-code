#!/usr/bin/python
# -*-coding: utf-8-*-

import requests, sys, os
from datetime import datetime


def current_time():
    """Return the current time.

    :return: (datetime) the current datetime
    """
    time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    return time


def warning(message, display_level=0, verbosity_level=0):
    """Display a warning message.

    :param message: (str) information about the current steps
    :param display_level: (int) level above which the message has to be displayed.
    :param verbosity_level: verbosity level
    :return: None
    """
    if verbosity_level >= display_level:
        verbosity("Warning\t" + message + "\n", display_level, verbosity_level)


def verbosity(message, display_level, verbosity_level=0):
    """Display the message and the time if verbosity level is superior to 1.

    :param message: (str) information about the current steps
    :param display_level: (int) level above which the message has to be displayed.
    :param verbosity_level: verbosity level
    :return: None
    """
    if verbosity_level >= display_level:
        sys.stderr.write(message + "\n")


def treat_one_request(database, ext, verbosity_level):
    """Treat one REST request on a given database ("ensemblgenomes" or "ensembl"),
    and return the result as a list or dictionary produced by applying decode on the JSON result.

    :param database:  (str) database name ("ensemblgenomes" or "ensembl").
    :param ext: URL extension (which specifies the type and parameters of the REST request.
    :param verbosity_level: verbosity level
    :return: a list or dictionary produced by applying the function json() on the result
    """

    if database == "ensembl":
        server = "http://rest.ensembl.org"

    else:
        server = "http://rest.ensemblgenomes.org"

    # Send the request to the REST server
    rest_url = server + ext + "?"

    verbosity(current_time() + "\t" + "Getting result from " + server +
              " REST URL (" + rest_url + "content-type=application/json" + ")" + "\n", 2, verbosity_level)

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
              " REST URL (" + rest_url + "content-type=application/json" + ")" + "\n", 2, verbosity_level)

    return r.json()


def open_output_handle(file_name=None, verbosity_level=0):
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


def report_command(command, fh, verbosity_level):
    """ Write the command line used for the request in a file handle.

    :param command: the command line used for the request in the outfile
    :param fh: a writable file handle
    :param verbosity_level: (int) verbosity level
    :return: None
    """
    if verbosity_level >= 1:
        fh.write(";" + str(" ".join(command)) + "\n")


def processing_time(str_start_time, str_end_time, fh, verbosity_level):
    """Display in the STDOUT the starting time, ending time and the time passed between the two.

    :param str_start_time: (datetime) the starting time
    :param str_end_time: (datetime) the ending time
    :param fh: a writable file handle
    :param verbosity_level: (int) level of verbosity from which the optional message should be displayed.
    :return: None
    """
    delta = datetime.strptime(str_end_time, '%Y-%m-%d %H:%M:%S') - datetime.strptime(str_start_time,
                                                                                     '%Y-%m-%d %H:%M:%S')
    if verbosity_level >= 1:
        fh.write("; Job started \t" + str_start_time + "\n")
        fh.write("; Job done \t" + str_end_time + "\n")
        fh.write("; Elapsed (seconds) \t" + str(delta.seconds) + "\n")


def create_headers(fh, verbosity_level, list_headers):
    """Display in the STDOUT the headers of every column

    :param fh:a writable file handle.
    :param verbosity_level: (int) level of verbosity from which the optional message should be displayed.
    :param list_headers: (list) List of the column's name
    :return: None
    """
    if verbosity_level >= 1:
        fh.write("; Column content \n")
        for i in range(len(list_headers)):
            fh.write("; " + str(i + 1) + "\t" + list_headers[i] + "\n")
        fh.write("; " + '\t'.join(list_headers) + "\n")


def parsing_file(input_file, column_nb):
    """Parse a .tab file and extract all the information from a column of this file.
    By default the function extract the second column.

    :param input_file: .tab file
    :param column_nb: (int) number of the column containing the information needed
    :return:(list) List of species name.
    """

    # Initialize list
    info_list = []

    file = open(input_file, "r")

    # Parse the given file
    while True:
        line = file.readline()
        # Stop the while at the end of the file
        if line == "":
            break
        # Ignore comment line beginning with ';' and empty line
        if line == '\n' or line.startswith(';') or line.startswith('";'):
            pass
        # Split each line at every tab into a list
        else:
            if '\t' in line:
                list_temp = line.strip("\n").split('\t')

            else:
                list_temp = [line.strip("\n")]
            # Select le 2nd column as species name
            try:
                # print(list_temp[column_nb-1])
                if list_temp[column_nb - 1] != "" and list_temp[column_nb - 1] != "\n":
                    info_list.append(list_temp[column_nb - 1])
            except IndexError:
                sys.stderr.write("IndexError: This column doesn't exist. Please use the option --column_nb to select a "
                                 "column containing the information requested.\n")
                sys.exit()

    return info_list


if __name__ == '__main__':
    print("This is a library for ensembl rest client.")
