
import os

from common.PWM import PWM

from utils.RSATUtils import RSATUtils
from utils.log.Log import Log
from utils.exception.ExecutionException import ExecutionException
from utils.exception.ParsingException import ParsingException


class MotifUtils:

    JASPAR_FLAT_DB_PATH = ""

    # ---------------------------------------------------------------------------------------------
    # Retrieve the size of the motifs in given Transfac database file
    @staticmethod
    def getMotifsSizesFromTransfac( database_file_path = None):

        sizes = {}

        if database_file_path == None:
            database_file_path = RSATUtils.RSAT_JASPAR_MOTIF_DATABASE

        try:
            database_file = open( database_file_path, "r")
            for line in database_file:
                # detect the transfac definition starting line
                if line[0:2] == "AC":
                    tokens = line.split()
                    motif_name = tokens[1]
                    # get the definition until the definition final line ("//")
                    for line in database_file:
                        if line[0:2] == "PO":
                            # read the values of the PWM and count the lines
                            size = 0
                            for line in database_file:
                                if line[0:2] != "XX":
                                    size += 1
                                else:
                                    break
                            break
                        elif line[0:2] == "//":
                            break
                    # assign the size to the corresponding motif name
                    if size != 0:
                        sizes[ motif_name] = size
        except IOError, io_exce:
            raise ExecutionException( "MotifUtils.getMotifsSizeFromTransfacDefinition : Unable to read motif definition from database file '" + database_file_path + "'. From:\n\t---> " + str( io_exce))

        return sizes


    # ---------------------------------------------------------------------------------------------
    # Retrieve the number of the motifs in given Transfac database file
    @staticmethod
    def getMotifsNumberFromTransfac( database_file_path = None):

        motif_count = 0

        if database_file_path == None:
            database_file_path = RSATUtils.RSAT_JASPAR_MOTIF_DATABASE

        try:
            database_file = open( database_file_path, "r")
            for line in database_file:
                # detect the transfac definition starting line
                if line[0:2] == "AC":
                    motif_count = motif_count + 1
        except IOError, io_exce:
            raise ExecutionException( "MotifUtils.getMotifsNumberFromTransfac : Unable to read motif definition from database file '" + database_file_path + "'. From:\n\t---> " + str( io_exce))

        return sizes


    # ---------------------------------------------------------------------------------------------
    # Retrieve the motif details (id, family, type, class) of the motifs in given Transfac database file
    @staticmethod
    def getMotifIDFromTransfac( searched_motif_name, database_file_path = None):
        
        if database_file_path == None:
            database_file_path = RSATUtils.RSAT_JASPAR_MOTIF_DATABASE

        try:
            database_file = open( database_file_path, "r")
            for line in database_file:
                # detect the transfac definition starting line
                if line[ 0:2] == "AC":
                    tokens = line.split()
                    motif_name = tokens[1]
                    if( searched_motif_name == motif_name):
                        # get the ID
                        for line in database_file:
                            if line[ 0:2] == "ID":
                                sub_tokens = line.split()
                                return sub_tokens[ 1]
                            if line[0:2] == "//":
                                break
        except IOError, io_exce:
            raise ExecutionException( "MotifUtils.getMotifIDFromTransfac : Unable to read motif definition from database file '" + database_file_path + "'. From:\n\t---> " + str( io_exce))

        return None


    # ---------------------------------------------------------------------------------------------
    # Retrieve the motif details (id, family, type, class) of the motifs in given Transfac database file
    @staticmethod
    def getMotifsDetailsFromTransfac( database_file_path = None):

        id = {}
        family = {}
        type = {}
        classe = {}
        
        if database_file_path == None:
            database_file_path = RSATUtils.RSAT_JASPAR_MOTIF_DATABASE

        try:
            database_file = open( database_file_path, "r")
            for line in database_file:
                # detect the transfac definition starting line
                if line[ 0:2] == "AC":
                    tokens = line.split()
                    motif_name = tokens[1]
                    # get the definition until the definition final line ("//")
                    for line in database_file:
                        if line[ 0:2] == "ID":
                            sub_tokens = line.split()
                            id[ motif_name] = sub_tokens[ 1]
                        if line[ 0:2] == "CC":
                            sub_tokens = line.split()
                            if sub_tokens[ 1].lower() == "family:":
                                family[ motif_name] = sub_tokens[ 2]
                            elif sub_tokens[ 1].lower() == "type:":
                                type[ motif_name] = sub_tokens[ 2]
                            elif sub_tokens[ 1].lower() == "class:":
                                classe[ motif_name] = sub_tokens[ 2]
                        elif line[0:2] == "//":
                            break
        except IOError, io_exce:
            raise ExecutionException( "MotifUtils.getMotifsDetailsFromTransfac : Unable to read motif definition from database file '" + database_file_path + "'. From:\n\t---> " + str( io_exce))

        return ( id, family, type, classe)


    # ---------------------------------------------------------------------------------------------
    # Retrieve the motif details (id, family, type, class) of the motifs in given Transfac database file
    @staticmethod
    def getMotifsDetailsFromJaspar():

        matrix_path = os.path.join( MotifUtils.JASPAR_FLAT_DB_PATH, "MATRIX.txt")
        matrix_annotation_path = os.path.join( MotifUtils.JASPAR_FLAT_DB_PATH, "MATRIX_ANNOTATION.txt")

        names = {}
        id = {}
        family = {}
        type = {}
        classe = {}
        
        try:
            matrix_file = open( matrix_path, "r")
            matrix_annotation_file = open( matrix_annotation_path,  "r")
            
            for line in matrix_file:
                tokens = line.split()
                if len( tokens) >= 5:
                    current_num = tokens[ 0]
                    current_name = tokens[ 2] + "." + tokens[3]
                    current_id = "".join( tokens[ 4:])
                    names[ current_num] = current_name
                    id[ current_name] = current_id
                else:
                    raise ParsingException( "MotifUtils.getMotifsDetailsFromJaspar : Matrix file is not correctly formatted: 5 columns required while " + str( len( tokens)) + " columns are found")
            
            for line in matrix_annotation_file:
                tokens = line.split()
                current_num = tokens[ 0]
                if current_num in names.keys():
                    current_key = tokens[ 1]
                    current_value = "".join( tokens[2:])
                    if current_key == "family":
                        family[ names[ current_num]] = current_value
                    elif current_key == "class":
                        classe[ names[ current_num]] = current_value
                    elif  current_key == "type":
                        type[ names[ current_num]] = current_value
                else:
                    Log.log( "MotifUtils.getMotifsDetailsFromJaspar : Motif number was not detected in matrix file : " + current_num)
            matrix_annotation_file.close()
            matrix_file.close()
        except (IOError, ParsingException),  exce:
            Log.log( "MotifUtils.getMotifsDetailsFromJaspar : unable to read motifs definition. From:\n\t---> " + str( exce))

        return ( id, family, type, classe)


    # ---------------------------------------------------------------------------------------------
    # Return the Transfac definition of each given motif
    @staticmethod
    def getMotifsDefinitionFromTF( source_motif_name_list, database_file_pathes = None):
        
        motif_name_list = []
        motif_name_list.extend( source_motif_name_list)
        
        if database_file_pathes == None:
            database_file_pathes = []
            database_file_pathes.append( RSATUtils.RSAT_JASPAR_MOTIF_DATABASE)

        motif_to_definition = {}
            
        try:
            while len( database_file_pathes) > 0 and len( motif_name_list) > 0:
                database_file_path = database_file_pathes[ 0]
                database_file = open( database_file_path, "r")
                for line in database_file:
                    # detect the transfac definition starting line
                    if line[0:2] == "AC":
                        tokens = line.split()
                        # test if one of the motif correspond to the current definition
                        motif_name = tokens[1]
                        if motif_name in motif_name_list:
                            motif_to_definition[ motif_name] = []
                            motif_to_definition[ motif_name].append(line)
                            # get the definition until the definition final line ("//")
                            for line in database_file:
                                if line[0:2] == "//":
                                    break
                                else:
                                    motif_to_definition[ motif_name].append(line)
                            motif_name_list.remove( motif_name)
                    
                database_file.close()
                database_file_pathes.remove( database_file_path)
        except IOError, io_exce:
            raise ExecutionException( "MotifUtils.getMotifsDefinitionFromTF : Unable to read motif definition from database file '" + database_file_path + "'. From:\n\t---> " + str( io_exce))

        return motif_to_definition

    # ---------------------------------------------------------------------------------------------
    # Retrieve the PWM of given motifs in given Transfac database file
    @staticmethod
    def getMotifsPWMFromJasparTF( motif_list, database_file_path):
        
        name_to_motif = {}
        for motif in motif_list:
            if not motif.name in name_to_motif.keys():
                name_to_motif[ motif.name] = []
            name_to_motif[ motif.name].append( motif)

        try:
            database_file = open( database_file_path, "r")
            motif_found = 0
            for line in database_file:
                # detect the transfac definition starting line
                if line[0:2] == "AC":
                    tokens = line.split()
                    # test if one of the motif correspond to the current definition
                    current_pwm = None
                    for motif_name in name_to_motif.keys():
                        if motif_name == tokens[1]:
                            # get the definition until the definition final line ("//")
                            for line in database_file:
                                if line[0:2] == "PO":
                                    # read the values of the PWM
                                    current_pwm = PWM()
                                    for line in database_file:
                                        if line[0:2] != "XX":
                                            tokens = line.split()
                                            current_pwm.addColumn( (int( tokens[1]), int( tokens[2]),int( tokens[3]), int( tokens[4])))
                                        else:
                                            break
                                    break
                            # assign the PWM to the corresponding motifs
                            for motif in name_to_motif[ motif_name]:
                                motif.pwm = current_pwm
                                motif_found += 1
                            break
                # Test if all the searched motifs have been found
                if motif_found >= len( motif_list):
                    break
            database_file.close()
        except IOError, io_exce:
            raise ExecutionException( "MotifUtils.getMotifsPWMFromJasparTF : Unable to read motif definition from database file '" + database_file_path + "'. From:\n\t---> " + str( io_exce))



    # ---------------------------------------------------------------------------------------------
    # Retrieve the PWM of given motifs in given MEME database file
    @staticmethod
    def getMotifsPWMFromMeme( motif_list, database_file_path):
        
        name_to_motif = {}
        for motif in motif_list:
            if not motif.name in name_to_motif.keys():
                name_to_motif[motif.name] = []
            name_to_motif[ motif.name].append(motif)

        #try:
        database_file = open( database_file_path, "r")
        motif_found = 0
        for line in database_file:
            # detect the transfac definition starting line
            if line[0:5] == "MOTIF":
                tokens = line.split()
                # test if one of the motif correspond to the current definition
                current_pwm = None
                for motif_name in name_to_motif.keys():
                    if motif_name == tokens[1]:
                        # get the definition until the definition final line ("//")
                        for line in database_file:
                            if line[0:25] == "letter-probability matrix":
                                infos = line.split()
                                nsites = 0
                                for index in range( len( infos)):
                                    if infos[index] == "nsites=":
                                        nsites = int( infos[ index+1])
                                # read the values of the PWM
                                current_pwm = PWM()
                                for line in database_file:
                                    if not line.isspace():
                                        tokens = line.split()
                                        current_pwm.addColumn( (int( float( tokens[0])*nsites), int( float( tokens[1])*nsites), int( float( tokens[2])*nsites), int( float( tokens[3])*nsites)))
                                    else:
                                        break
                                break
                        # assign the PWM to the corresponding motifs
                        for motif in name_to_motif[ motif_name]:
                            motif.pwm = current_pwm
                            motif_found += 1
                        break
            # Test if all the searched motifs have been found
            if motif_found >= len( motif_list):
                break
        #except IOError, io_exce:
        #    raise ExecutionException( "MotifUtils.getMotifsPWMFromMeme : Unable to read motif definition from database file '" + database_file_path + "'. From:\n\t---> " + str( io_exce))



    # ---------------------------------------------------------------------------------------------
    # Detect if the two given motif have a non  null intersection
    @staticmethod
    def intersect( motif1, motif2, ratio):
        
        if motif1.indexStart < motif2.indexEnd and motif1.indexEnd > motif2.indexStart:
            intersect_start = max( motif1.indexStart, motif2.indexStart)
            intersect_end = min( motif1.indexEnd, motif2.indexEnd)
            if (intersect_end - intersect_start) / float( motif1.indexEnd - motif2.indexStart) >= ratio:
                return True
            
        return False
