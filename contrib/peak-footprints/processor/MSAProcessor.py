from utils.FileUtils import FileUtils

import os, commands,  shutil

from processor.Processor import Processor
from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct

from utils.log.Log import Log
from utils.exception.ExecutionException import ExecutionException
from utils.Constants import Constants

from common.SequenceAlignment import SequenceAlignment

# This processor aims to align mutliple sequences using standard algorithm like Clus
# Parameters:
#   Method : the method used to obtain the result. Value should be "ClustalW" or "MAFFT"
#   For method ClustalW:
#       DesiredSpeciesList (optional): the species to take into account in the multiple alignments
#       CommandOptions (optional): the extra options to set in the ClustalW command line
#   For method MAFFT:
#       DesiredSpeciesList (optional): the species to take into account in the multiple alignments
#       CommandOptions (optional): the extra options to set in the MAFFT command line

class MSAProcessor(Processor):
    
    METHOD_PARAM = "Method"
    METHOD_VALUE_CLUSTALW = "clustalw"
    METHOD_VALUE_MAFFT = "mafft"
    
    # Parameters for ClustalW
    DESIRED_SPECIES_LIST_PARAM = "DesiredSpeciesList"
    COMMAND_OPTIONS_PARAM = "CommandOptions"
    
    
    # --------------------------------------------------------------------------------------
    def __init__( self):
        Processor.__init__( self)


    # --------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as input 
    # (None if no input CommStruct)
    @staticmethod
    def getInputCommStructClass():
        
        return ( BedSeqAlignmentStatsCommStruct, )


    # --------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as output
    # (None if no output CommStruct)
    @staticmethod
    def getOutputCommStructClass():
        
        return ( BedSeqAlignmentStatsCommStruct, )



    # --------------------------------------------------------------------------------------
    # Returns a name that will be used as display name in the user friendly outputs
    @staticmethod
    def getDisplayName():
        
        return "Alignment of multiple sequences"



    #------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return ( MSAProcessor.METHOD_PARAM, )



    # --------------------------------------------------------------------------------------
    # Execute the processor
    def execute( self, input_commstructs):
        
        if input_commstructs == None or len( input_commstructs) == 0:
            raise ExecutionException( "MAFProcessor.execute : No inputs")
        
        input_commstruct = input_commstructs[0]
        
        # Retrieve the Processor parameters
        method = self.getParameter( MSAProcessor.METHOD_PARAM).lower()
        
        #Select the method to use
        if method == MSAProcessor.METHOD_VALUE_CLUSTALW:
            self.executeClustalW( input_commstruct)
        elif method == MSAProcessor.METHOD_VALUE_MAFFT:
            self.executeMAFFT( input_commstruct)
        else:
            raise ExecutionException( "MSAProcessor.execute : required method is not known : " + method)


        return input_commstruct
        


    # ---------------------------------------------------------------------------------------------
    # Execute ClustalW to generate the MSA
    def executeClustalW( self, input_commstruct):
        
        #Retrieve the method parameters
        
        desired_species_line = self.getParameter( MSAProcessor.DESIRED_SPECIES_LIST_PARAM, False)
        if desired_species_line != None:
            desired_species_list = desired_species_line.split()
        else:
            desired_species_list = []
        
        command_options_line = self.getParameter( MSAProcessor.COMMAND_OPTIONS_PARAM, False)
        if command_options_line == None:
            command_options = ""
        else:
            command_options = command_options_line
        
        # Prepare the outputdir for FASTA file export
        file_info = self.prepareOutputDir()
        dir_path = file_info[0]
        file_name = file_info[1]
        file_path = os.path.join( dir_path, file_name + ".fasta")
        

        # Change directory to output dir
        working_dir = os.getcwd()
        os.chdir( dir_path)
        
        command = self.component.getParameter( Constants.CLUSTALW_COMMAND_PARAM)
         
        # Compose the ClustalW command line with all required options
        output_filepath = file_path + "result.txt"
        cmd = command
        cmd += " -INFILE=" + file_path 
        cmd += " -ALIGN"
        cmd += " -TYPE=DNA"
        cmd += " -OUTFILE=" + output_filepath
        cmd += " " + command_options
        
        for bed_sequence in input_commstruct.bedToMA.keys():
            final_result = []
            for alignment in input_commstruct.bedToMA[ bed_sequence]:
                #output the alignment to FASTA file
                self.outputAlignmentToFASTAFile( alignment, file_path, desired_species_list)
                # Execute the command
                cmd_result = commands.getstatusoutput( cmd)
                if cmd_result[0] != 0:
                    Log.log( "MSAProcessor.executeClustalW : status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
                    Log.log( "MSAProcessor.executeClustalW : command output is = \n" + str( cmd_result[1]))
                    continue
                
                # Parse the result of the compare-matrices command to get the result list
                final_result.append( self.parseClustalWResult( output_filepath, desired_species_list))
            if final_result != None:
                input_commstruct.bedToMA[ bed_sequence] = final_result
        
        
        # Change dir to previous working dir
        os.chdir( working_dir)
        
        
        
    # ---------------------------------------------------------------------------------------------
    # Prepare the output directory where FASTA files will be exported
    def prepareOutputDir(self):
        
        try:
            dir_path = os.path.join( self.component.outputDir, self.component.getComponentPrefix())
            shutil.rmtree( dir_path, True)
            FileUtils.createDirectory( dir_path, 0777)
            file_name = "motif"
            file_path = os.path.join( dir_path, file_name)
            return (dir_path, file_path)
        except IOError, io_exce:
            raise ExecutionException( "MSAProcessor.prepareOutputDir : Unable to create output directory for FASTA file export : '" + dir_path + "'. From \n\t" +str( io_exce))



    # ---------------------------------------------------------------------------------------------
    # Output the list of motif to a file in MEME format
    def outputAlignmentToFASTAFile(self, alignment, file_path, desired_species_list):

        try:
            file =open( file_path, "w")
            file.write( alignment.convertToFASTA( desired_species_list))
            file.close()
        except IOError, io_exce:
            raise ExecutionException( "MSAProcessor.outputAlignmentToFASTAFile : Unable to save alignment to FASAT file : '" + file_path + "'. From \n\t" +str( io_exce))



    # ---------------------------------------------------------------------------------------------
    # Parse the output of ClustalW command
    def parseClustalWResult( self, file_path, desired_species_list):
        
        try:
            result = {}
            length = 0
            file = open( file_path, "r")
            for line in file:
                tokens=line.split()
                if tokens != None and len( tokens) == 2:
                    species = tokens[0]
                    if desired_species_list == None or len( desired_species_list) == 0 or species in desired_species_list:
                        if not species in result.keys():
                            result[species] = []
                        result[species].extend( tuple( tokens[1]))
                        length = len( result[species])
            
            result = self.removeFirstAndLastNoInfoColumns( result, length)
                    
            alignment = SequenceAlignment()
            for species in result:
                alignment.addSequence( species, result[species])
            file.close()
            return alignment
            
        except IOError, io_exce:
            raise ExecutionException( "MSAProcessor.parseClustalWResult : Unable to open the ClustalW result file : '" + file_path + "'. From:\n\t---> " +str( io_exce))



    # ---------------------------------------------------------------------------------------------
    # Execute MAFFT to generate the MSA
    def executeMAFFT( self, input_commstruct):
        
        #Retrieve the method parameters
        
        desired_species_line = self.getParameter( MSAProcessor.DESIRED_SPECIES_LIST_PARAM, False)
        if desired_species_line != None:
            desired_species_list = desired_species_line.split()
        else:
            desired_species_list = []
        
        command_options_line = self.getParameter( MSAProcessor.COMMAND_OPTIONS_PARAM, False)
        if command_options_line == None:
            command_options = ""
        else:
            command_options = command_options_line
        
        # Prepare the outputdir for FASTA file export
        file_info = self.prepareOutputDir()
        dir_path = file_info[0]
        file_name = file_info[1]
        file_path = os.path.join( dir_path, file_name + ".fasta")
        
        # Change directory to output dir
        working_dir = os.getcwd()
        os.chdir( dir_path)
        
        # Compose the MAFFT command line with all required options
        output_filepath = file_path + "result.txt"
        command = self.component.getParameter( Constants.MAFFT_COMMAND_PARAM, False)
        if command == None or len( command) == 0:
            command = "mafft --auto"
        if command_options == None or len( command_options) == 0:
            cmd = command + " --anysymbol"
            cmd += " " + file_path 
            cmd += " > "+ output_filepath
        else:
            cmd = command + " --anysymbol"
            cmd += " " + command_options + " " + file_path 
            cmd += " > " + output_filepath

        for bed_sequence in input_commstruct.bedToMA.keys():
            final_result = []
            for alignment in input_commstruct.bedToMA[ bed_sequence]:
                #output the alignment to FASTA file
                self.outputAlignmentToFASTAFile( alignment, file_path, desired_species_list)
                # Execute the command
                cmd_result = commands.getstatusoutput( cmd)
                if cmd_result[0] != 0:
                    Log.log( "MSAProcessor.executeMAFFT : status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
                    Log.log( "MSAProcessor.executeMAFFT : command output is = \n" + str( cmd_result[1]))
                    continue
                
                # Parse the result of the compare-matrices command to get the result list
                final_result.append( self.parseMAFFTResult( output_filepath, desired_species_list))
                
            if final_result != None:
                input_commstruct.bedToMA[ bed_sequence] = final_result
        
        # Change dir to previous working dir
        os.chdir( working_dir)


    # ---------------------------------------------------------------------------------------------
    # Parse the output of ClustalW command
    def parseMAFFTResult( self, file_path, desired_species_list):
        
        try:
            result = {}
            length = 0
            file = open( file_path, "r")
            line = file.readline()
            while True:
                if len( line) == 0:
                    break
                else:
                    tokens=line.split()
                    if tokens != None and len( tokens) == 2 and tokens[0] == ">":
                        species = tokens[1]
                        if desired_species_list == None or len( desired_species_list) == 0 or species in desired_species_list:
                            if not species in result.keys():
                                result[species] = []
                            while True:
                                line = file.readline()
                                if len( line) == 0:
                                    break
                                else:
                                    tokens = line.split()
                                    if tokens != None and len( tokens) > 0:
                                        if tokens[0] == ">":
                                            break
                                        else:
                                            result[species].extend( tuple( tokens[0]))
                                            length = len( result[species])
                    else:
                        line = file.readline()

            result = self.removeFirstAndLastNoInfoColumns( result, length)

            alignment = SequenceAlignment()
            for species in result:
                alignment.addSequence( species, result[species])
                
            file.close()
            
            return alignment
            
        except IOError, io_exce:
            raise ExecutionException( "MSAProcessor.parseMAFFTResult : Unable to open the MAFFT result file : '" + file_path + "'. From:\n\t---> " +str( io_exce))



    # ---------------------------------------------------------------------------------------------
    # Remove alignment only dashes (or dot) columns from the firsts and lasts positions
    def removeFirstAndLastNoInfoColumns(self, align, length ):

        for index in range( length):
            dash_residu = True
            for species in align:
                if align[ species][ index] != Constants.SEQUENCE_INSERTION_CHAR and align[ species][ index] != Constants.SEQUENCE_INIT_CHAR:
                    dash_residu = False
            if dash_residu == False:
                for species in align:
                    align[ species] = align[ species][index:]
                break
                
        for index in range( length-1):
            dash_residu = True
            for species in align:
                if align[ species][ -index-1] != Constants.SEQUENCE_INSERTION_CHAR and align[ species][ -index-1] != Constants.SEQUENCE_INIT_CHAR:
                    dash_residu = False
            if dash_residu == False:
                if index > 0:
                    for species in align:
                        align[ species] = align[ species][:-index]
                break
                    
        return align
