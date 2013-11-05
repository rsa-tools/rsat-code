
from processor.Processor import Processor

from manager.ProgressionManager import ProgressionManager

from utils.FileUtils import FileUtils
from utils.SequenceUtils import SequenceUtils
from utils.log.Log import Log
from utils.Constants import Constants
from utils.exception.ExecutionException import ExecutionException
from utils.exception.ParsingException import ParsingException

class MAFIndexerProcessor( Processor):
    
    INPUT_MAF_FILE_PARAM = "MAFFile"
    REFERENCE_SPECIES_PARAM = "ReferenceSpecies"
    
    _lineType_col = 0
    _speciesChrom_col = 1
    _startindex_col = 2
    _textlength_col = 3
    _strand_col = 4
    _source_size_col = 5
    _text_col = 6


    # --------------------------------------------------------------------------------------
    def __init__(self):
        Processor.__init__( self)
        self.referenceSpecies = ""


    # --------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as input 
    # (None if no input CommStruct)
    @staticmethod
    def getInputCommStructClass():
        
        return None


    # --------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as output
    # (None if no output CommStruct)
    @staticmethod
    def getOutputCommStructClass():
        
        return None



    # --------------------------------------------------------------------------------------
    # Returns a name that will be used as display name in the user friendly outputs
    @staticmethod
    def getDisplayName():
        
        return "Indexation of MAF files"
        


    #------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return ( MAFIndexerProcessor.INPUT_MAF_FILE_PARAM, MAFIndexerProcessor.REFERENCE_SPECIES_PARAM)
        


    # --------------------------------------------------------------------------------------
    # Execute the processor
    def execute( self, input_commstructs):
        
        source_maffile = self.getParameter( MAFIndexerProcessor.INPUT_MAF_FILE_PARAM)
        self.referenceSpecies = self.getParameter( MAFIndexerProcessor.REFERENCE_SPECIES_PARAM)
        
        # look for MAF files to parse
        maf_file_list = FileUtils.getFileList( source_maffile, "maf", self.referenceSpecies)
        if maf_file_list == None:
            raise ExecutionException( "MAFIndexerProcessor.execute : The path '" + source_maffile + "' does not point to a MAF file or a directory containing MAF files and does not contain a subdirectory '" + self.referenceSpecies + "' containing MAF files.")
        
        count_file = 0
        for maf_file_path in maf_file_list:
            Log.trace( "MAFIndexerProcessor.execute : Indexing " + maf_file_path)
            self.parseFile( maf_file_path)
            count_file += 1
            ProgressionManager.setComponentProgression( self.component, count_file/float( len( maf_file_list)))


    # --------------------------------------------------------------------------------------
    # Parse the MAF file
    def parseFile(self, maf_file_path):

        try:
            input_file = open( maf_file_path, "r")
        except IOError,  io_exec:
            raise ParsingException( "MAFIndexerProcessor.parseFile : Unable to open file '" + maf_file_path + "'. From:\n\t---> " + str(io_exec))
         
        if input_file != None:
            
            # Verify if the token '##maf' indicating a MAF file is found in the first lines
            is_maf_file = False
            while 1:
                line = input_file.readline()
                if len( line) == 0:
                    break
                elif not line.isspace():
                    tokens = line.split()
                    if tokens != None and len( tokens) > 0 and tokens[0] == "##maf":
                        is_maf_file = True
                        break
            
            if is_maf_file == True:
                output = []
                try:

                    # search for the next line starting with 'a' (meaning new alignment lbock)
                    counter = 0
                    specialized = True
                    ordered = True
                    previous_indexing = None
                    while 1:
                        line = input_file.readline()
                        if len( line) == 0:
                            break
                        elif not line.isspace():
                            tokens = line.split()
                            if tokens != None and len( tokens) > 0 and tokens[ MAFIndexerProcessor._lineType_col] == "a":
                                counter += 1
                                if counter % 100000 == 0:
                                    Log.trace( "MAFIndexerProcessor.execute : Number of MSA already indexed : " + str( counter))
                                line_number = input_file.tell()
                                indexing = self.indexBlock( input_file, previous_indexing)
                                if indexing != None:
                                    if previous_indexing != None:
                                        specialized = specialized and (indexing[1] == previous_indexing[1])
                                        ordered = ordered and (indexing[2] >= previous_indexing[3])
                                    output.append( indexing[0] + "\t" + str( line_number))
                                    previous_indexing = indexing
                    
                    #Write the result of indexing in file
                    output_path = maf_file_path + "index"
                    output_file = open( output_path, "w")
                    output_file.write( Constants.COMMENT_CHAR)
                    if specialized == True:
                        output_file.write( "\t" + previous_indexing[ 1])
                        if ordered == True:
                            output_file.write( "\t" + Constants.ORDERED)
                    else:
                        output_file.write( "\t" + Constants.MIXED)
                    output_file.write( "\n")
                    for indexing in output:
                        output_file.write( indexing + "\n")
                    output_file.flush()
                    self.closeFile( input_file)
                    self.closeFile( output_file)
                    return
                except IOError, io_exec:
                    raise ParsingException( "MAFIndexerProcessor.parseFile : Enable to create/write file '" + output_path + "'. From:\n\t---> " + str( io_exec))
                
                
            else:
                self.closeFile( input_file)
                raise ParsingException( "MAFIndexerProcessor.parseFile : The file '" + maf_file_path + "' is not a MAF file")
                
        else:
            raise ParsingException( "MAFIndexerProcessor.parseFile : Opened file is null '" + maf_file_path + "'")


    # --------------------------------------------------------------------------------------
    # Retrieve the indexation of the MSA block 
    def indexBlock(self, input_file, previous_indexing):

        indexing = []
        while 1:
            line = input_file.readline()
            # Check if the line is not void
            if len( line) == 0:
                return None
            elif not line.isspace():
                tokens = line.split()
                # Check if the line contains enough tokens
                if tokens != None and len( tokens) > MAFIndexerProcessor._text_col:
                    if tokens[ MAFIndexerProcessor._lineType_col] == 's':
                        # Verify if current sequence species match with reference species
                        spec_chrom_token =  tokens[ MAFIndexerProcessor._speciesChrom_col]
                        spec_chrom = SequenceUtils.getSpeciesAndChrom( spec_chrom_token)
                        species = spec_chrom[0]
                        if species == self.referenceSpecies:
                            strand = tokens[ MAFIndexerProcessor._strand_col]
                            start = self.computeStartIndex( tokens, strand)
                            #start = self.getIntValue( tokens[ MAFIndexerProcessor._startindex_col])
                            end = start + self.getIntValue( tokens[ MAFIndexerProcessor._textlength_col])
                            indexing.append( tokens[ MAFIndexerProcessor._speciesChrom_col] + "\t" + str( start) + "\t" + str( end))
                            indexing.append( spec_chrom_token)
                            indexing.append( start)
                            indexing.append( end)
                            return indexing
                        else:
                            return None
            else:
                # This block is void and don't have to be indexed
                return None


    # --------------------------------------------------------------------------------------
    # Compute the start index of the sequence according to the strand
    # if the strand is "-", the coordinates must be inversed
    def computeStartIndex(self, tokens,  strand):
        
        # If current sequence direction is inserved, coordinates must be transformed
        if strand != Constants.POSITIVE_STRAND:
            Log.trace( "MAFIndexer : Negative strand detected")
            source_size = self.getIntValue( tokens[ MAFIndexerProcessor._source_size_col])
            text_length = self.getIntValue( tokens[ MAFIndexerProcessor._textlength_col])
            rev_start = self.getIntValue( tokens[ MAFIndexerProcessor._startindex_col])
            bp_start = source_size + 1 - (text_length + rev_start)
        else:
            bp_start = self.getIntValue( tokens[ MAFIndexerProcessor._startindex_col])

        return bp_start


    # --------------------------------------------------------------------------------------
    # Return the token as an int value if possible
    def getIntValue(self, token):
        
        try:
            return int( token)
        except ValueError,  val_exce:
            raise ParsingException( "MAFIndexerProcessor : Unable to get integer value of '" + token + "'. From:\n\t---> " + str( val_exce))


    # --------------------------------------------------------------------------------------
    # Close the given File
    def closeFile(self, file):
        
        try:
            file.close()
        except IOError,  exce:
            Log.log( "MAFIndexerProcessor.closeFile : Enable to close file '" + file + "'. From:\n\t---> " +  str(exce))
