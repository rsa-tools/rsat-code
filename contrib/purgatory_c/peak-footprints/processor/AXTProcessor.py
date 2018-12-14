
import os

from processor.Processor import Processor
from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct
from processor.io.BedSeqCommStruct import BedSeqCommStruct

from common.SequenceAlignment import SequenceAlignment

from utils.FileUtils import FileUtils
from utils.Constants import Constants
from utils.log.Log import Log
from utils.exception.ExecutionException import ExecutionException
from utils.exception.ParsingException import ParsingException

# Parameters:
#   AXTFile : the path to the AXT file(s)
#   ReferenceSpecies : the species taken as reference in the AXT file(s)
#   DesiredSpeciesList (optionnal): the list of aligned species to take into account

class AXTProcessor( Processor):
    
    INPUT_AXT_FILE_PARAM = "AXTFile"
    REFERENCE_SPECIES_PARAM = "ReferenceSpecies"
    DESIRED_SPECIES_LIST_PARAM = "DesiredSpeciesList"
    
    # ---------------------------------------------------------------------------------------------
    def __init__( self):
        Processor.__init__( self)
        self.desiredSpeciesList = []
        self.referenceSpecies = None
        self.bedSequencesDict = None
        self.axtBlockDic = {}


    # ---------------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as input 
    # (None if no input CommStruct)
    @staticmethod
    def getInputCommStructClass():
        
        return ( BedSeqCommStruct, )


    # ---------------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as output
    # (None if no output CommStruct)
    @staticmethod
    def getOutputCommStructClass():
        
        return ( BedSeqAlignmentStatsCommStruct, )


    # --------------------------------------------------------------------------------------
    # Returns a name that will be used as display name in the user friendly outputs
    @staticmethod
    def getDisplayName():
        
        return "Research for matching pairwise alignements"
    
    
    # --------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return ( AXTProcessor.INPUT_AXT_FILE_PARAM, AXTProcessor.REFERENCE_SPECIES_PARAM)
        
    
    # ---------------------------------------------------------------------------------------------
    # Execute the processor
    def execute( self, input_commstructs):
        
        if input_commstructs == None or len( input_commstructs) == 0:
            raise ExecutionException( "AXTProcessor.execute : No inputs")
        
        input_commstruct = input_commstructs[0]

        # Retrieve processor parameters
        source_axtfile = self.getParameter( AXTProcessor.INPUT_AXT_FILE_PARAM)
            
        desired_species_line = self.getParameter( AXTProcessor.DESIRED_SPECIES_LIST_PARAM, False)
        if desired_species_line != None:
            self.desiredSpeciesList = desired_species_line.split()
            
        self.referenceSpecies = self.getParameter( AXTProcessor.REFERENCE_SPECIES_PARAM)        
        
        # Retrieve BED sequences from the input CommStruct
        self.bedSequencesDict = input_commstruct.bedSequencesDict
        
        if self.bedSequencesDict == None or len( self.bedSequencesDict) == 0:
            raise ExecutionException( "AXTProcessor.execute : No BEDSequence provided as input")
        
        # look for MAF files to parse
        axt_file_list = {}
        source_axtfile = os.path.join( source_axtfile,  self.referenceSpecies)
        subdir_list = FileUtils.getDirectoryList( source_axtfile)
        for dir in subdir_list:
            if len( self.desiredSpeciesList) == 0 or (dir in self.desiredSpeciesList):
                sub_list = FileUtils.getFileList( os.path.join( source_axtfile, dir), "axt")
                if sub_list != None and len( sub_list) > 0:
                    if not dir in axt_file_list.keys():
                        axt_file_list[ dir] = []
                    axt_file_list[dir] = sub_list
                else:
                    Log.log( "AXTProcessor.execute : There is no AXT file for query species : '" + dir + "' in path + '" + os.path.join( source_axtfile,  dir) + "'")
        
        if len( axt_file_list) == 0:
            raise ExecutionException( "AXTProcessor.execute : The path '" + source_axtfile + "' does not point to a directory named '" + self.referenceSpecies + "' containing query species directories with AXT files")
        
        # parse the list of MAF files
        for species in axt_file_list.keys():
            for file_path in axt_file_list[ species]:
                self.parseFile( file_path, species)
        
        # Assign the whole list of parsed species if a desired list was not set
        if len( self.desiredSpeciesList) == 0:
            self.desiredSpeciesList = axt_file_list.keys()

        # TO REMOVE
        for sequence in self.axtBlockDic.keys():
            print sequence.toString()
            for block in self.axtBlockDic[ sequence]:
                print block.toString( self.referenceSpecies)
        # END TO REMOVE

        output_commstruct = BedSeqAlignmentStatsCommStruct()
        output_commstruct.processorName = self.component.processorName
        output_commstruct.baseSpecies = input_commstruct.baseSpecies
        output_commstruct.bedSequencesDict = input_commstruct.bedSequencesDict

        for bed_sequence in self.axtBlockDic.keys():
            alignment = self.composeSequenceAlignment( bed_sequence)
            output_commstruct.addSequenceAlignment( bed_sequence, alignment)
        
        return output_commstruct


    # --------------------------------------------------------------------------------------
    # For each BEDSequence, try to recompose the corresponding DNA sequence from the
    # reference sequences in the associated AXTblocks and build at the same time the
    # multiple alignment
    def composeSequenceAlignment( self, bed_sequence):
        
        axt_blocs = self.axtBlockDic[ bed_sequence]
        for index in range( bed_sequence.getLength()):
            dna_index = bed_sequence.indexStart
            no_insertion = True
            for axt_bloc in axt_blocs:
                no_insertion = no_insertion and axt_bloc.hasOnlyResidu( dna_index)
            
            if no_insertion:
                for axt_bloc in axt_blocs:
                    pass
                    
        
        
        final_seq_align = SequenceAlignment()
        final_seq_align.name =  bed_sequence.name + "_1"
        final_seq_align.initializeWithDots( bed_sequence.getLength(), self.referenceSpecies, self.desiredSpeciesList)
        for axt_block in self.axtBlockDic[ bed_sequence]:
            ref_dna_seq = axt_block.sequences[ self.referenceSpecies]

            # Determine the start and end indexes (related to the string) of the intersection between
            # the BED sequence and the reference DNA sequence
            dna_index_start = max( bed_sequence.indexStart - ref_dna_seq.indexStart,  0)
            dna_index_end = min( bed_sequence.indexEnd - ref_dna_seq.indexStart, ref_dna_seq.textLength)
            
            # Compute the start and end index the block must be placed to in the final alignenment text
            bed_index_start = max( ref_dna_seq.indexStart - bed_sequence.indexStart,  0)
            bed_index_end = bed_index_start + (dna_index_end - dna_index_start)
            
            # Modify the start and end indexes according to the number of "-" present
            # in the sequence text to catch all the required letters
            current_index = 0
            while current_index <= dna_index_start:
                if ref_dna_seq.text[current_index] == Constants.SEQUENCE_INSERTION_CHAR:
                    dna_index_start += 1
                    dna_index_end += 1
                current_index += 1
            
            while current_index < dna_index_end:
                if ref_dna_seq.text[current_index] == Constants.SEQUENCE_INSERTION_CHAR:
                    dna_index_end += 1
                current_index += 1
            
            # Modify the bed_index_start and bed_index_end according to 
            # the number of "-" present in the sequence text
            # of the reference species currently present in the alignment
            count = 0
            for i in range( bed_index_start+1):
                if final_seq_align.sequences[ self.referenceSpecies][ i] == Constants.SEQUENCE_INSERTION_CHAR:
                    count += 1
            bed_index_start += count
            bed_index_end += count

            #print "bed_sequence.indexStart = " + str( bed_sequence.indexStart)
            #print "ref_dna_seq.indexStart = " + str( ref_dna_seq.indexStart)
            #print "dna_index_start = " + str( dna_index_start)
            #print "bed_index_start = " + str( bed_index_start)
            #print "bed_sequence.indexEnd = " + str( bed_sequence.indexEnd)
            #print "ref_dna_seq.indexStart = " + str( ref_dna_seq.indexStart)
            #print "ref_dna_seq.textLength = " + str( ref_dna_seq.textLength)
            #print "dna_index_end = " + str( dna_index_end)
            #print "bed_index_end = " + str( bed_index_end)

            # Insert the block of sequence of each species in the final alignement 
            for species in axt_block.sequences.keys():
                final_seq_align.insertSequenceBlock( species, bed_index_start, bed_index_end, axt_block.sequences[species].text[ dna_index_start : dna_index_end]   )
            
        final_seq_align.finalizeSequences()
            
        return final_seq_align



    # --------------------------------------------------------------------------------------
    # Parse the file at the given path if it is an AXT file.
    # 
    # Returns a dictionnary with the BEDSequence as key and the list of corresponding AXTBlock as value
    def parseFile(self, file_name, query_species):
        
        try:
            input_file = open( file_name, 'r')
            
            # if it is a axt file, verify if an index file exists
            indexed = False
            try:
                index_path =  file_name + "index"
                input_index_file = open( index_path,  "r")
                indexed = True
            except IOError:
                pass
            
            if indexed == True:
                Log.trace( "AXTProcessor.parseFile : parsing file '" + file_name + "' using index '" + index_path + "'")
                self.parseBlockListWithIndex( input_index_file, input_file, query_species)
                self.closeFile( input_index_file)
            else:
                Log.trace( "AXTProcessor.parseFile : parsing file '" + file_name + "'")
                self.parseBlockListWithoutIndex( input_file, query_species)
                
            self.closeFile( input_file)
            return
                
        except IOError,  io_exec:
            raise ParsingException( "AXTProcessor.parseFile : Enable to open file '" + file_name + "'. From:\n\t---> " + str(io_exec))



    # --------------------------------------------------------------------------------------
    # Parse the alignment block of the file linearly
    def parseBlockListWithoutIndex(self, input_file, query_species):
        
        # search for the next line starting with 'a' (meaning new alignment lbock)
        counter = 0
        while 1:
            line = input_file.readline()
            if len( line) == 0:
                break
            elif not line.isspace() and line[0] != Constants.COMMENT_CHAR:
                tokens = line.split()
                if tokens != None and len( tokens) ==9:
                    counter += 1
                    if counter % 100000 == 0:
                        Log.trace( "AXTIndexerProcessor.execute : Number of MSA already parsed : " + str( counter))
                    parsed = self.parseBlock( input_file, tokens, query_species)
                    if not parsed:
                        return



    # --------------------------------------------------------------------------------------
    # Parse a alignment block
    def parseBlock( self, input_file, block_def, query_species, indexed = False):
                
        new_block = None
        
        chrom_col = 1
        ref_chromosom = block_def[ chrom_col].lower()
        
        bed_sequences = self.getAssociatedBEDSequences( self.referenceSpecies + "." + ref_chromosom)
        if bed_sequences != None and len( bed_sequences) > 0:
            ref_start = self.getIntValue( block_def[ chrom_col+1])
            ref_end = self.getIntValue( block_def[ chrom_col+2])
            query_chromosom = block_def[ chrom_col+3].lower()
            query_start = self.getIntValue( block_def[ chrom_col+4])
            query_end = self.getIntValue( block_def[ chrom_col+5])
            new_block = self.findMatchingBEDSequences( bed_sequences, ref_start, ref_end)
            if new_block != None:
                line = input_file.readline()
                if line != None and not line.isspace():
                    new_block.addReferenceSequence( self.referenceSpecies, ref_chromosom, ref_start, ref_end - ref_start, line)
                else:
                    Log.log( "AXTProcessor.parseBlock : Wrongly formatted line in AXT block. Should be reference species sequence text : '" + line + "'")
                line = input_file.readline()
                if line != None and not line.isspace():
                    new_block.addQuerySequence( query_species, query_chromosom, query_start, query_end - query_start, line)
                else:
                    Log.log( "AXTProcessor.parseBlock : Wrongly formatted line in AXT block. Should be query species sequence text : '" + line + "'")
        else:
            # There is no BED sequence for this chromosom so the file has to be skipped
            return False

        return True


    # --------------------------------------------------------------------------------------
    # Return the list of BEDsequences having the same species and chromosom
    # than the given ones
    def getAssociatedBEDSequences(self, spec_chrom):
        
        if( self.bedSequencesDict.has_key( spec_chrom)):
            return self.bedSequencesDict[spec_chrom]
        else:
            return None



    # --------------------------------------------------------------------------------------
    # Detect if the given sequence intersect with some of the given BED sequence
    # if so, a new instance of AXTblock is created and attached to the suitable BED Sequence
    def findMatchingBEDSequences(self, bed_sequences, dna_start, dna_end):
        
        new_block = None

        for bed_sequence in bed_sequences:
            #print str( dna_start) + "-" + str( dna_end) + " vs " + str( bed_sequence.indexStart) + "-" +  str( bed_sequence.indexEnd)
            if dna_end > bed_sequence.indexStart and dna_start < bed_sequence.indexEnd:
                
                if new_block == None:
                    new_block = AXTBlock()
                        
                if not self.axtBlockDic.has_key( bed_sequence):
                    self.axtBlockDic[ bed_sequence] = []
                    
                self.axtBlockDic[ bed_sequence].append( new_block)
        
        return new_block



    # --------------------------------------------------------------------------------------
    # Close the given File
    def closeFile(self, file):
        
        try:
            file.close()
        except IOError,  exce:
            Log.log( "AXTProcessor.closeFile : Enable to close file '" + file + "'. From:\n\t--> " +  str(exce))



    # --------------------------------------------------------------------------------------
    # Return the int value of the given token
    def getIntValue(self, token):
        
        try:
            return int( token)
        except ValueError,  val_exce:
            raise ParsingException( "AXTProcessor : Unable to get integer value of '" + token + "'. From:\n\t---> " + str( val_exce))



# ############################################################

#               AXTBlock Class

# ############################################################

from common.DNASequence import DNASequence

class AXTBlock:
    
    # --------------------------------------------------------------------------------------
    def __init__(self):
        
        self.referenceSpecies = ""
        self.referenceSequence = ""
        self.querySpecies = ""
        self.querySequence = ""
        self.inverseSequences = False
        
        self.parsingIndex = 0
        self.parsingStarted = False
        self.onlyResidu = True


    # --------------------------------------------------------------------------------------
    def addReferenceSequence(self, species, chromosom, bp_start, text_lenght,  text):
        
        if text_lenght != len( text):
            print "ATXBlock.addReferenceSequence : incohence between text length = " + str( text_lenght) + " and length of text string = " + str( len(text))
        
        if self.inverseSequences:
            text = self.inverseText( text)
        sequence = DNASequence( species, chromosom, bp_start, text_lenght, text)
        self.referenceSpecies = species
        self.referenceSequence = sequence
        self.parsingIndex = bp_start



    # --------------------------------------------------------------------------------------
    def addQuerySequence(self, species, chromosom, bp_start, text_lenght,  text):
        
        if text_lenght != len( text):
            print "ATXBlock.addReferenceSequence : incohence between text length = " + str( text_lenght) + " and length of text string = " + str( len(text))
        
        if self.inverseSequences:
            text = self.inverseText( text)
        sequence = DNASequence( species, chromosom, bp_start, text_lenght, text)
        self.querySpecies = species
        self.querySequence = sequence



    # --------------------------------------------------------------------------------------
    # Inverse the order of the characters in the given text
    def inverseText(self,  text):
        
        result = ""
        length = len( text)
        for index in range( length):
            result += text[ length - index -1]
            
        return result



    # --------------------------------------------------------------------------------------
    def hasOnlyResidu(self, index):
        
            if self.parsingStarted == False and index == self.parsingIndex:
                self.parsingStarted = True
            
            if self.parsingStarted == True:
                if self.referenceSequence[ self.parsingIndex] == Constants.SEQUENCE_INSERTION_CHAR or self.querySequence[ self.parsingIndex] == Constants.SEQUENCE_INSERTION_CHAR:
                    self.onlyResidu = False
                    return False
            
            self.onlyResidu = True
            return True
            
            
    
    # --------------------------------------------------------------------------------------
    def toString(self,  reference_species=""):
        
        result = "AXT Block : \n"
        if reference_species == "":
            for species in sorted( self.sequences.keys()):
                sequence = self.sequences[ species]
                result += "|--" + sequence.toString() + "\n"
        else:
            sequence = self.sequences[ reference_species]
            result += "|--" + sequence.toString() + "\n\t\t" + sequence.text + "\n"
            for species in sorted( self.sequences.keys()):
                if species != reference_species:
                    sequence = self.sequences[ species]
                    result += "|--" + sequence.toString() + "\n"
        return result

# eflag: FileType = Python2

 
