
import os, shutil

from processor.Processor import Processor
from processor.io.BedSeqCommStruct import BedSeqCommStruct

from utils.log.Log import Log
from utils.parser.BEDParser import BEDParser
from utils.RSATUtils import RSATUtils
from utils.exception.ParsingException import ParsingException

# This processor takes as input a BED format file and returns a BedSeqAlignmentStatsCommStruct
# containing the definition of the BEDSequences retrieved in the file
#
# Parameters:
#   BEDFile : the path to the BED file
#   Species : the species of the described sequences
#   PeakFile (optional): the path to the file giving information on the maximum of the peak 
#                        for each sequence in the BED file. Used in case of BED file represent ChIP-seq peaks result.
class BEDProcessor( Processor):
    
    INPUT_BED_FILE_PARAM = "BEDFile"
    INPUT_PEAK_FILE = "PeakFile"
    SPECIES_PARAM = "Species"
    
    # --------------------------------------------------------------------------------------
    def __init__(self):
        Processor.__init__( self)
    
    
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
        
        return ( BedSeqCommStruct,  )
        

    # --------------------------------------------------------------------------------------
    # Returns a name that will be used as display name in the user friendly outputs
    @staticmethod
    def getDisplayName():
        
        return "Parsing of BED format input file"



    # --------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return ( BEDProcessor.INPUT_BED_FILE_PARAM, BEDProcessor.SPECIES_PARAM)
        
        

    # --------------------------------------------------------------------------------------
    # Execute the processor
    def execute(self, input_comm_structs):
        
        # Retrieve the processor parameters
        bed_filepath = self.getParameter( BEDProcessor.INPUT_BED_FILE_PARAM)
        species = self.getParameter( BEDProcessor.SPECIES_PARAM)
        peak_filepath = self.getParameter( BEDProcessor.INPUT_PEAK_FILE, False)
        
        # Parse the BED file and get the BED sequences ordered by species and chromosom
        bedseq_dictionnary = BEDParser.getBEDSequenceDictionnary( species, bed_filepath)
        
        # Parse the peak info file if exists
        if peak_filepath != None and len( peak_filepath) > 0:
            chrom_col = 0
            max_peak_col = 10
            id_col = 13
            
            try:
                input_file = open( peak_filepath)
                for line in input_file:
                    tokens = line.split()
                    if len( tokens) > id_col:
                        chrom = tokens[ chrom_col]
                        max_peak = self.getTokenAsint( tokens[ max_peak_col])
                        id = tokens[ id_col]
                        if chrom != None and max_peak != None and id != None:
                            for bed_seq in bedseq_dictionnary[ species + "." + chrom]:
                                if bed_seq.id == id:
                                    bed_seq.referenceIndex = max_peak
            except ParsingException, par_exce:
                raise ParsingException( "BEDProcessor.execute : An error occured while parsing peak information file : '" + peak_filepath + "'. From:\n\t---> " + str( par_exce))
            except IOError,  io_exce:
                raise ParsingException( "BEDProcessor.execute : Unable to open peak information file : '" + peak_filepath + "'. From:\n\t---> " + str( io_exce))
        
        # Generate the output CommStruct
        output_commstruct = BedSeqCommStruct()
        output_commstruct.processorName = self.component.processorName
        output_commstruct.baseSpecies = species
        output_commstruct.bedSequencesDict = bedseq_dictionnary
        
        # Count the number of BED sequences
        bedseq_number = 0
        for chrom in output_commstruct.bedSequencesDict.keys():
            bedseq_number += len( output_commstruct.bedSequencesDict[chrom])
        
        output_commstruct.paramStatistics[ BedSeqCommStruct.BEDSEQUENCE_NUMBER] = bedseq_number
        Log.trace( "BEDProcessor.execute : Total number of BED Sequences = " + str( bedseq_number))
        
        # output the sequences size histogram
        self.outputSequenceSizeHistogram( bedseq_dictionnary, output_commstruct)
        
        return output_commstruct



    # --------------------------------------------------------------------------------------
    # Compute and output to file and graph the sequence size histogram
    def outputSequenceSizeHistogram(self, bedseq_dictionnary, output_comm_struct):
        
        sizes = []
        
        for chrom in bedseq_dictionnary.keys():
            for sequence in bedseq_dictionnary[chrom]:
                sizes.append( sequence.indexEnd - sequence.indexStart)
        
        # Prepare output directory
        dir_path = os.path.join( self.component.outputDir, self.component.getComponentPrefix())
        shutil.rmtree( dir_path,  True)
        os.mkdir( dir_path)
        
        # Output the histogram
        file_infos = RSATUtils.outputHistogram( sizes, 10, dir_path, "sequenceSize", self.component.pipelineName, "",  "Peak Size",  "Number of peaks", ('5', '6'))

        output_comm_struct.paramStatistics[ BedSeqCommStruct.BED_SEQUENCES_SIZE_PATH] = file_infos[0]
        output_comm_struct.paramStatistics[ BedSeqCommStruct.BED_SEQUENCES_SIZE_GRAPH_PATH] = file_infos[1]

