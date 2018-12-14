
import os, shutil, random, commands

from manager.ProgressionManager import ProgressionManager

from utils.RSATUtils import RSATUtils
from utils.Constants import Constants
from utils.log.Log import Log
from utils.exception.ExecutionException import ExecutionException

from processor.Processor import Processor
from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct
from processor.io.BedSeqCommStruct import BedSeqCommStruct

from common.SequenceAlignment import SequenceAlignment
from common.BEDSequence import BEDSequence

# The purpose of this Processor is to create controled random dataset in order to bu used as test data
# of the complete workflow
#
# Parameters:
#   PeakNumber : the number of BED Sequence to create
#   PeakMediumSize : the medium size of peaks (sizes will be randomly chosen between this value +/-30%)
#   MSASize : the number of species represented in the MSA
#   MotifList : the list of motif to be inserted (spaces separated list)
#   OptimizeMotif : If set to True, the site generated from the motif definition will correspond to the
#                   max occurence residu at each position. If set to False, random sites will be generated
#                   from the motif PWM.
#   DatabaseFilePath : The full path to the Transfac file where motif PSSM are stored
#   SiteNumber : The total number of binding site to insert in the sequences
#   DistributionMode : The mode used for the sites distribution. Should be "Uniform" or "Centered"


class GenerateMSAProcessor(Processor):
    
    PEAK_NUMBER_PARAM = "PeakNumber"
    PEAK_MEDIUM_SIZE_PARAM = "PeakMediumSize"
    MSA_SIZE_PARAM = "MSASize"
    TRIVIAL_SEQUENCES_PARAM = "TrivialSequences"
    INSERTION_NUMBER_PARAM = "InsertionNumber"

    
    
    # ---------------------------------------------------------------------------------------------
    def __init__( self):
        Processor.__init__( self)


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
        
        return "Generation of random MSA"
        
 

    #------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return ( GenerateMSAProcessor.PEAK_NUMBER_PARAM, GenerateMSAProcessor.INSERTION_NUMBER_PARAM, GenerateMSAProcessor.MSA_SIZE_PARAM)



    
    # ---------------------------------------------------------------------------------------------
    # Execute the processor
    def execute( self, input_commstructs):

        if input_commstructs == None or len( input_commstructs) == 0:
            raise ExecutionException( "ImplantSitesProcessor.execute : No inputs")
        
        input_commstruct = input_commstructs[0]

        # Retrieve the processor parameters
        bedseq_number = self.getParameterAsint( GenerateMSAProcessor.PEAK_NUMBER_PARAM)
        
        insertion_number = self.getParameterAsint( GenerateMSAProcessor.INSERTION_NUMBER_PARAM)
        
        bedseq_medium_length = self.getParameterAsint( GenerateMSAProcessor.PEAK_MEDIUM_SIZE_PARAM, False)
        
        msa_length = self.getParameterAsint( GenerateMSAProcessor.MSA_SIZE_PARAM)
        
        trivial_msa = self.getParameter( GenerateMSAProcessor.TRIVIAL_SEQUENCES_PARAM, False)
        if trivial_msa == None:
            trivial_msa = False
        else:
            trivial_msa = (trivial_msa.lower() == "true")

        # Prepare the processor output dir
        out_path = os.path.join( self.component.outputDir, self.component.getComponentPrefix())
        shutil.rmtree( out_path, True)
        os.mkdir( out_path)

        # Build the output CommStruct
        output_commstruct = BedSeqAlignmentStatsCommStruct()
        output_commstruct.baseSpecies = input_commstruct.baseSpecies

        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.REFERENCE_SPECIES] = output_commstruct.baseSpecies
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.ALIGNED_SPECIES] = ""
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.BEDSEQUENCE_NUMBER] = bedseq_number
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.BEDSEQUENCE_WITH_MSA_NUMBER] = bedseq_number

        if bedseq_medium_length == None:
            # Get the required number of sequence from BED sequence list
            count_peak = 0
            max_length = 0
            for chrom in input_commstruct.bedSequencesDict:
                output_commstruct.bedSequencesDict[chrom] = []
                for bedseq in input_commstruct.bedSequencesDict[ chrom]:
                    output_commstruct.bedSequencesDict[chrom].append( bedseq)
                    length = bedseq.indexEnd - bedseq.indexStart
                    if length > max_length:
                        max_length = length
                    count_peak += 1
                    if count_peak >= bedseq_number:
                        break
                if count_peak >= bedseq_number:
                    break
        else:
            # Generate random BED sequences
            self.generateBedSequences( bedseq_number, bedseq_medium_length, output_commstruct)
        
        # Export the new bedsequence size histogram and graph
        self.outputSequenceSizeHistogram( output_commstruct)        
        
        # Generate MSA for each BED Sequence
        Log.trace( "GenerateMSAProcessor.execute : Generating MSA")
        ProgressionManager.setTaskProgression( "Generating MSA", self.component, 0.0)
        if trivial_msa:
            self.generateTrivialMSA( msa_length, bedseq_number, output_commstruct)
        else:
            self.generateRandomMSA( msa_length, bedseq_number, bedseq_medium_length, output_commstruct)
        ProgressionManager.setTaskProgression( "Generating MSA", self.component, 1.0)
        
        # Implant insertion characters into the MSA Sequences
        Log.trace( "GenerateMSAProcessor.execute : Implanting insertions")
        ProgressionManager.setTaskProgression( "Implanting insertions", self.component, 0.0)
        self.implantInsertions( output_commstruct, insertion_number)
        ProgressionManager.setTaskProgression( "Implanting insertions", self.component, 1.0)
        
        # Export the new bedsequence size histogram and graph
        self.outputMSALenghtHistogram( output_commstruct)

        
        return output_commstruct


    # ---------------------------------------------------------------------------------------------
    # Generate a dictionnary of random BED sequences ordered by species and chromosom
    def generateBedSequences( self, bedseq_number, max_length, output_commstruct):
        
        chrom = "chr1"
        output_commstruct.bedSequencesDict[chrom] = []
        
        start = 0
        random.seed()
        for count in range( bedseq_number):
            jump = int( random.uniform( 50, 150))
            start = start + jump
            length = int( random.uniform( max_length*0.7, max_length*1.3))
            end = start + length
            output_commstruct.bedSequencesDict[chrom].append( BEDSequence( output_commstruct.baseSpecies, chrom, start, end))
            start = end


    # ---------------------------------------------------------------------------------------------
    # Build a MSA composed of repeated random DNA sequences generated by RSAT random-seq
    def generateRandomMSA( self, msa_length, bedseq_number, max_length, output_commstruct):
        
        # Retrieve method required parameters
        RSAT_PATH = self.component.getParameter( Constants.RSAT_DIR_PARAM)
        dir_path = os.path.join( self.component.outputDir, self.component.getComponentPrefix())
        file_path = os.path.join( dir_path, "random_sequences.txt")
        
        try:
            # Execute the RSAT random-seq command
            cmd = os.path.join( RSAT_PATH , "perl-scripts/random-seq")
            cmd += " -l " + str( int( max_length * 1.5))
            cmd += " -n " + str( bedseq_number)
            cmd += " -a a:t 0.3 c:g 0.2"
            cmd += " -type DNA"
            cmd += " -format multi"
            cmd += " -o " + file_path
        
            Log.info( "GenerateMSAProcessor.generateMSA : starting random sequence generation. Command used is : " + cmd)
            
            # Execute the command
            cmd_result = commands.getstatusoutput( cmd)
            if cmd_result[0] != 0:
                Log.log( "GenerateMSAProcessor.generateMSA : status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
                Log.log( "GenerateMSAProcessor.generateMSA : command output is = \n" + str( cmd_result[1]))
                raise ExecutionException( "GenerateMSAProcessor.generateMSA : Cannot execute random-seq commands. See logs for more details")
            
            # Read the output file to get the random sequences
            sequence_list = []
            sequence_file = open( file_path, "r")
            for line in sequence_file:
                sequence_list.append( line.split()[0])
            
            # Generate  the species list
            species_list = []
            species_list.append( output_commstruct.baseSpecies)
            for index in range( msa_length-1):
                species_list.append( "Species" + str( index +1))
            
            # Create and fill the MSA for each BED sequence
            count_seq = 0
            for chrom in output_commstruct.bedSequencesDict.keys():
                for bedseq in output_commstruct.bedSequencesDict[ chrom]:
                    msa = SequenceAlignment()
                    msa.name =  bedseq.name + "_1"
                    msa.referenceSpecies = output_commstruct.baseSpecies
                    seq_length = bedseq.indexEnd - bedseq.indexStart
                    sequence = list( sequence_list[count_seq][: seq_length])
                    for index in range( msa_length):
                         msa.addSequence( species_list[index], sequence)
                         #msa.addSequence( species_list[index], list(['.'] * len( sequence)))
                    msa.finalizeSequences()
                    output_commstruct.addSequenceAlignment( bedseq, msa)
                    count_seq += 1

        except IOError, io_exce:
            raise ExecutionException( "GenerateMSAProcessor.generateMSA : Unable to save/read random sequences file. From:\n\t---> " + str( io_exce))


    # ---------------------------------------------------------------------------------------------
    # Build a MSA composed of repeated trivial DNA sequences (sequences with only dots)
    def generateTrivialMSA( self, msa_length, bedseq_number, output_commstruct):
    
        # Generate  the species list
        species_list = []
        species_list.append( output_commstruct.baseSpecies)
        for index in range( msa_length-1):
            species_list.append( "Species" + str( index +1))
        
        # Create and fill the MSA for each BED sequence
        for chrom in output_commstruct.bedSequencesDict.keys():
            for bedseq in output_commstruct.bedSequencesDict[ chrom]:
                msa = SequenceAlignment()
                msa.name =  bedseq.name + "_1"
                seq_length = bedseq.indexEnd - bedseq.indexStart
                sequence = list(['.'] * seq_length)
                for index in range( msa_length):
                     msa.addSequence( species_list[index], sequence)
                msa.finalizeSequences()
                output_commstruct.addSequenceAlignment( bedseq, msa)


    # ---------------------------------------------------------------------------------------------
    # Add insertion characters to the MSA in order to create break in the conservation sequences
    def implantInsertions(self, output_commstruct, insertion_number):

        bedseq_list = output_commstruct.bedToMA.keys()
        bedseq_list_length = len( bedseq_list)
        
        count_insertion = 0
        while count_insertion < insertion_number:
            chosen_bedseq = bedseq_list[ int( random.uniform( 0, bedseq_list_length))]
            msa_numbers = len( output_commstruct.bedToMA[ chosen_bedseq])
            chosen_msa = output_commstruct.bedToMA[ chosen_bedseq][ int( random.uniform( 0, msa_numbers))]
            chosen_index = random.randint( 0, chosen_msa.totalLength-1)
            for sequence in chosen_msa.sequences.values():
                sequence[chosen_index] = Constants.SEQUENCE_INSERTION_CHAR
            count_insertion += 1


    # --------------------------------------------------------------------------------------
    # Compute and output to file and graph the sequence size histogram
    def outputSequenceSizeHistogram(self, output_comm_struct):
        
        sizes = []
        
        for chrom in output_comm_struct.bedSequencesDict.keys():
            for sequence in output_comm_struct.bedSequencesDict[chrom]:
                sizes.append( sequence.indexEnd - sequence.indexStart)

        dir_path = os.path.join( self.component.outputDir, self.component.getComponentPrefix())
        
        # Output the histogram
        file_infos = RSATUtils.outputHistogram( sizes, 10, dir_path, "sequenceSize", self.component.pipelineName, "",  "Peak size",  "Number of peaks", ('5', '6'))

        output_comm_struct.paramStatistics[ BedSeqCommStruct.BED_SEQUENCES_SIZE_PATH] = file_infos[0]
        output_comm_struct.paramStatistics[ BedSeqCommStruct.BED_SEQUENCES_SIZE_GRAPH_PATH] = file_infos[1]



    # --------------------------------------------------------------------------------------
    # Export the MSA lengths histogram and graph it
    def outputMSALenghtHistogram( self, output_commstruct):
        
        msa_lenghts = []
        for bedseq in output_commstruct.bedToMA.keys():
            for msa in output_commstruct.bedToMA[ bedseq]:
                msa_lenghts.append( msa.totalLength)
        
        dir_path = os.path.join( self.component.outputDir, self.component.getComponentPrefix())
        
        # Output the histogram
        file_infos = RSATUtils.outputHistogram( msa_lenghts, 10, dir_path, "MSASize", self.component.pipelineName, "", "Conserved region size", "Number of regions", ('5', '6'))

        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.MSA_SIZE_PATH] = file_infos[0]
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.MSA_SIZE_GRAPH_PATH] = file_infos[1]
