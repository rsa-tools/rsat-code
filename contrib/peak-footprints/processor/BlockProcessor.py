from processor.Processor import Processor

from utils.Constants import Constants
from utils.RSATUtils import RSATUtils
from utils.log.Log import Log
from utils.exception.ExecutionException import ExecutionException


from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct

from common.PWM import PWM
from common.Motif import Motif

# This processor generates a PWM (position weight matrix) of the provided sequence alignements
# and searches through this matrix for conserved blocks of text (using the corresponding ratio matrix)
#
# Parameters:
#   WindowSize : the size of the window used to parse the sequences when searching conserved blocks
#   ResiduConservationLimit : the ratio limit to consider a residue is conserved
#   WindowConservationLimit : the ratio limit to consider a block is conserved
#   DesiredSpeciesList (optional): the list of species used in the multiple alignement to research conserved blocks
#   Algorithm : the chosen algorithm. Could be "OccurenceRatio", "InformationRatio" or "null". null mean that the whole
#               peak region is considered as a unique conserved region

class BlockProcessor( Processor):
    
    WINDOW_SIZE_PARAM = "WindowSize"
    RESIDU_CONSERVATION_LIMIT_PARAM = "ResiduConservationLimit"
    WINDOW_CONSERVATION_LIMIT_PARAM = "WindowConservationLimit"
    REFERENCE_SPECIES_PARAM = "ReferenceSpecies"
    DESIRED_SPECIES_LIST_PARAM = "DesiredSpeciesList"
    
    ALGORITHM_PARAM = "Algorithm"
    ALGORITHM_INFORMATION_RATIO_VALUE = "informationratio"
    ALGORITHM_OCCURENCE_RATIO_VALUE = "occurenceratio"
    ALGORITHM_NONE_VALUE = "null"
    
    # --------------------------------------------------------------------------------------
    def __init__( self):
        
        Processor.__init__( self)
        self.windowSize = 0
        self.conservationLimit = 0.0
        self.desiredSpeciesList = []
        self.algorithm = BlockProcessor.ALGORITHM_OCCURENCE_RATIO_VALUE


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
        
        return "Detection of conserved regions in MSA"


    # --------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return ( BlockProcessor.WINDOW_SIZE_PARAM,\
                BlockProcessor.RESIDU_CONSERVATION_LIMIT_PARAM,\
                BlockProcessor.WINDOW_CONSERVATION_LIMIT_PARAM,\
                BlockProcessor.REFERENCE_SPECIES_PARAM)



    # --------------------------------------------------------------------------------------
    # Execute the processor
    def execute( self, input_commstructs):

        if input_commstructs == None or len( input_commstructs) == 0:
            raise ExecutionException( "BlockProcessor.execute : No inputs")
        
        input_commstruct = input_commstructs[0]

        # retrieve the processor parameters
        self.windowSize = self.getParameterAsint( BlockProcessor.WINDOW_SIZE_PARAM)
        self.residuConservationLimit = self.getParameterAsfloat( BlockProcessor.RESIDU_CONSERVATION_LIMIT_PARAM)
        self.windowConservationLimit = self.getParameterAsfloat( BlockProcessor.WINDOW_CONSERVATION_LIMIT_PARAM)
        algo = self.getParameter( BlockProcessor.ALGORITHM_PARAM, False)
        if algo != None:
            self.algorithm = algo.lower()
        
        referenceSpecies = self.getParameter( BlockProcessor.REFERENCE_SPECIES_PARAM)
        
        desired_species_line = self.getParameter( BlockProcessor.DESIRED_SPECIES_LIST_PARAM, False)
        Log.trace( "BlockProcessor.execute : Chosen Algorithm is '" + self.algorithm + "'")
        
        self.desiredSpeciesList = []
        self.desiredSpeciesList.append( referenceSpecies)
        if desired_species_line != None:
            self.desiredSpeciesList.extend( desired_species_line.split())
        
        # Analyze the conserved region in each MSA
        # If 'None' algorithm is chosen, the entire MSA is considered as conserved
        for bed_seq in input_commstruct.bedToMA.keys():
            for alignment in input_commstruct.bedToMA[ bed_seq]:
                pwm = PWM()
                pwm.initFromAlignment( alignment, self.desiredSpeciesList)
                if self.algorithm != BlockProcessor.ALGORITHM_NONE_VALUE:
                    self.analyzeConservedBlocks( pwm, alignment)
                else:
                    new_block = Motif( 0, alignment.totalLength, "", pwm)
                    new_block.composeName( alignment.name)
                    alignment.addMotif( new_block, True)
        
        return input_commstruct


    # --------------------------------------------------------------------------------------
    # Search for conserved block in the PWM
    def analyzeConservedBlocks( self, pwm, alignment):
        
        block_lenghts = []
        
        index_start = 0
        index_end = index_start + self.windowSize
        
        finished = False
        left_limit = index_start
        while not finished and index_start >= 0 and index_end < pwm.totalLength:
            ratio = self.computeBlockRatio( index_start, index_end, pwm)
            if ratio >= self.windowConservationLimit:
                block = self.extendBlock( index_start, index_end, ratio, left_limit, pwm)
                block.composeName( alignment.name)
                alignment.addMotif( block, True)
                block_lenghts.append( block.pwm.totalLength)
                index_start = block.indexEnd
                index_end = index_start + self.windowSize
                left_limit = index_start
            else:
                index_start +=1
                index_end +=1
        
        return block_lenghts


    # --------------------------------------------------------------------------------------
    # Compute the conservation ratio of the given PWM block
    def computeBlockRatio(self, index_start, index_end, pwm):
      
        # Algorithm DirectRatio : Window ratio = ratio of position having a minimum max ratio
        if self.algorithm == BlockProcessor.ALGORITHM_OCCURENCE_RATIO_VALUE:
            window_ratio = 0
            for index in range( index_start, index_end):
                letter_max = pwm.getMostConservedResidu( index)
                if letter_max != None:
                    max_ratio =  pwm.ratioMatrix[ letter_max][index]
                    if  max_ratio >= self.residuConservationLimit:
                        window_ratio += 1
                    else:
                        # If the number of "-" is greater than the number of occurence of the most conserved letter
                        # the window is considered as not conserved
                        sum = 0
                        for letter in Constants.DNA_ALPHABET:
                            sum += pwm.matrix[letter][index]
                        if (pwm.nbSequences - sum) > pwm.matrix[ letter_max][index]:
                            return 0.0
                else:
                    return 0.0

            return window_ratio / float( index_end - index_start)
        
        elif self.algorithm == BlockProcessor.ALGORITHM_INFORMATION_RATIO_VALUE:
            window_ratio = 0
            for index in range( index_start, index_end):
                letter_max = pwm.informationMatrix[ Constants.MAX_INDEX][ index]
                max_info = pwm.informationMatrix[ letter_max][index]
                info_ratio = (max_info - pwm.informationLimits[ letter_max][0]) / float( pwm.informationLimits[ letter_max][1] - pwm.informationLimits[ letter_max][0])
                if info_ratio > self.residuConservationLimit : 
                    window_ratio += 1
                else:
                    # If the number of "-" is greater than the number of occurence of the most conserved letter
                    # the window is considered as not conserved
                    letter_max = pwm.getMostConservedResidu( index)
                    if letter_max != None:
                        sum = 0
                        for letter in Constants.DNA_ALPHABET:
                            sum += pwm.matrix[letter][index]
                        # If the number of "-" is greater than the number of occurence of the most conserved letter
                        # the window is considered as not conserved
                        if (pwm.nbSequences - sum) > pwm.matrix[ letter_max][ index]:
                            return 0.0
                    else:
                        return 0.0
                        
            return window_ratio / float( index_end - index_start)
                
        else:
            raise ExecutionException("BlockProcessor.computeBlockRatio : No known algorithm with name : " + self.algorithm)


    # --------------------------------------------------------------------------------------
    # Try to extend the block keeping the conservation ratio above the limit
    def extendBlock( self, index_start, index_end, ratio, left_limit, pwm):
        
        # try to extend on the right
        # note: since the intervals are semi-open, the value at position index_end was not
        # taken into account in the previous ratio. That's why, extending on the right means
        # adding the value at position index_end
        finished = False
        while not finished and index_end <= pwm.totalLength -1:
            if pwm.ratioMatrix[ Constants.MAX_INDEX][ index_end] >= 0:
                new_ratio = self.computeBlockRatio( index_start, index_end +1, pwm)
                if new_ratio < self.windowConservationLimit:
                    finished = True
                else:
                    index_end += 1
            else:
                finished = True
        
        # try to extend on the left but not more then left_limit
        finished = False
        while not finished and index_start >= 1 and index_start > left_limit:
            if pwm.ratioMatrix[ Constants.MAX_INDEX][ index_start -1] >= 0:
                new_ratio = self.computeBlockRatio( index_start -1,  index_end, pwm)
                if new_ratio < self.windowConservationLimit:
                    finished = True
                else:
                    index_start -= 1
            else:
                finished = True

        # build the resulting block
        
        block_pwm = pwm.getPWM( index_start, index_end)
        
        block = Motif( index_start, index_end, "", block_pwm)
        
        return block


    # --------------------------------------------------------------------------------------
    # Export the block length histogram and graph it
    def outputBlockLenghtHistogram( self, block_lengths, input_commstruct, out_path):
        
        # Output the histogram
        file_infos = RSATUtils.outputHistogram( block_lengths, 10, out_path, "conservedRegionSize", self.component.pipelineName, "", "Conserved block size", "Number of blocks", ('5', '6'), False)

        input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_SIZE_PATH] = file_infos[0]
        input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_SIZE_GRAPH_PATH] = file_infos[1]
    

# eflag: FileType = Python2
