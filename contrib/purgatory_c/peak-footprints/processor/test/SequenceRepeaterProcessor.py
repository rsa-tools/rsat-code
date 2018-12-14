
from processor.Processor import Processor
from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct

from utils.exception.ExecutionException import ExecutionException

class SequenceRepeaterProcessor( Processor):
    
    
    # ---------------------------------------------------------------------------------------------
    def __init__( self):
        Processor.__init__( self)


    # ---------------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as input 
    # (None if no input CommStruct)
    @staticmethod
    def getInputCommStructClass():
        
        return ( BedSeqAlignmentStatsCommStruct, )


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
        
        return "Generation of MSA by sequence repetition"



    #------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return None


    
    # ---------------------------------------------------------------------------------------------
    # Execute the processor
    def execute( self, input_commstructs):

        if input_commstructs == None or len( input_commstructs) == 0:
            raise ExecutionException( "SequenceRepeaterProcessor.execute : No inputs")
        
        input_commstruct = input_commstructs[0]


        # In each MSA, repeat the sequence of the reference species in all the other species sequences
        for bedseq in input_commstruct.bedToMA.keys():
            for alignment in input_commstruct.bedToMA[ bedseq]:
                sequence = alignment.sequences[ alignment.referenceSpecies]
                for species in alignment.sequences:
                    if species != alignment.referenceSpecies:
                        alignment.sequences[ species] = []
                        alignment.sequences[ species].extend( sequence)
                alignment.totalLength = 0
                alignment.finalizeSequences()
                        
        return input_commstruct
                    
                
