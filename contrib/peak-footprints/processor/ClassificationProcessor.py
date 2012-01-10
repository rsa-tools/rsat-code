
from processor.Processor import Processor
from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct
from processor.io.BedSeqAlignmentStatsCommStruct import MotifStatistics

from utils.exception.ExecutionException import ExecutionException


# The aim of this Processor is to classify the detected motif according to their family and hit scores
# The Processor export an XML file presenting this classification both with information on motif (family, type, class, logo)
#
# Parameters:


class ClassificationProcessor( Processor):
    
    
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
        
        return "Classification of identified motifs"


    # --------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return None


    # --------------------------------------------------------------------------------------
    # Execute the processor
    def execute( self, input_commstructs):
    
        if input_commstructs == None or len( input_commstructs) == 0:
            raise ExecutionException( "ClassificationProcessor.execute : No inputs")
        
        input_commstruct = input_commstructs[0]

        # Classify the motif by score and family
        self.classifyFamiliesAndMotifs( input_commstruct)
                
        return input_commstruct



    # --------------------------------------------------------------------------------------
    # Classify the motifs by score and family
    def classifyFamiliesAndMotifs(self, input_commstruct):
        
        identified_motif_names = input_commstruct.motifStatistics.keys()
        
         # Classify the motif by families and scores and normalize the motif score with homo-location score
        motif_stats_list = []
        for motif_name in identified_motif_names:
            motif_stats =  input_commstruct.motifStatistics[ motif_name]
            if motif_stats != None:
                motif_stats_list.append( motif_stats)
        
        # Order the motif according to choices made in MotifStatistics.compare() method
        motif_stats_list.sort( MotifStatistics.compare)
        
        # set the family and motif ranks
        family_order = []
        motif_order_in_family = {}
        for motif_stats in motif_stats_list:
            motif_stats.setAttribute( MotifStatistics.MOTIF_RANK, motif_stats_list.index( motif_stats) + 1)
            motif_family = motif_stats.motifFamily
            motif_name = motif_stats.motifName
            # add the family in the family ordered list
            if not motif_family in family_order:
                family_order.append( motif_family)
            # set the rank of the family as the index of the family in the ordered list
            motif_stats.setAttribute( MotifStatistics.MOTIF_FAMILY_RANK, family_order.index( motif_family) + 1)
            # add the motif in the ordered list of its family
            if not motif_family in motif_order_in_family.keys():
                motif_order_in_family[ motif_family] = []
            motif_order_in_family[ motif_family].append( motif_name)
            # set the rank of the motif in its family as the index of the motif in its family ordered list
            motif_stats.setAttribute( MotifStatistics.MOTIF_RANK_IN_FAMILY, motif_order_in_family[ motif_family].index( motif_name) + 1)

