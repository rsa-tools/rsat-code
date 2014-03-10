
from processor.Processor import Processor
from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct
from processor.io.BedSeqAlignmentStatsCommStruct import MotifStatistics

from utils.exception.ExecutionException import ExecutionException


# The aim of this Processor is to classify the detected motif according to their family and hit scores
# The Processor export an XML file presenting this classification both with information on motif (family, type, class, logo)
#
# Parameters:
#   MaxMotifNumber (optional): the maximum number of motifs to display in the result
#   MaxMotifByFamily (optional): the maximum number of motifs to display in the same family
#   MaxHypergeometricEValue (optional): the maximum value of the hypergeometric e-value the displayed motif must have
#   MaxChi2EValue (optional): the maximum value of the chi2 e-value the displayed motif must have

class ClassificationProcessor( Processor):
    
    MAX_MOTIF_NUMBER = "MaxMotifNumber"
    MAX_MOTIF_BY_FAMILY = "MaxMotifByFamily"
    MAX_HYP_EVALUE = "MaxHypergeometricEValue"
    MAX_CHI2_EVALUE = "MaxChi2EValue"
    
    
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

        max_motif_number = self.getParameterAsint( ClassificationProcessor.MAX_MOTIF_NUMBER, False)
        max_motif_by_family = self.getParameterAsint( ClassificationProcessor.MAX_MOTIF_BY_FAMILY, False)
        max_hyp_evalue = self.getParameterAsfloat( ClassificationProcessor.MAX_HYP_EVALUE, False)
        max_chi2_evalue = self.getParameterAsfloat( ClassificationProcessor.MAX_CHI2_EVALUE, False)

        # Classify the motif by score and family
        self.classifyFamiliesAndMotifs( input_commstruct, max_motif_number, max_motif_by_family, max_hyp_evalue, max_chi2_evalue)
                
        return input_commstruct



    # --------------------------------------------------------------------------------------
    # Classify the motifs by score and family
    def classifyFamiliesAndMotifs(self, input_commstruct, max_motif_number, max_motif_by_family, max_hyp_evalue, max_chi2_evalue):
        
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
        rank = 1
        count_motif = 1
        old_hyp_eval = 0
        old_chi2_eval = 0
        for motif_stats in motif_stats_list:
            motif_name = motif_stats.motifName
            motif_family = motif_stats.motifFamily
            # test the number of result already listed if required
            if max_motif_number != None and count_motif > max_motif_number:
                input_commstruct.motifStatistics.pop( motif_name)
                continue
            # test the value of the hypergeometric e-value if required
            hyp_evalue = motif_stats.getAttributeAsfloat( MotifStatistics.MOTIF_HYP_PVALUE)
            if max_hyp_evalue != None:
                if hyp_evalue == None or hyp_evalue > max_hyp_evalue:
                    input_commstruct.motifStatistics.pop( motif_name)
                    continue
            # test the value of the chi2 e-value if required
            chi2_evalue = motif_stats.getAttributeAsfloat( MotifStatistics.MOTIF_CHI2_PVALUE)
            if max_chi2_evalue != None:
                if chi2_evalue == None or chi2_evalue > max_chi2_evalue:
                    input_commstruct.motifStatistics.pop( motif_name)
                    continue
                            
            # add the motif in the ordered list of its family if the family is not full
            if not motif_family in motif_order_in_family.keys():
                motif_order_in_family[ motif_family] = []
            if max_motif_by_family != None and len( motif_order_in_family[ motif_family]) >= max_motif_by_family:
                input_commstruct.motifStatistics.pop( motif_name)
                continue
            motif_order_in_family[ motif_family].append( motif_name)
            # set the rank of the motif in its family as the index of the motif in its family ordered list
            motif_stats.setAttribute( MotifStatistics.MOTIF_RANK_IN_FAMILY, motif_order_in_family[ motif_family].index( motif_name) + 1)
            # add the family in the family ordered list
            if not motif_family in family_order:
                family_order.append( motif_family)
            # set the rank of the family as the index of the family in the ordered list
            motif_stats.setAttribute( MotifStatistics.MOTIF_FAMILY_RANK, family_order.index( motif_family) + 1)
            
            # Test if current motif is exaequo with previous motif
            if count_motif > 1:
                if old_hyp_eval != hyp_evalue or old_chi2_eval != chi2_evalue:
                    rank = rank + 1

            # Set the global rank of the motif
            motif_stats.setAttribute( MotifStatistics.MOTIF_RANK, rank)
            
            # Increase the motif counter
            count_motif = count_motif +1
            
            # Memorize the motif E-values
            old_hyp_eval = hyp_evalue
            old_chi2_eval = chi2_evalue
                

