
import os, shutil, math

from sets import Set

from manager.ProgressionManager import ProgressionManager

from processor.Processor import Processor
from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct
from processor.io.BedSeqAlignmentStatsCommStruct import MotifStatistics

from utils.Constants import Constants
from utils.RSATUtils import RSATUtils
from utils.exception.ExecutionException import ExecutionException
from utils.FileUtils import FileUtils

# This Processor 

class ContingencyProcessor(Processor):

    REFERENCE_MOTIF_PARAM = "ReferenceMotif"
    
    # ---------------------------------------------------------------------------------------------
    def __init__(self):
        Processor.__init__(self)



    # ---------------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as input 
    # (None if no input CommStruct)
    @staticmethod
    def getInputCommStructClass():
        
        return (BedSeqAlignmentStatsCommStruct,)



    # ---------------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as output
    # (None if no output CommStruct)
    @staticmethod
    def getOutputCommStructClass():
        
        return (BedSeqAlignmentStatsCommStruct,)


    # --------------------------------------------------------------------------------------
    # Returns a name that will be used as display name in the user friendly outputs
    @staticmethod
    def getDisplayName():
        
        return "Generation of motif binding sites contingency table"



    # --------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return (ContingencyProcessor.REFERENCE_MOTIF_PARAM,)
        

    
    # ---------------------------------------------------------------------------------------------
    # Execute the processor
    def execute(self, input_commstructs):

        if input_commstructs == None or len(input_commstructs) == 0:
            raise ExecutionException("ContingencyProcessor.execute : No inputs")
        
        input_commstruct = input_commstructs[0]

        # Retrieve the processor parameters
        reference_motif = self.getParameter(ContingencyProcessor.REFERENCE_MOTIF_PARAM)
        
        # Prepare the processor output dir
        dir_path = os.path.join(self.component.outputDir, self.component.getComponentPrefix())
        shutil.rmtree(dir_path, True)
        FileUtils.createDirectory(dir_path, 0777)
        
        # Initialize motif contingency statistics
        for motif_name in input_commstruct.motifStatistics.keys():
            motif_statistics = input_commstruct.motifStatistics[ motif_name]
            motif_statistics.setAttribute( MotifStatistics.CONTIGENCY_MOTIF_COOCCURENCE, 0)
            motif_statistics.setAttribute( MotifStatistics.CONTIGENCY_REFERENCE_MOTIF_BEDSEQ, 0)
            motif_statistics.setAttribute( MotifStatistics.CONTINGENCY_TOTAL_BEDSEQ, 0)

        # Compute the contingency table
        input_commstruct = self.computeContingencyTable( input_commstruct, reference_motif)

        return input_commstruct
        
        
    # ---------------------------------------------------------------------------------------------
    # Compute for each motif the number of BED sequence in which the
    # considered motif and the reference motif are found together
    def computeContingencyTable(self, input_commstruct, reference_motif):
        
        count_bedseq = 0
        contingency_dict = {}
        # Count the number of BED sequence common with reference motif for each other motifs
        # -Look at each chromosome
        for chrom in input_commstruct.bedSequencesDict.keys():
            # Look at each BED sequence in the current chromosome
            for bedseq in input_commstruct.bedSequencesDict[chrom]:
                count_bedseq = count_bedseq + 1
                # If the current BED sequence has alignment, look at the TFBS
                # found in the BED sequence
                if input_commstruct.bedToMA.has_key( bedseq):
                    motif_name_list = Set()
                    motif_id_list = Set()
                    for alignment in input_commstruct.bedToMA[ bedseq]:
                        for motif in alignment.motifs:
                            motif_name_list.add( motif.name)
                            if motif.id != None and len( motif.id) > 0:
                                motif_id_list.add( motif.id)
                    
                    if reference_motif in motif_name_list:
                        for motif_name in motif_name_list:
                            if contingency_dict.has_key( motif_name):
                                contingency_dict[ motif_name] = contingency_dict[ motif_name] + 1 
                            else:
                                contingency_dict[ motif_name] = 1
        
        # Assign the count to the input_commstruct
        for motif_name in contingency_dict.keys():
            if motif_name in input_commstruct.motifStatistics.keys():
                motif_statistics = input_commstruct.motifStatistics[ motif_name]
                motif_statistics.setAttribute( MotifStatistics.CONTIGENCY_MOTIF_COOCCURENCE, contingency_dict[ motif_name])
                motif_statistics.setAttribute( MotifStatistics.CONTIGENCY_REFERENCE_MOTIF_BEDSEQ, contingency_dict[ reference_motif]) 
                motif_statistics.setAttribute( MotifStatistics.CONTINGENCY_TOTAL_BEDSEQ, count_bedseq)
    
        return input_commstruct      
                    
                    
                    



