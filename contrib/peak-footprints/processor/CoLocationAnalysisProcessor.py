
import os, shutil, math

from manager.ProgressionManager import ProgressionManager

from processor.Processor import Processor
from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct
from processor.io.BedSeqAlignmentStatsCommStruct import MotifStatistics

from utils.RSATUtils import RSATUtils
from utils.exception.ExecutionException import ExecutionException
from utils.FileUtils import FileUtils

# This Processor computes the distance between each motif hit and the nearest site of the reference motif. Those distances
# are reported in histograms indenpendently for each identified motif.
#
# Parameters:
#   ReferenceMotif = Motif that can be used as reference (should be the motif of the TF studied by the chIP-seq)
#   HistogramInterval = Class size of the computed histograms
#   MaximalDistance = The distance between the reference motif and the considered motif above which the hit is not reported

class CoLocationAnalysisProcessor(Processor):
    
    REFERENCE_MOTIF_PARAM = "ReferenceMotif"
    HISTOGRAM_INTERVAL_PARAM = "HistogramInterval"
    MAXIMAL_DISTANCE_PARAM = "MaximalDistance"
    
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
        
        return "Analysis of identified motifs co-locations"



    # --------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return (CoLocationAnalysisProcessor.REFERENCE_MOTIF_PARAM, CoLocationAnalysisProcessor.HISTOGRAM_INTERVAL_PARAM, CoLocationAnalysisProcessor.MAXIMAL_DISTANCE_PARAM)
        

    
    # ---------------------------------------------------------------------------------------------
    # Execute the processor
    def execute(self, input_commstructs):

        if input_commstructs == None or len(input_commstructs) == 0:
            raise ExecutionException("CoLocationAnalysisProcessor.execute : No inputs")
        
        input_commstruct = input_commstructs[0]

        # Retrieve the processor parameters
        reference_motif = self.getParameter(CoLocationAnalysisProcessor.REFERENCE_MOTIF_PARAM)
        
        histogram_interval = self.getParameterAsint(CoLocationAnalysisProcessor.HISTOGRAM_INTERVAL_PARAM)
        
        distance_max = self.getParameterAsint(CoLocationAnalysisProcessor.MAXIMAL_DISTANCE_PARAM)
        
        # compute the distance of motifs to the given reference motif
        final_distances = self.computeMotifDistances(reference_motif, input_commstruct, distance_max)

        # Prepare the processor output dir
        dir_path = os.path.join(self.component.outputDir, self.component.getComponentPrefix())
        shutil.rmtree(dir_path, True)
        FileUtils.createDirectory(dir_path, 0777)
        
        # Output the distance distribution histogram of each identified motif
        ProgressionManager.setTaskProgression("Building histograms", self.component, 0.0)
        total_motif = len(final_distances.keys())
        count_motif = 0
        for motif_name in final_distances.keys():
            if motif_name in input_commstruct.motifStatistics.keys():
                motif_stats = input_commstruct.motifStatistics[ motif_name]
                if len(final_distances[ motif_name]) > 1:
                    prefix = motif_name + "_CoLocation"
                    title1 = self.component.pipelineName
                    title2 = "Distribution of distances between " + motif_name + " and " + reference_motif
                    legendx = "Distance to " + reference_motif
                    legendy = "Number of occurence"
                    result_info = RSATUtils.outputHistogram(final_distances[ motif_name], histogram_interval, dir_path, prefix, title1, title2, legendx, legendy, None, True)
                    motif_stats.setAttribute(MotifStatistics.MOTIF_COLOCATION_HISTOGRAM, result_info[ 0])
                    motif_stats.setAttribute(MotifStatistics.MOTIF_COLOCATION_HISTOGRAM_GRAPH, result_info[ 1])
            count_motif += 1
            ProgressionManager.setTaskProgression("Building histograms", self.component, count_motif / float(total_motif))
            
        ProgressionManager.setTaskProgression("Building histograms", self.component, 1.0)

        return input_commstruct
        
        

    # ---------------------------------------------------------------------------------------------
    # Compute the distance of motifs to the given reference motif
    def computeMotifDistances(self, reference_motif, input_commstruct, distance_max):
        
        final_distances = {}
        
        ProgressionManager.setTaskProgression("Computing distances", self.component, 0.0)
        
        total_chrom = len(input_commstruct.bedSequencesDict.keys())
        count_chrom = 0
        for chromosom in input_commstruct.bedSequencesDict.keys():
            # print "---------------------------------------------"
            # print "Chromosom " + chromosom
            current_ref_motif = None
            current_studied_motifs = []
            motifs_middle_positions = {}
            bed_sequence_list = input_commstruct.bedSequencesDict[ chromosom]
            bed_sequence_list.sort(self.bedSequenceComparator)
            bedtoma_keys = input_commstruct.bedToMA.keys()
            for bed_sequence in bed_sequence_list:
                # print bed_sequence.toString()
                if bed_sequence in bedtoma_keys:
                    motif_list = []
                    for msa in input_commstruct.bedToMA[ bed_sequence]:
                        for motif in msa.motifs:
                            motifs_middle_positions[ motif] = self.getMiddlePosition(motif, bed_sequence.indexStart + msa.fixIndex(motif.indexStart))
                            motif_list.append(motif)
                    motif_list.sort(self.motifComparator)
                    for motif in motif_list:
                        # print motif.toString()
                        # print " ===> middle position = " + str( motifs_middle_positions[ motif])
                        if motif.name != reference_motif:
                            current_studied_motifs.append(motif)
                        else:
                            # print " ---> Reference Motif"
                            distances = self.computeLocalDistances(current_studied_motifs, motifs_middle_positions, current_ref_motif, motif, distance_max)
                            for motif_name in distances.keys():
                                if not motif_name in final_distances.keys():
                                    final_distances[ motif_name] = []
                                final_distances[ motif_name].extend(distances[ motif_name])
                            current_studied_motifs = []
                            current_ref_motif = motif
            count_chrom += 1
            ProgressionManager.setTaskProgression("Computing distances", self.component, count_chrom / float(total_chrom))
        
        ProgressionManager.setTaskProgression("Computing distances", self.component, 1.0)
        
        return final_distances
                        


    # ---------------------------------------------------------------------------------------------
    # Compute the distances between the motifs of the given list and the nearest reference motif site
    def computeLocalDistances(self, current_studied_motifs, motifs_middle_positions, current_ref_motif, new_ref_motif, distance_max):
    
        distances = {}
        
        for motif in current_studied_motifs:
            motif_name = motif.name
            motif_length = motif.indexEnd - motif.indexStart;
            motif_middle_position = motifs_middle_positions[ motif]
            
            distance = 0
            # If a reference motif site has been located previously, compare its distance
            # to the one with to the new motif site
            if current_ref_motif != None:
                distance1 = motif_middle_position              
                distance1 -= motifs_middle_positions[ current_ref_motif]

                distance2 = motif_middle_position
                distance2 -= motifs_middle_positions[ new_ref_motif]
                
                if math.fabs(distance1) > math.fabs(distance2):
                    ref_motif_length = new_ref_motif.indexEnd - new_ref_motif.indexStart
                    distance = distance2
                else:
                    ref_motif_length = current_ref_motif.indexEnd - current_ref_motif.indexStart
                    distance = distance1
            # If no reference motif site was found before, uses the distance to the new reference motif site
            else:
                ref_motif_length = new_ref_motif.indexEnd - new_ref_motif.indexStart
                distance2 = motif_middle_position
                distance2 -= motifs_middle_positions[ new_ref_motif]
                distance = distance2
            
            # Consider only the distance under the distance max
            if distance >= -distance_max and distance <= distance_max:
                # #(commented)Consider only the distances that show motifs are not overlapping
                # if math.fabs(distance) > max( motif_length / float(2), (ref_motif_length / float(2))) :
                if not motif_name in distances.keys():
                    distances[ motif_name] = []
                distances[ motif_name].append(distance)
            
            # Remove studied motif from the list of motif to analyse
            del motifs_middle_positions[ motif]
        
        return distances


    # ---------------------------------------------------------------------------------------------
    # Compute the position index of the given motif that start at the given position
    def getMiddlePosition(self, motif, start_pos):
        
        return start_pos + int((motif.indexEnd - motif.indexStart) / float(2))



    # ---------------------------------------------------------------------------------------------
    # Compare two BED Sequence classifying them according to their start index
    def bedSequenceComparator(self, bed_seq1, bed_seq2):
        
        if bed_seq1.species != bed_seq2.species:
            raise ExecutionException("CoLocationAnalysisProcessor.bedSequenceComparator : Unable to compare two BED sequences of differents species : '" + bed_seq1.species + "' != '" + bed_seq2.species + "'") 
            
        if bed_seq1.chromosom != bed_seq2.chromosom:
            raise ExecutionException("CoLocationAnalysisProcessor.bedSequenceComparator : Unable to compare two BED sequences of differents chromosom : '" + bed_seq1.chromosom + "' != '" + bed_seq2.chromosom + "'") 
        
        return bed_seq1.indexStart - bed_seq2.indexStart


    # ---------------------------------------------------------------------------------------------
    # Compare two motifs classifying them according to their start index
    def motifComparator(self, motif1, motif2):
        
        return motif1.indexStart - motif2.indexStart
