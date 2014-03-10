
import os, shutil, commands, math
from collections import defaultdict

from manager.ProgressionManager import ProgressionManager

from processor.Processor import Processor
from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct
from processor.io.BedSeqAlignmentStatsCommStruct import MotifStatistics

from utils.Constants import Constants
from utils.RSATUtils import RSATUtils
from utils.MotifUtils import MotifUtils
from utils.log.Log import Log
from utils.exception.ExecutionException import ExecutionException
from utils.exception.ParsingException import ParsingException
from utils.FileUtils import FileUtils

#import site
#site.addsitedir( "/usr/share/pyshared")
#from scipy.stats import chi2

# This processor aims to create the histograms of statistics obtained from previous analysis
# It compute also a chi2 test on each histogram in comparison with the histogram that would be
# obtained in case of a uniform distribution of the motifs
#
# Parameters:
#   HistogramInterval : the interval length of value used to compute the histogram
#   ReferenceMotif    : the motif used as reference in the peaks


class HistogramProcessor( Processor):
    
    HISTOGRAM_INTERVAL_PARAM = "HistogramInterval"
    REFERENCE_MOTIF = "ReferenceMotif"


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
        
        return "Analysis of identified motifs scaterring"



    #------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return ( HistogramProcessor.HISTOGRAM_INTERVAL_PARAM, HistogramProcessor.REFERENCE_MOTIF)
        


    # --------------------------------------------------------------------------------------
    # Execute the processor
    def execute( self, input_commstructs):
    
        if input_commstructs == None or len( input_commstructs) == 0:
            raise ExecutionException( "HistogramProcessor.execute : No inputs")
        
        input_commstruct = input_commstructs[0]

        # Retrieve processor parameters
        histogram_interval = self.getParameterAsint( HistogramProcessor.HISTOGRAM_INTERVAL_PARAM)
        
        reference_motif = self.getParameter( HistogramProcessor.REFERENCE_MOTIF)
        
        # Store the common data at the stats param level
        input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.REFERENCE_MOTIF] = reference_motif
        
        # Compute the histogram and the chi2 test and draw the histogram graph
        self.buildHistogramsAndGraphs( input_commstruct, histogram_interval)
        
        # Compute the percentage of motif sites overlapping a reference motif site
        self.computeOverlappingScores( input_commstruct, reference_motif)
        
        return input_commstruct
        


    # --------------------------------------------------------------------------------------
    # For each motif, compute the histogram corresponding to the motif statistics 
    # using RSAT classfreq command, execute a chi2 test against uniform distribution and
    # draw the histogram graph using RSAT XYgraph command
    def buildHistogramsAndGraphs(self, input_commstruct, histogram_interval):

        # Retrieve the algorithm parameters
        RSAT_PATH = self.component.getParameter( Constants.RSAT_DIR_PARAM)

        # Compute the statistics of the motifs
        Log.info( "HistogramProcessor.buildHistogramsAndGraphs : collecting motifs statistics")
        statistics = self.computeMotifStatistics( input_commstruct, )
        
        hits_distances = statistics[ 0]
        motif_size_min = statistics[ 1]
        motif_size_max = statistics[ 2]
        hits_peakscore = statistics[ 3]
        
        #print "motif_size_max = " + str( motif_size_max)
        
        # Compute the uniform distribution probabilities
        Log.info( "HistogramProcessor.buildHistogramsAndGraphs : computing uniform distribution")
        uniform_distributions = self.computeUniformDistributions( input_commstruct, histogram_interval, motif_size_min, motif_size_max)
        
        # Build the output CommStruct
        Log.info( "HistogramProcessor.buildHistogramsAndGraphs : building histogram and graphs")
        
        # Execute the RSAT commands and computations
        try:
            # Prepare the output directories
            dir_path = os.path.join( self.component.outputDir, self.component.getComponentPrefix())
            shutil.rmtree( dir_path, True)
            FileUtils.createDirectory( dir_path, 0777)

            # Parse the motif list and execute the computations and commands for each of them 
            ProgressionManager.setTaskProgression( "Building motifs histogram",  self.component, 0.0)
            total_motif_number = len( hits_distances.keys())
            count_motif = 0
            for motif_name in hits_distances.keys():
                count_motif += 1
                motif_stats = input_commstruct.motifStatistics[motif_name]
                motif_id = motif_stats.motifID
                motif_size = motif_stats.motifSize
                hit_number = motif_stats.getAttributeAsint( MotifStatistics.MOTIF_HIT_SCORE)
                
                # Initialize the motif prefix ID
                #if motif_id != None and len( motif_id) > 0:
                #    prefix_id = "_" + motif_id
                #else:
                #    prefix_id = ""
                prefix_id = "";
                
                # save the stats to a tabbed file for classfreq command
                input_path = os.path.join( dir_path, motif_name + prefix_id + "_Distances.tab")
                self.outputMotifStatistics( hits_distances[ motif_name], input_path)

                # execute the classfreq command
                histo_path = os.path.join( dir_path, motif_name + prefix_id + "_Distances_histogram.tab")  
                
                cmd = os.path.join( RSAT_PATH, "perl-scripts/classfreq")
                cmd += " -i '" + input_path + "'"
                cmd += " -col 1"
                cmd += " -ci " + str( histogram_interval)
                cmd += " -o '" + histo_path + "'"
                
                cmd_result = commands.getstatusoutput( cmd)
                if cmd_result[0] != 0:
                    Log.log( "HistogramProcessor.buildHistogramsAndGraphs : status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
                    Log.log( "  Command output is = \n" + str( cmd_result[1]))
                    continue

                # retrieve the classfreq results from output file
                motif_distribution = self.parseClassfreqResults( histo_path)

                # compute the homogen distribution for the current motif
                null_distribution = self.computeMotifNullDistribution( uniform_distributions[ motif_size], hit_number)
                
                # Save both histograms to same file to create a common graph
                all_histo_path = os.path.join( dir_path, motif_name + prefix_id + "_Distances_histograms.tab")
                label1 = motif_name
                label2 = "Homogeneous model"
                self.outputAllHistograms( motif_distribution, label1, null_distribution, label2, histogram_interval, all_histo_path)
                motif_stats.setAttribute( MotifStatistics.MOTIF_DISTANCE_HISTOGRAM, all_histo_path)
                
                # Execute a chi2 test on the motif distribution against the motif homogen distribution
                chi2_test = RSATUtils.executeChi2Test( all_histo_path, 4, 5)
                if chi2_test != None:
                    motif_stats.setAttribute( MotifStatistics.MOTIF_CHI2, chi2_test[0])
                    motif_stats.setAttribute( MotifStatistics.MOTIF_CHI2_PVALUE, chi2_test[1])
                else:
                    motif_stats.setAttribute( MotifStatistics.MOTIF_CHI2, "0.0")
                    motif_stats.setAttribute( MotifStatistics.MOTIF_CHI2_PVALUE, "1.0")                    
                
                # Build the PNG graph corresponding to all histograms using RSAT XYGraph command
                graph_path = os.path.join( dir_path, motif_name + prefix_id + "_Distances.png")
                cmd = os.path.join( RSAT_PATH, "perl-scripts/XYgraph")
                cmd += " -i '" + all_histo_path + "'"
                cmd += " -title1 '" + self.component.pipelineName + "'" 
                cmd += " -title2 ''" 
                #cmd += " -xcol 3 -ycol 4,5"
                cmd += " -xcol 3 -ycol 4"
                cmd += " -xleg1 'Distance to peak maximum'"
                cmd += " -yleg1 'Number of motif hits'"
                cmd += " -legend -header -format png -fhisto"
                cmd += " -o '" + graph_path + "'"
                
                cmd_result = commands.getstatusoutput( cmd)
                if cmd_result[0] != 0:
                    Log.log( "HistogramProcessor.buildHistogramsAndGraphs : status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
                    Log.log( "  Command output is = \n" + str( cmd_result[1]))
                    continue
                    
                motif_stats.setAttribute( MotifStatistics.MOTIF_DISTANCE_HISTOGRAM_GRAPH, graph_path)
                
                # Build the PDF graph corresponding to all histograms using RSAT XYGraph command
                graph_path_pdf = os.path.join( dir_path, motif_name + prefix_id + "_Distances.pdf")
                cmd = os.path.join( RSAT_PATH, "perl-scripts/XYgraph")
                cmd += " -i '" + all_histo_path + "'"
                cmd += " -title1 '" + self.component.pipelineName + "'" 
                cmd += " -title2 ''" 
                #cmd += " -xcol 3 -ycol 4,5"
                cmd += " -xcol 3 -ycol 4"
                cmd += " -xleg1 'Distance to peak maximum'"
                cmd += " -yleg1 'Number of motif hits'"
                cmd += " -legend -header -format pdf -fhisto"
                cmd += " -o '" + graph_path_pdf + "'"
                
                cmd_result = commands.getstatusoutput( cmd)
                if cmd_result[0] != 0:
                    Log.log( "HistogramProcessor.buildHistogramsAndGraphs : status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
                    Log.log( "  Command output is = \n" + str( cmd_result[1]))
                    continue
                    
                motif_stats.setAttribute( MotifStatistics.MOTIF_DISTANCE_HISTOGRAM_GRAPH_PDF, graph_path_pdf)
                
                # Output the histogram of motif peak scores
                if len( hits_peakscore[ motif_name]) > 1:
                    valuable = False
                    for value in hits_peakscore[ motif_name]:
                        if value != 0:
                            valuable = True
                            break
                    if valuable:
                        score_histo_prefix = motif_name + prefix_id + "_PeakScores"
                        title1 = self.component.pipelineName
                        title2 = "Distribution of peak score for " + motif_name + prefix_id
                        legendx = "Peak Score"
                        legendy = "Number of occurence"
                        pathes = RSATUtils.outputHistogram( hits_peakscore[ motif_name], histogram_interval, dir_path, score_histo_prefix, title1, title2, legendx, legendy, None, True)
                        motif_stats.setAttribute( MotifStatistics.MOTIF_PEAK_SCORE_HISTOGRAM, pathes[0])
                        motif_stats.setAttribute( MotifStatistics.MOTIF_PEAK_SCORE_HISTOGRAM_GRAPH, pathes[1])

                # Update the progression
                if count_motif % 10 == 0:
                    ProgressionManager.setTaskProgression( "Building motifs histogram",  self.component, count_motif/float( total_motif_number))
                
        except IOError, io_exce:
            raise ExecutionException( "HistogramProcessor.buildHistogramsAndGraphs : Unable to build histogram and graph. From:\n\t---> " +str( io_exce))


    # --------------------------------------------------------------------------------------
    # Collect the position statistics of all recognized motifs
    def computeMotifStatistics(self, input_commstruct):
        
        hits_distances = {}
        hits_size = {}
        hits_peakscore = {}
        
        identified_motif_names = input_commstruct.motifStatistics.keys()
        
        motif_size_min = 10000
        motif_size_max = 0
        for bedseq in input_commstruct.bedToMA.keys():
            reference_index = bedseq.referenceIndex
            for alignment in input_commstruct.bedToMA[ bedseq]:
                #filtered_motifs = self.filterMotifs( alignment.motifs)
                filtered_motifs = alignment.motifs
                for motif in filtered_motifs:
                    motif_name = motif.name
                    if motif_name in identified_motif_names:
                        if not motif_name in hits_distances.keys():
                            hits_distances[ motif_name] = []
                            hits_peakscore[ motif_name] = []
                            hits_size[ motif_name] = 0
                        motif_fix_index_start = bedseq.indexStart + alignment.fixIndex( motif.indexStart)
                        motif_fix_index_end = bedseq.indexStart + alignment.fixIndex( motif.indexEnd)
                        motif_fix_index_center = int( (motif_fix_index_start + motif_fix_index_end) / float( 2))
                        distance = motif_fix_index_center - reference_index
                        hits_distances[ motif_name].append( distance)
                        hits_peakscore[ motif_name].append( bedseq.score)
                        motif_length = input_commstruct.motifStatistics[ motif_name].motifSize
                        if motif_length > motif_size_max:
                            motif_size_max = motif_length
                        if motif_length < motif_size_min:
                            motif_size_min = motif_length
        
        return (hits_distances, motif_size_min, motif_size_max, hits_peakscore)



    # --------------------------------------------------------------------------------------
    # Consider motifs of the same name and strand in a peak and remove motif intersection, keeping only the
    # longest motif alignment among those which intersect
    def filterMotifs(self, motifs):
        
        classified_motifs = {}
        
        for motif in motifs:
            if not motif.name in classified_motifs.keys():
                classified_motifs[ motif.name] =  []
                classified_motifs[ motif.name].append( motif)
            else:
                other_motifs = classified_motifs[ motif.name]
                conflict = False
                conflicted = []
                for other_motif in other_motifs:
                    if motif.strand == other_motif.strand and motif.indexStart < other_motif.indexEnd and motif.indexEnd > other_motif.indexStart:
                        conflict = True
                        conflicted.append( other_motif)
                
                if conflict == True:
                    to_add = True
                    for conflict_motif in conflicted:
                        if conflict_motif.indexEnd - conflict_motif.indexStart >= motif.indexEnd - motif.indexStart:
                            to_add = False
                    if to_add == True:
                        for conflict_motif in conflicted:
                            classified_motifs[ motif.name].remove( conflict_motif)
                        classified_motifs[ motif.name].append( motif)
                else:
                     classified_motifs[ motif.name].append( motif)
                        
            
        filtered_motifs = []
        for name in classified_motifs:
            filtered_motifs.extend( classified_motifs[name])
        
        return filtered_motifs
        
        


    # --------------------------------------------------------------------------------------
    # Compute the frequency distribution considering a uniform distrubtion on each sequence
    def computeUniformDistributions(self, input_commstruct, histogram_interval, motif_size_min, motif_size_max):
        
        # Initialize the probability table
        index_number_max = 0
        index_total_number = 0 
        for bedseq in input_commstruct.bedToMA.keys():
            reference_index = bedseq.referenceIndex
            indexes_left = reference_index - bedseq.indexStart
            index_total_number += indexes_left
            if indexes_left > index_number_max:
                index_number_max = indexes_left
            indexes_right = bedseq.indexEnd - reference_index
            index_total_number += indexes_right
            if indexes_right > index_number_max:
                index_number_max = indexes_right

        index_distribution = []* (motif_size_max +1)
        for site_size in range ( motif_size_max +1):
            index_distribution.append( defaultdict(lambda: 0))
    
        
        # Compute the probability values for each index
        bedseq_number = len( input_commstruct.bedToMA.keys())
        for bedseq in input_commstruct.bedToMA.keys():
            reference_index = bedseq.referenceIndex
            for site_size in range ( motif_size_min, motif_size_max + 1):
                index_start = bedseq.indexStart + int( site_size /float(2.0))
                index_end = bedseq.indexEnd - int( math.ceil( site_size /float(2.0)))
                for index in range( index_start, index_end):
                    index_distribution[ site_size][ index - reference_index] += 1 / float( index_total_number - site_size * bedseq_number)

        # Compute the probability values for each histogram class
        histo_distribution = []* ( motif_size_max + 1)
        for site_size in range ( motif_size_max + 1):
            histo_distribution.append( defaultdict(lambda: 0.0))
        
        for site_size in range ( motif_size_min, motif_size_max + 1):
            class_number = 0
            index = 0
            # Compute histogram for positive classes
            while index < index_number_max:
                for ind in range( index, min( index + histogram_interval, index_number_max)):
                    histo_distribution[ site_size][ class_number] += index_distribution[site_size][ ind]
                class_number += histogram_interval
                index += histogram_interval
            
            # Compute histogram for negative classes
            class_number = -histogram_interval
            index = -1
            while index > -index_number_max:
                for ind in range( index, max( index - histogram_interval, -index_number_max), -1):
                    histo_distribution[ site_size][ class_number] += index_distribution[site_size][ ind]
                class_number -= histogram_interval
                index -= histogram_interval            
        
        return histo_distribution



    # --------------------------------------------------------------------------------------
    # Ajust the uniform distribution to the motif hit number
    def computeMotifNullDistribution(self, uniform_distrib, hit_number):
        
        result = {}
        for classe in uniform_distrib.keys():
            result[ classe] = uniform_distrib[ classe] * hit_number
            
        return result

        
    # --------------------------------------------------------------------------------------
    # Parse the result of the RSAT classfreq command
    def parseClassfreqResults( self, file_path):
        
        result = {}
        
        total_number_colums = 9
        class_col = 1
        frequency_col = 3
        
        try:
            file = open( file_path, "r")
            for line in file:
                tokens = line.split()
                if len(tokens) == total_number_colums:
                    try:
                        result[ int( tokens[ class_col])] = int( tokens[ frequency_col])
                    except (TypeError, ValueError), exce :
                        raise ParsingException( "HistogramProcessor.parseClassfreqResults : Unable to get int value from histogram file : '" + file_path + "'. From:\n\t---> " + str( exce))
                else:
                    raise ParsingException( "HistogramProcessor.parseClassfreqResults : The histogram file is not correct formatted. Number of column is abnormal : '" + file_path)
            file.close()
        except IOError,  io_exce:
            raise ParsingException( "HistogramProcessor.parseClassfreqResults : Unable to read histogram file : '" + file_path + "'. From:\n\t---> " + str( io_exce))
        
        return result

    
    # --------------------------------------------------------------------------------------
    # output the given statistics to file
    def outputMotifStatistics( self, statistics, path):

        try:        
            file = open( path, "w")
            for number in statistics:
                file.write( str( number) + "\n")
                file.flush()
            file.close()
        except IOError, io_exce:
            raise ExecutionException( "HistogramProcessor.outputMotifStatistics : Unable to build statistics file. From:\n\t---> " +str( io_exce))


    # --------------------------------------------------------------------------------------
    # output the given histogram to tabbed file
    def outputAllHistograms(self, motif_histo, label1, null_histo, label2, histo_interval, path):
    
        try:
            file = open( path, "w")
            file.write("# x_min \t x_max \t x_mean \t " + label1 + "\t" + label2 + "\n")
            for classe in sorted( null_histo.keys()):
                file.write( str( classe) + "\t")
                file.write( str( classe + histo_interval) + "\t")
                file.write( str( classe + (histo_interval/2.0)) + "\t")
                if classe in motif_histo.keys():
                    file.write( str( motif_histo[ classe]) + "\t")
                else:
                    file.write( "0.0\t")
                file.write( str( null_histo[ classe]))
                file.write( "\n")
                file.flush()
            file.close()
        except IOError, io_exce:
            raise ExecutionException( "HistogramProcessor.outputAllHistograms : Unable to build all histograms file. From:\n\t---> " +str( io_exce))




    # --------------------------------------------------------------------------------------
    # Compute the percentage of motif site overlapping a reference motif sitre
    def computeOverlappingScores(self, input_commstruct, reference_motif):
        
        homo_location = defaultdict(lambda: defaultdict(lambda: 0))

        identified_motif_names = input_commstruct.motifStatistics.keys()

        # Compute the homo-location score
        for bedseq in input_commstruct.bedToMA.keys():
            for alignment in input_commstruct.bedToMA[ bedseq]:
                current_motifs = alignment.motifs
                current_motifs_length = len (current_motifs)
                for index_a in range( current_motifs_length -1):
                    motif = current_motifs[ index_a]
                    motif_name = motif.name
                    if motif_name in identified_motif_names:
                        # Compute the overlapping scores
                        for index_b in range( index_a + 1, current_motifs_length):
                            motif_b = current_motifs[ index_b]
                            if MotifUtils.intersect( motif, motif_b, 0.5):
                                homo_location[ motif_name][ motif_b.name] += 1
                            if MotifUtils.intersect( motif_b, motif, 0.5):
                                homo_location[ motif_b.name][ motif_name] += 1
        
        # Normalize the motif score with homo-location score
        for motif_name in identified_motif_names:
            motif_stats =  input_commstruct.motifStatistics[ motif_name]
            
            if motif_name != reference_motif:
                motif_hit_score = motif_stats.getAttributeAsint( MotifStatistics.MOTIF_HIT_SCORE)
                if homo_location[ reference_motif][ motif_name] != 0:
                    hits_normalized_score = min( homo_location[ motif_name][ reference_motif], homo_location[ reference_motif][ motif_name]) / float( motif_hit_score)
                else:
                    hits_normalized_score = 0.0
            else:
                hits_normalized_score = 1.0
                
            motif_stats.setAttribute( MotifStatistics.MOTIF_RATIO_HOMOLOCATION, hits_normalized_score)
