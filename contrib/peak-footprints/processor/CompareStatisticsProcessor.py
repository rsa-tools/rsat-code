
from processor.Processor import Processor
from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct

from utils.Constants import Constants
from utils.log.Log import Log
from utils.exception.ExecutionException import ExecutionException

import os, shutil, commands
from utils.FileUtils import FileUtils

class CompareStatisticsProcessor( Processor):

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
        
        return None



    # --------------------------------------------------------------------------------------
    # Returns a name that will be used as display name in the user friendly outputs
    @staticmethod
    def getDisplayName():
        
        return "Comparison of branches statistics"
        


    #------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return None



    # --------------------------------------------------------------------------------------
    # Execute the processor
    def execute( self, input_commstructs):
    
        if input_commstructs == None or len( input_commstructs) == 0:
            raise ExecutionException( "CompareStatisticsProcessor.execute : No inputs")
        
        self.compareMotifHistogram( input_commstructs)
        
    
    # --------------------------------------------------------------------------------------
    # Build histogram grouping motif histograms
    def compareMotifHistogram(self, input_commstructs):
        
        # Retrieve the required parameters
        RSAT_PATH = self.component.getParameter( Constants.RSAT_DIR_PARAM)
        
        # List the motifs that are commons to all input CommStructs
        common_motifs = self.getCommonMotifs( input_commstructs)
        
        # List the label of each origins
        labels = self.getLabels( input_commstructs)
        
        number_inputs = len( input_commstructs)
        
        # Prepare the output directories
        dir_path = os.path.join( self.component.outputDir, self.component.getComponentPrefix())
        shutil.rmtree( dir_path, True)
        FileUtils.createDirectory( dir_path, 0777)
        
        for motif_name in common_motifs.keys():
            motifs = common_motifs[motif_name]
            full_histogram = {}
            full_null_histogram = {}
            for index in range( len( motifs)):
                motif = motifs[index]
                if motif != None:
                    # Add entries for the motif histogram for the current previous output
                    motif_histogram = input_commstructs[index].motifToStatistics[ motif].histogram
                    if motif_histogram != None:
                        for token in motif_histogram:
                            xy = token.split( MotifStatisticsCommStruct.HISTOGRAM_VALUE_SEPARATOR_CHAR)
                            x = xy[0]
                            y = self.getTokenAsfloat( xy[1])
                            if not x in full_histogram.keys():
                                full_histogram[x] = [0.0]*number_inputs
                            full_histogram[x][index] = y
                    # Add entries for the motif null histogram for the current previous output
                    motif_null_histogram = input_commstructs[index].motifToStatistics[ motif].nullHistogram
                    if motif_null_histogram != None:
                        for token in motif_null_histogram:
                            xy = token.split( MotifStatisticsCommStruct.HISTOGRAM_VALUE_SEPARATOR_CHAR)
                            x = xy[0]
                            y = self.getTokenAsfloat( xy[1])
                            if not x in full_null_histogram.keys():
                                full_null_histogram[x] = [0.0]*number_inputs
                            full_null_histogram[x][index] = y
            
            if len( full_histogram) > 0:
                try:
                    # Output the histogram values to file
                    file_name = motif_name + "_full_histogram"
                    file_path = os.path.join( dir_path, file_name + ".tab")
                    file = open( file_path, "w")
                    # Write the headers in the file
                    file.write("# x")
                    for label in labels:
                         file.write( "\t" + motif_name + " (" + label + ")")
                    for label in labels:
                         file.write( "\tHomogeneous model (" + label + ")")
                    file.write("\n")
                    # Write the data in the file
                    while len( full_histogram) > 0:
                        # Search for the littlest x in the dictionnary keys
                        x_val_min = 10000
                        for x in full_histogram.keys():
                            x_val = self.getTokenAsfloat(x)
                            if x_val < x_val_min:
                                x_min = x
                                x_val_min = x_val
                        file.write( x_min)
                        # write the y values of motif histograms corresponding to the min x found
                        for y in full_histogram[ x_min]:
                            file.write( "\t" + str( y))
                        del full_histogram[ x_min]
                        # write the y values of null histograms corresponding to the min x found
                        if x_min in full_null_histogram.keys():
                            for y in full_null_histogram[ x_min]:
                                file.write( "\t" + str( y))
                            del full_null_histogram[ x_min]
                        else:
                            for y in range( number_inputs):
                                file.write( "\t0.0")
                        file.write( "\n")

                    file.flush()
                    file.close()
                    
                    # Draw the histogram graph
                    graph_path = os.path.join( dir_path, file_name + "_graph.png")
                    value_cols = ""
                    for index in range( number_inputs):
                        value_cols += str( index + 2) + ","
                    for index in range( number_inputs):
                        value_cols += str( number_inputs + index + 2) + ","
                    value_cols = value_cols[:-1]
                    cmd = os.path.join( RSAT_PATH, "perl-scripts/XYgraph")
                    cmd += " -i " + file_path
                    cmd += " -title1 'Global distribution over peaks for " + motif_name +"'"
                    cmd += " -xcol 1 -ycol " + value_cols
                    cmd += " -xleg1 'Position against peak maximum' -lines"
                    cmd += " -yleg1 'Number of occurence'"
                    cmd += " -legend -header -format png -histo"
                    cmd += " -o " + graph_path
                    
                    cmd_result = commands.getstatusoutput( cmd)
                    if cmd_result[0] != 0:
                        Log.log( "CompareStatisticsProcessor.compareMotifHistogram : status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
                        Log.log( "  Command output is = \n" + str( cmd_result[1]))
                        continue
                except IOError, io_exce:
                    raise ExecutionException( "CompareStatisticsProcessor.compareMotifHistogram : Unable to save histogram to tab file : '" + file_path + "'. From:\n\t---> " +str( io_exce))


            
        
        
    # --------------------------------------------------------------------------------------
    # Returns a list the motifs that are common to all the given inpu CommStructs
    def getCommonMotifs(self, input_commstructs):
                
        common_motifs =  {}
        length = len( input_commstructs)
        for index in range( length):
            for motif in input_commstructs[index].motifList:
                if not motif.name in common_motifs.keys():
                    common_motifs[ motif.name] = [ None]*length
                common_motifs[ motif.name][index] = motif
                
        return common_motifs


    # --------------------------------------------------------------------------------------
    # Return the list of lable to be used in the histogram graphs
    def getLabels( self, input_commstructs):
        
        labels =[]
        for input_commstruct in input_commstructs:
            label = input_commstruct.origin.processorBranch
            labels.append( label)
        
        return labels
