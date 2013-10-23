
import os, commands, shutil, Queue, time, threading, gc, math, random

from manager.ProgressionManager import ProgressionManager

from utils.log.Log import Log
from utils.Constants import Constants
from utils.FileUtils import FileUtils
from utils.MotifUtils import MotifUtils
from utils.RSATUtils import RSATUtils
from utils.exception.ParsingException import ParsingException

from processor.Processor import Processor
from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct
from processor.io.BedSeqAlignmentStatsCommStruct import MotifStatistics

from common.PWM import PWM
from common.Motif import Motif

from utils.exception.ExecutionException import ExecutionException

# This processor aims to identify the motifs hidden in the provided PWM matrices coming from multiple alignements
#
# Parameters:
#   Method : the method used to obtain the result. Value should be "RSAT compare-matrices" or "TOMTOM"
#   For method "RSAT compare-matrices":
#       MotifDatabasePath : the path to the motifs files
#       MotifDatabaseFileList : list of the name of the motif files to be used as comparison
#       MotifDatabaseFormatList : list of the format of the motif files listed in the MotifDatabaseFileList parameter (one to one)
#       CustomMotifDatabaseFile (optional) : path to the motif database file provided by the user
#       CustomMotifDatabaseFormat (optional) : format of the motif database file provided by the user
#       DesiredSpeciesList (optional): the species to take into account in the multiple alignments
#       CorrelationLimit (optional) : the selection minimum threshold of the normalized correlation between reference and query motif (float number beetween 0 and 1)
#       CommandOptions (optional) : specify option to pass to the compare-matrices command
#   For method "TOMTOM":
#       MotifDatabasePath : the path to the motifs files
#       MotifDatabaseFileList : list of the name of the motif files to be used as comparison
#       CustomMotifDatabaseFile (optional) : path to the motif database file provided by the user
#       DesiredSpeciesList (optional): the species to take into account in the multiple alignments
#       CommandOptions (optional) : specify option to pass to the TOMTOM command

class MotifProcessor( Processor):
    
    METHOD_PARAM = "Method"
    METHOD_VALUE_RSAT = "rsat compare-matrices"
    METHOD_VALUE_TOMTOM = "tomtom"
    
    # Parameters for RSAT compare-matrices and TOMTOM
    MOTIF_DATABASE_PATH_PARAM = "MotifDatabasePath"
    MOTIF_DATABASE_FILE_LIST_PARAM = "MotifDatabaseFileList"
    MOTIF_PERMUTED_DATABASE_FILE_LIST_PARAM = "MotifPermutedDatabaseFileList"
    CUSTOM_MOTIF_DATABASE_FILE_PARAM = "CustomMotifDatabaseFile"
    REFERENCE_SPECIES_PARAM = "ReferenceSpecies"
    DESIRED_SPECIES_LIST_PARAM = "DesiredSpeciesList"
    COMMAND_OPTIONS_PARAM = "CommandOptions"
    THREAD_NUMBER_PARAM = "ThreadNumber"
    REPORT_MOTIF_PWM = "ReportMotifPWM"

    # Parameters for RSAT compare-matrices only
    MOTIF_DATABASE_FORMAT_LIST_PARAM = "MotifDatabaseFormatList"
    CUSTOM_MOTIF_DATABASE_FORMAT_PARAM = "CustomMotifDatabaseFormat"
    CORRELATION_LIMIT_PARAM = "CorrelationLimit"

    # Parameter for multi-threading
    THREAD_CHECK_DELAY = 2

    # Parameter for shifting in comparison of matrices
    MAXIMAL_SHIFT = 3
    
    # Parameter indicating the max number of blocks sent to a single processor.
    # In case of multi-threads, this number is devided over all threads 
    TOTAL_CONCURENT_BLOCKS = 128
    

    # ---------------------------------------------------------------------------------------------
    def __init__( self):
        Processor.__init__( self)
        self.threadLock = threading.Lock()


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
        
        return "Identification of motifs in conserved regions"
        


    #------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return ( MotifProcessor.METHOD_PARAM,\
                MotifProcessor.MOTIF_DATABASE_PATH_PARAM,\
                MotifProcessor.MOTIF_DATABASE_FILE_LIST_PARAM,\
                MotifProcessor.REFERENCE_SPECIES_PARAM)


    
    # ---------------------------------------------------------------------------------------------
    # Execute the processor
    def execute( self, input_commstructs):
    
        if input_commstructs == None or len( input_commstructs) == 0:
            raise ExecutionException( "MotifProcessor.execute : No inputs")
        
        input_commstruct = input_commstructs[0]
    
        # Retrieve the processor method to use
        method = self.getParameter( MotifProcessor.METHOD_PARAM).lower()
        
        # Retrieve the method parameters
        arguments = self.getMethodParameters( method)
        
        # Prepare the processor output dir
        self.out_path = os.path.join( self.component.outputDir, self.component.getComponentPrefix())
        shutil.rmtree( self.out_path, True)
        FileUtils.createDirectory( self.out_path, 0777)
        
        # build the output CommStruct
        output_commstruct = BedSeqAlignmentStatsCommStruct()
        output_commstruct.baseSpecies = input_commstruct.baseSpecies
        output_commstruct.bedSequencesDict = input_commstruct.bedSequencesDict
        output_commstruct.bedToMA = input_commstruct.bedToMA
        output_commstruct.motifStatistics = input_commstruct.motifStatistics
        output_commstruct.paramStatistics = input_commstruct.paramStatistics
        
        # Execute the analysis using the chosen method
        data_stats = self.executeTool( method, arguments, output_commstruct)
    
        # Finalize motif statistics
        self.finalizeStatistics( data_stats, arguments, output_commstruct)
        
        return output_commstruct
        
        
    # ---------------------------------------------------------------------------------------------
    # Retrieved the required parameters for the RSAT compare-matrices tool
    def getMethodParameters(self, method):

        arguments = {}

        # Retrieve the main path the motif databases
        db_base_path = self.getParameter( MotifProcessor.MOTIF_DATABASE_PATH_PARAM)
        
        # Retrieve the list of motif database files to use
        database_file_line = self.getParameter( MotifProcessor.MOTIF_DATABASE_FILE_LIST_PARAM)
        if database_file_line != None and not database_file_line.isspace():
            file_list = database_file_line.split()
            arguments[ MotifProcessor.MOTIF_DATABASE_FILE_LIST_PARAM] = []
            for file_path in file_list:
                arguments[ MotifProcessor.MOTIF_DATABASE_FILE_LIST_PARAM].append( os.path.join( db_base_path, file_path))
        else:
            raise ExecutionException( "MotifProcessor.getMethodParameters : No motif database file specified in parameter '" + MotifProcessor.MOTIF_DATABASE_FILE_LIST_PARAM + "'")

        # Add the custom motif database files if any
        custom_database_file_line = self.getParameter( MotifProcessor.CUSTOM_MOTIF_DATABASE_FILE_PARAM, False)
        if custom_database_file_line != None and not custom_database_file_line.isspace():
            arguments[ MotifProcessor.MOTIF_DATABASE_FILE_LIST_PARAM].append( custom_database_file_line)

        # Retrieve the list of desired species
        referenceSpecies = self.getParameter( MotifProcessor.REFERENCE_SPECIES_PARAM)
        desired_species_line = self.getParameter( MotifProcessor.DESIRED_SPECIES_LIST_PARAM, False)
        desiredSpeciesList = []
        desiredSpeciesList.append( referenceSpecies)
        if desired_species_line != None:
            desiredSpeciesList.extend( desired_species_line.split())
        arguments[ MotifProcessor.DESIRED_SPECIES_LIST_PARAM] = desiredSpeciesList

        # Retrieve the tool command line options
        command_options_line = self.getParameter( MotifProcessor.COMMAND_OPTIONS_PARAM, False)
        if command_options_line == None:
            arguments[ MotifProcessor.COMMAND_OPTIONS_PARAM] = ""
        else:
            arguments[ MotifProcessor.COMMAND_OPTIONS_PARAM] = command_options_line

        # Retrieve the Thread Number parameter
        arguments[ MotifProcessor.THREAD_NUMBER_PARAM] = self.getParameterAsint( MotifProcessor.THREAD_NUMBER_PARAM, False)
        if arguments[ MotifProcessor.THREAD_NUMBER_PARAM] == None or arguments[ MotifProcessor.THREAD_NUMBER_PARAM] <= 0:
            arguments[ MotifProcessor.THREAD_NUMBER_PARAM] = 1

        # Retrieve the parameter telling if whether or not reporting of motif PWM is required
        report_pwm = self.getParameterAsint( MotifProcessor.REPORT_MOTIF_PWM, False)
        if report_pwm == None:
            arguments[  MotifProcessor.REPORT_MOTIF_PWM] = False
        else:
            arguments[  MotifProcessor.REPORT_MOTIF_PWM] = (report_pwm.lower == "true")

        # Retrieve parameters specific to RSAT compare-matrices tool
        if method == MotifProcessor.METHOD_VALUE_RSAT:
            
            # retrieve the list of motif database file format
            database_format_line = self.getParameter( MotifProcessor.MOTIF_DATABASE_FORMAT_LIST_PARAM)
            if database_format_line != None and not database_format_line.isspace():
                arguments[ MotifProcessor.MOTIF_DATABASE_FORMAT_LIST_PARAM]  = database_format_line.split()
            else:
                raise ExecutionException( "MotifProcessor.getMethodParameters : No motif database format specified in parameter '" + MotifProcessor.MOTIF_DATABASE_FORMAT_LIST_PARAM + "'")

            # add the format of the custom motif database if any
            custom_database_format_line = self.getParameter( MotifProcessor.CUSTOM_MOTIF_DATABASE_FORMAT_PARAM, False)
            if custom_database_format_line != None and not custom_database_format_line.isspace():
                arguments[ MotifProcessor.MOTIF_DATABASE_FORMAT_LIST_PARAM].append( custom_database_format_line)

            if len( arguments[ MotifProcessor.MOTIF_DATABASE_FILE_LIST_PARAM]) != len( arguments[ MotifProcessor.MOTIF_DATABASE_FORMAT_LIST_PARAM] ):
                for added_format in range( len( arguments[ MotifProcessor.MOTIF_DATABASE_FILE_LIST_PARAM]) - len( arguments[ MotifProcessor.MOTIF_DATABASE_FORMAT_LIST_PARAM])): 
                    arguments[ MotifProcessor.MOTIF_DATABASE_FORMAT_LIST_PARAM].append( "tf")

            # retrieve the correlation threshold value
            arguments[ MotifProcessor.CORRELATION_LIMIT_PARAM] = self.getParameterAsfloat( MotifProcessor.CORRELATION_LIMIT_PARAM, False)
            if arguments[ MotifProcessor.CORRELATION_LIMIT_PARAM] == None:
                arguments[ MotifProcessor.CORRELATION_LIMIT_PARAM] = 0.5
        
        # Retrieve parameters specific to MEME TOMTOM
        if method == MotifProcessor.METHOD_VALUE_TOMTOM:
            arguments[ MotifProcessor.MOTIF_DATABASE_FORMAT_LIST_PARAM] = ['meme'] * len( arguments[ MotifProcessor.MOTIF_DATABASE_FILE_LIST_PARAM])
            
    
        return arguments

    
    # ---------------------------------------------------------------------------------------------
    # Search for the motifs using the RSAT compare-matrices tool
    def executeTool(self, method, arguments, output_commstruct):

        # Prepare the processor output dir
        out_path = os.path.join( self.component.outputDir, self.component.getComponentPrefix())
        shutil.rmtree( out_path, True)
        FileUtils.createDirectory( out_path, 0777)

        # Retrieve the list of motif databases to use
        database_file_list = arguments[ MotifProcessor.MOTIF_DATABASE_FILE_LIST_PARAM]
        database_format_list = arguments[ MotifProcessor.MOTIF_DATABASE_FORMAT_LIST_PARAM]        
        
        # Permute the motif matrices in the database files for the hypergeometric test
        permuted_database_file_list = {}
        for index in range( len( database_file_list)):
            permuted_list = RSATUtils.permuteMatrix( database_file_list[ index], database_format_list[ index], out_path, 100)
            if permuted_list != None:
                permuted_database_file_list[ database_file_list[ index]] = permuted_list
                
        arguments[ MotifProcessor.MOTIF_PERMUTED_DATABASE_FILE_LIST_PARAM] = permuted_database_file_list

        # Generates the list of unrecognized motifs
        motif_list_infos = self.getMotifList( output_commstruct, arguments[ MotifProcessor.DESIRED_SPECIES_LIST_PARAM])
        motif_list = motif_list_infos[ 0]
        
        # Initialize parameters
        motif_pack_length = int( MotifProcessor.TOTAL_CONCURENT_BLOCKS / float( arguments[ MotifProcessor.THREAD_NUMBER_PARAM]))
        old_working_dir = os.getcwd()
        
        # Build the motif list queue (by pack)
        motif_pack_queue =Queue.Queue()
        queue_initial_length = 0
        for index in range( 0, len( motif_list), motif_pack_length):
            index_end = index + motif_pack_length
            if index_end > len( motif_list):
                index_end = len( motif_list)
            motif_pack_queue.put( motif_list[index:index_end])
            queue_initial_length+=1

        # Launch the firsts Threads
        ProgressionManager.setTaskProgression( "Identifying Motifs",  self.component, 0.0)
        thread_list =[]
        count_thread = 0
        while not motif_pack_queue.empty() and len( thread_list) < arguments[ MotifProcessor.THREAD_NUMBER_PARAM]:
            new_thread = self.startNewThread( count_thread, method, arguments, motif_pack_queue, output_commstruct, out_path)
            if new_thread != None:
                thread_list.append( new_thread)
                count_thread += 1

        # Manage the Threads list in order to empty the file queue
        progress = 0
        while not motif_pack_queue.empty() or len( thread_list) > 0:
            threads_to_remove = []
            threads_to_add = []
            for thread in thread_list:
                if not thread.is_alive():
                    # remove the thread from the thread list and update progression
                    threads_to_remove.append( thread)
                    # start a new thread for the next motif pack
                    new_thread = self.startNewThread( count_thread, method, arguments, motif_pack_queue, output_commstruct, out_path)
                    if new_thread != None:
                        threads_to_add.append( new_thread)
                        count_thread += 1
                    progress += 1
                    ProgressionManager.setTaskProgression( "Identifying Motifs", self.component, progress/float( queue_initial_length))
            
            # remove from the list the threads that were detected as not alive
            if len( threads_to_remove) > 0:
                for th in threads_to_remove:
                    thread_list.remove( th)
                    th = None
                gc.collect()
            
            # add to the list the thread that were started
            if len( threads_to_add) > 0:
                for th in threads_to_add:
                    thread_list.append( th)          
                Log.trace( "MAFProcessor.executeTool : Starting new thread to identify motifs. Number of active Thread = " + str( len( thread_list)))

            # wait before testing again the threads status
            time.sleep( MotifProcessor.THREAD_CHECK_DELAY)

        # Change dir to previous working dir
        os.chdir( old_working_dir)

        return motif_list_infos



    # ---------------------------------------------------------------------------------------------
    # Set the motif details (family, type, class) in the MotifStatistics and 
    # compute the p-value with the hyperbolic test
    def finalizeStatistics(self, data_stats, arguments, output_commstruct):
        
        # Get the motif details from Jaspar database
        motif_details = MotifUtils.getMotifsDetailsFromJaspar()
        
        motif_number_in_db = 0
        for database_path in arguments[ MotifProcessor.MOTIF_DATABASE_FILE_LIST_PARAM]:
            current_number = MotifUtils.getMotifsNumberFromTransfac( database_path)
            if current_number != None:
                motif_number_in_db = motif_number_in_db + current_number
        
        for motif_name in output_commstruct.motifStatistics.keys():
            motif_stats = output_commstruct.motifStatistics[ motif_name]
            
            # set the motif details (family, type, class)
            self.setMotifStatsDetails( motif_stats, motif_details)
            
            # compute hyperbolic test p-value
            self.computeHypergeometricPValue( motif_stats, data_stats, motif_number_in_db)



    # ---------------------------------------------------------------------------------------------
    # Set the motif details (family, type, class) in the given MotifStatistic
    def setMotifStatsDetails(self, motif_stats, motif_details):
        
        motif_ids = motif_details[ 0]
        motif_families = motif_details[ 1]
        motif_types = motif_details[ 2]
        motif_classes = motif_details[ 3]
        
        motif_name = motif_stats.motifName
        
        if not motif_name in motif_ids.keys():
            return
        
        motif_id = motif_stats.motifID
        if motif_id != motif_ids[ motif_name]:
            Log.log( "MotifStatisticsProcessor.buildStatistics : A motif have an invalid ID : Motif ID = " + motif_id + " Motif Jaspar ID = " + motif_ids[ motif_name])
            return
        
        if not motif_name in motif_families.keys():
            Log.log( "MotifStatisticsProcessor.buildStatistics : A motif has no family defined in Jaspar: Motif name = " + motif_name)
        else:
            motif_stats.motifFamily = motif_families[ motif_name]
        
        if not motif_name in motif_types.keys():
            Log.log( "MotifStatisticsProcessor.buildStatistics : A motif has no type defined in Jaspar: Motif name = " + motif_name)
        else:
            motif_stats.motifType = motif_types[ motif_name]

        if not motif_name in motif_classes.keys():
            Log.log( "MotifStatisticsProcessor.buildStatistics : A motif has no class defined in Jaspar: Motif name = " + motif_name)
        else:
            motif_stats.motifClass = motif_classes[ motif_name]



    # ---------------------------------------------------------------------------------------------
    # Prepare the computation of the p-value with the hyperbolic test
    def computeHypergeometricPValue( self, motif_stats, data_stats, motif_number_in_db): 

        if not motif_stats.hasAttribute( MotifStatistics.MOTIF_HIT_SCORE):
            Log.log( "MotifProcessor.computeHypergeometricPValue : the motif " + motif_stats.motifName + " has not been identified with non-permuted matrix")
            return
            
        if not motif_stats.hasAttribute( MotifStatistics.MOTIF_HYP_HIT_SCORE):
            motif_stats.setAttribute( MotifStatistics.MOTIF_HYP_HIT_SCORE, 1.0)
            
        input_motif_size_list = data_stats[ 1]
        
        output_motif_length = motif_stats.motifSize
        
        # Compute the total number of "balls" drawned :
        # k = number of motif hit with non permuted matrix + number of motif hit with permuted matrix
        # if the motif has no hits with no permuted matrix, bypass it
        k = motif_stats.getAttributeAsint( MotifStatistics.MOTIF_HIT_SCORE)
        if k == 0:
            return
        k += motif_stats.getAttributeAsint( MotifStatistics.MOTIF_HYP_HIT_SCORE)
        
        # compute the total number of success "balls"
        # m = total number of possible motif start index in the conserved regions (input motifs)
        # This number has to take into account position that have been avoided for algorithm reasons (uncount)
        possible_position = 0
        for input_size in input_motif_size_list:
            length_to_add = input_size + output_motif_length - 2*MotifProcessor.MAXIMAL_SHIFT + 1
            if length_to_add > 0:
                possible_position += length_to_add

	    m = (2 * possible_position) - motif_stats.getAttributeAsint( MotifStatistics.MOTIF_UNCOUNT)

        # compute the total number of failure "balls"
        # n = total number of possible motif start index in the conserved regions (input motifs)
        # This number has to take into account position that have been avoided for algorithm reasons (uncount)
        n = (2 * possible_position) - motif_stats.getAttributeAsint( MotifStatistics.MOTIF_HYP_UNCOUNT)
        
        # Compute the total number of succes "balls" drawned :
        # x = number of motif hit with non permuted matrix
        x = motif_stats.getAttributeAsint( MotifStatistics.MOTIF_HIT_SCORE)
        
        # Compute the hypergeometric p-value
        #print "MOTIF = " + motif_stats.motifName
        #print "total number of balls in urn (m+n)=" + str( m+n)
        #print "black ball in urn (m)=" + str( m)
        #print "number of balls drawned (k)=" + str( k)
        #print "number of drawned black balls (x)=" + str( x) 
        #print "TO min( m, k)= " + str( min( m, k))
        p_value = self.hypergeometric( m, m+n, k, x, min( m, k))
        e_value = p_value * motif_number_in_db
        
        motif_stats.setAttribute( MotifStatistics.MOTIF_HYP_PVALUE, "%(number)2e" %{ "number" : p_value})
        motif_stats.setAttribute( MotifStatistics.MOTIF_HYP_EVALUE, "%(number)2e" %{ "number" : e_value})
        

    # ---------------------------------------------------------------------------------------------
    # Compute the p-value with the hypergeometric test
    def hypergeometric(self, var_m, var_n, var_k, var_x, var_to):

        var_w = var_n - var_m ### white balls in the urn

        # initialization
        var_proba = 0
        var_log_proba = 0

        # incompatible parameter values
        if var_k > var_n:
            Log.log("Sample (var_k) cannot be larger than number of balls in the urn (var_n)")
            return 1.0

        if var_m > var_n:
            Log.log("Number of black balls in the urn (var_m) cannot be larger than number of balls in the urn (var_n)")
            return 1.0

        if var_w > var_n:
            Log.log("Number of white balls in the urn (var_m) cannot be larger than number of balls in the urn (var_n)")
            return 1.0

        if var_x > var_k:
            Log.log("Number of black balls in the sample (var_x) cannot be larger than sample size (var_k)")
            return 1.0

        if var_x < 0:
            Log.log( "Number of black balls in the sample (var_x) must be strictly positive\n")
            return 1.0

        # In the case of impossible events, log is returned
        if var_x > var_m:
            Log.log("Number of black balls in the sample (var_x) cannot be larger than number of black balls in the urn (var_m)")
            return 1.0

        if var_k - var_x > var_w:
            Log.log( "".join( ("", "Number of white balls in the sample (", var_k-var_x,") cannot be larger than number of white balls in the urn (", var_w,")")))
            return 1.0

        if var_w < var_k:
            # If the number of non-marked balls (var_w) is lower than the size
            # of the selection (var_k), there will be at least var_k - var_m marked
            # balls in the selection (var_x >= var_min_x = var_k - var_m) -> the proba
            # of all values from 0 to (var_min_x - 1) is null.
            
            var_min_x = var_k - var_w
            var_log_proba = 0
            for var_i in range( var_n - var_k + 1, var_n + 1):
                var_log_proba -= math.log(var_i)

            for var_i in range( 1, var_min_x+1):
                var_log_proba -= math.log(var_i)
              
            for var_i in range( var_m - var_min_x + 1, var_m + 1):
                var_log_proba += math.log(var_i)

            for var_i in range(var_k - var_min_x + 1, var_k + 1):
                var_log_proba += math.log(var_i)

            for var_i in range(var_w -var_k + var_min_x + 1, var_w +1):
                var_log_proba += math.log(var_i)
    
            var_start = var_min_x + 1

        else:
            # calculate value for 0 successes
            var_log_proba = 0
            for var_i in range(var_w - var_k + 1, var_n - var_m +1):
                var_log_proba += math.log(var_i)
                
            for var_i in range(var_n - var_k + 1, var_n +1):
                var_log_proba -= math.log(var_i)
                
            var_start = 1

        var_proba = self.logToEng( var_log_proba)

        # recursive calculation of the hypergeometric density for var_x
        # (if cumulative, var_x is the first value of the sum)
        if var_start <= var_x:
            for var_i in range( var_start, var_x +1):
                var_log_proba += math.log(var_m - var_i + 1)
                var_log_proba -= math.log(var_i)
                var_log_proba += math.log(var_k - var_i + 1)
                var_log_proba -= math.log(var_w - var_k + var_i)

            var_proba = self.logToEng(var_log_proba)

        # If var_to > var_x (cumulative density has to be computed),
        # pursue the recursive computation and add up the values to the sum
        for var_i in range(var_x + 1, var_to + 1):
            var_last_proba = var_proba
            var_log_proba += math.log(var_m - var_i + 1)
            var_log_proba -= math.log(var_i)
            var_log_proba += math.log(var_k - var_i + 1)
            var_log_proba -= math.log(var_w - var_k + var_i)
            var_new_proba = self.logToEng(var_log_proba)
            var_proba = var_last_proba + var_new_proba
            # beyond precision limit, it is worthless adding more elements
            if var_proba == var_last_proba and var_proba >0:
                break 

        var_proba = min( var_proba, 1) ### floating point calculation errors

        # var_proba = exp(var_log_proba);
        return var_proba



    # ---------------------------------------------------------------------------------------------
    # Convert log to exponentiel expression
    def logToEng( self, var_log):
        ''''
        var_base = 10
        var_log_base = math.log( var_base)
        var_log /= math.log(10)
        var_eng = 10**(1+var_log - int(var_log))
        var_eng += "E"
        if int(var_log)-1 > 0:
            var_eng += "+";
        var_eng += int(var_log)-1
        '''
        var_base = 10
        var_log_base = math.log( var_base)
        var_log /= math.log(10)
        var_eng = 10**(1+var_log - int(var_log))
        str_var_eng = str( var_eng) + "E"
        if int(var_log)-1 > 0:
            str_var_eng += "+"
        str_var_eng += str(int(var_log)-1)
        
        return float( str_var_eng)


    # #######################################
    #
    #     RSAT COMPARE-MATRICES FUNCTIONS
    #
    # #######################################    


    # ---------------------------------------------------------------------------------------------
    # Search for the motifs using the RSAT compare-matrices tool
    def launchCompareMatrices(self, arguments, sub_motif_list, output_commstruct, out_path):

        thread_name = threading.currentThread().getName()

        # Prepare the input file and output directory for rsat compare-matrices
        file_info = self.outputMotifListToTransfacFile( sub_motif_list, out_path, thread_name)
        
        # Execute the motif identification
        Log.trace( "MotifProcessor.launchCompareMatrices : " + thread_name + " : Sending " + str( len( sub_motif_list)) + " conserved regions to identification")
        self.executeCompareMatrices( False, output_commstruct, "results", file_info, arguments)
        
        # Execute the motif identification for the hypergeometric test
        Log.trace( "MotifProcessor.launchCompareMatrices : " + thread_name + " : Sending " + str( len( sub_motif_list)) + " conserved regions to hypergeometric test")
        self.executeCompareMatrices( True, output_commstruct, "hyp_results", file_info, arguments)
        
        # Remove the RSAT compare-matrices result dir
        os.chdir( out_path)
        shutil.rmtree( file_info[0], True)

    
    
    # ---------------------------------------------------------------------------------------------
    # Execute the comparison with all the listed databases
    def executeCompareMatrices( self, hypergeometric_test, output_commstruct, prefix, file_info, arguments):
        
        # Retrieve the method arguments
        if hypergeometric_test == False:
            database_file_list = []
            for database_path in arguments[ MotifProcessor.MOTIF_DATABASE_FILE_LIST_PARAM]:
                database_file_list.append( database_path)
        else:
            # In case of hypergeometric test, choose one of the database of permuted motifs previously computed
            permuted_databases = arguments[ MotifProcessor.MOTIF_PERMUTED_DATABASE_FILE_LIST_PARAM]
            database_file_list = []
            for initial_database in permuted_databases.keys():
                permuted_databases_list = permuted_databases[ initial_database]
                rand_index = int( random.uniform( 0, len( permuted_databases_list)))
                database_file_list.append( permuted_databases_list[ rand_index])
                
        database_format_list = arguments[ MotifProcessor.MOTIF_DATABASE_FORMAT_LIST_PARAM]
        command_options = arguments[ MotifProcessor.COMMAND_OPTIONS_PARAM]
        correlation_limit = arguments[ MotifProcessor.CORRELATION_LIMIT_PARAM]
        
        # get the transfac file of conserved region matrices information
        dir_path = file_info[0]
        file_path = file_info[1]
        
        # Change directory to output dir
        os.chdir( dir_path)

        RSAT_PATH = self.component.getParameter( Constants.RSAT_DIR_PARAM)
        # Parse the list of motif databases and apply each of them to the compare-matrices tool
        # on the given file of conserved blocks
        all_input_motif_to_remove = {}
        for db_index in range( len( database_file_list)):
            sub_prefix = prefix + "_" + str( db_index)
            # Compose the compare-matrices command line with all required options
            Log.trace( "MotifProcessor.executeCompareMatrices : " + threading.currentThread().getName() + " sub_prefix before execution = " + sub_prefix)
            
            # Command for the PERL version of compare-matrices 
	    #cmd = os.path.join( RSAT_PATH , "perl-scripts/compare-matrices")
            #cmd += " -file1 " + file_path + " -format1 tf"
            #cmd += " -file2 " + database_file_list[db_index] + " -format2 " + database_format_list[db_index]
            #cmd += " -mode profiles"
            #cmd += " -lth w " + str( MotifProcessor.MAXIMAL_SHIFT)
            #cmd += " -lth ncor2 0.7"
            #cmd += " -return matrix_id,cor,Ncor2,w,consensus,offset" 
            #cmd += " -o " + sub_prefix
            #cmd += " " + command_options

            # Command for the C version of compare-matrices
            cmd = os.path.join( RSAT_PATH , "contrib/peak-footprints/tools/compare-matrices")
            cmd += " -file1 " + file_path
            cmd += " -file2 " + database_file_list[db_index]
            cmd += " -lth_w " + str( MotifProcessor.MAXIMAL_SHIFT)
            cmd += " -lth_ncor2 0.7"
            cmd += " -o " + sub_prefix + ".tab"	
            cmd += " > " + sub_prefix + ".log"	
            
            Log.info( "MotifProcessor.executeCompareMatrices : Starting comparison with database : " + database_file_list[db_index])
            Log.info( "MotifProcessor.executeCompareMatrices : command used is : " + cmd)
            
            # Execute the command
            cmd_result = commands.getstatusoutput( cmd)
            Log.trace( "MotifProcessor.executeCompareMatrices : " + threading.currentThread().getName() + " (" + sub_prefix + ") : status returned is :" + str( cmd_result[0]))
            if cmd_result[0] != 0:
                Log.log( "MotifProcessor.executeCompareMatrices : status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
                Log.log( "MotifProcessor.executeCompareMatrices : command output is = \n" + str( cmd_result[1]))
                continue
                
            # Convert the tab result to html file
            commands.getstatusoutput( os.path.join( RSAT_PATH , "perl-scripts/text-to-html") +" -i " + sub_prefix + ".tab -o " + sub_prefix + ".html")
            
            # Verify if the compare-matrices is present and have finished to be written
            result_file_path = os.path.join( dir_path, sub_prefix + ".tab")
            Log.trace( "MotifProcessor.executeCompareMatrices : " + threading.currentThread().getName() + " sub_prefix after execution = " + sub_prefix)
            if not self.verifyResultFile( result_file_path):
                Log.log( "MotifProcessor.executeCompareMatrices : Compare-matrices result file is not accessible at " + result_file_path)
                continue
            
            # Parse the result of the compare-matrices command to get the result list
            result = self.parseCompareMatricesResult( os.path.join( dir_path, sub_prefix + ".tab"))
            if result == None:
                Log.log( "MotifProcessor.executeCompareMatrices : Compare-matrices result file gave a null result :" + result_file_path)
                continue
	    else:
		Log.trace("MotifProcessor.executeCompareMatrices : " + threading.currentThread().getName() + " : Number of motif read = " + str(len(result))) 
                
            # Analyze the result list and store the conserved elements in the thread_result of the analysis
            if hypergeometric_test == False:
                input_motif_to_remove = self.analyseCompareMatricesResult( output_commstruct, result, correlation_limit)
                # Update the list of input motifs to remove with the list provided by the last analysis 
                for alignment in input_motif_to_remove.keys():
                    if alignment in all_input_motif_to_remove.keys():
                        all_input_motif_to_remove[ alignment].extend( input_motif_to_remove[ alignment])
                    else:
                        all_input_motif_to_remove[ alignment] = input_motif_to_remove[ alignment]
            else:
                self.analyseCompareMatricesHypResult( output_commstruct, result, correlation_limit)
        
        #Remove the input motifs that were associated to identified reference motifs
        for alignment in all_input_motif_to_remove.keys():
            motifs_to_remove = all_input_motif_to_remove[ alignment]
            if motifs_to_remove != None:
                for motif in motifs_to_remove:
                    if motif in alignment.motifs:
                        alignment.motifs.remove( motif)
        


    # ---------------------------------------------------------------------------------------------
    # Output the motif list to file and prepare the output directory
    def outputMotifListToTransfacFile(self, motif_list, out_path, folder_prefix):
        
        try:
            dir_path = os.path.join( out_path, folder_prefix)
            shutil.rmtree( dir_path, True)
            FileUtils.createDirectory( dir_path, 0777)
            file_name = "motifs"
            file_path = os.path.join( dir_path, file_name + ".tf")
            tf_file = FileUtils.openFile( file_path, "w", 0666)
            for motif in motif_list:
                tf_file.write( "AC\t" + motif.name + "\n")
                tf_file.write( "XX\n")
                tf_file.write( "ID\t" + motif.name + "\n")
                tf_file.write( "XX\n")
                tf_file.write( motif.pwm.convertToTransfac())
                tf_file.write( "XX\n")
                tf_file.write("//\n")
            tf_file.flush()
            tf_file.close()
            return (dir_path, file_path)
        except IOError, io_exce:
            raise ExecutionException( "MotifProcessor.outputMotifListToTransfacFile : Unable to save motifs to tab tf_file : '" + file_path + "'. From:\n\t---> " +str( io_exce))



    # ---------------------------------------------------------------------------------------------
    # Parse the tab file of the compare-matrices result directory
    def parseCompareMatricesResult(self, file_path):
        
        result = {}
        
        try:
            tf_file = open( file_path,  "r")
           
            # define some constants for the parsing
            input_motif_col = 2
            
            # load the result tab tf_file in a dictionary : key = input_motif / Value = dictionary of tf_file column header/value
            headers_ok = False
            headers = []
            count_line = 0
            for line in tf_file:
                count_line += 1
                # Read the headers of the columns
                if line[0] == "#":
                    headers = line.split()
                    headers[0] = headers[0][1:]
                    headers_ok = True
                else:
                    if headers_ok and line[0] != ";":
                        # read each column and assign it to the corresponding header
                        tokens = line.split()
                        input_motif = tokens[ input_motif_col]
                        if not input_motif in result.keys():
                            result[ input_motif] = []
                        output_dic = {}
                        if len( tokens) == len( headers):
                            for index in range( len( tokens)):
                                token = tokens[index]
                                header = headers[ index]
                                output_dic[ header] = token
                            result[ input_motif].append( output_dic)
                        else:
                            Log.log( "MotifProcessor.parseCompareMatricesResult: wrongly formatted line : #tokens = " + str( len( tokens)) + " #headers = " + str( len( headers)) + " in line " + str( count_line) + " of tf_file " + file_path)
                            
            tf_file.close()
        except (KeyError, IndexError, IOError),  io_exce:
            Log.log( "MotifProcessor.parseCompareMatricesResult : Unable to open/parse the compare-matrices result : '" + file_path + "'. From:\n\t---> " +str( io_exce))
            return None
            #raise ExecutionException( "MotifProcessor.parseCompareMatricesResult : Unable to open the compare-matrices result : '" + file_path + "'. From:\n\t---> " +str( io_exce))
 
        return result



    # ---------------------------------------------------------------------------------------------
    # Analyse the compare-matrices result and store the information in the output CommStruct
    def analyseCompareMatricesResult( self, output_commstruct, result, correlation_limit):
        
        # define some constants for the analysis
        name2_header      = "name2"
        id2_header        = "id2"
        ncor2_header      = "Ncor2"
        w2_header         = "w2"
        offset_header     = "offset"
        consensus2_header = "consensus2"
        strand_header     = "strand"
	uncount_header    = "uncounted"
        
	total_kept = 0
        motif_list = []
        input_motif_to_remove = {}
        # parse the bedseq and the bedseq alignment to get the motifs
        for bedseq in output_commstruct.bedToMA.keys():
            for alignment in output_commstruct.bedToMA[ bedseq]:
                if alignment.motifs != None or len( alignment.motifs) > 0:
                    analysis_result = {}
                    for input_motif in alignment.motifs:
                        # if the motif name is in the results, retrieve the list of result series
                        if input_motif.name in result.keys():
			    #Log.trace( "FOUND MOTIF = " + input_motif.name)
                            motif_dics = result[ input_motif.name]
                            #for each entry in the result series, create an output motif if the result ncor2 > correlation_limit
                            Log.info( "MotifProcessor.analyseCompareMatricesResult : motif_dics size for " + input_motif.name + " = " + str( len( motif_dics)))
                            for motif_dic in motif_dics:
                                try:
                                    n_cor2 = self.getTokenAsfloat( motif_dic[ ncor2_header])
                                    Log.info( "MotifProcessor.analyseCompareMatricesResult :     n_cor2 for " + str( motif_dic[ id2_header]) + " = " + str( n_cor2))
                                    if n_cor2 >= correlation_limit:
                                        # retrieve the data
                                        output_motif_offset = self.getTokenAsint( motif_dic[ offset_header])
                                        output_motif_length = self.getTokenAsint( motif_dic[ w2_header])
                                        output_motif_uncount = self.getTokenAsint( motif_dic[ uncount_header])
                                        output_motif_name = motif_dic[ id2_header]
                                        Log.info( "Detected motif = " + output_motif_name)
                                        # build the new motif
                                        output_motif = Motif( input_motif.indexStart + output_motif_offset, input_motif.indexStart + output_motif_offset + output_motif_length, output_motif_name, None )
                                        #output_motif.consensus = self.removeNoInfoChar( motif_dic[ consensus2_header])
					output_motif.consensus = ""
                                        output_motif.id = motif_dic[ name2_header]
                                        output_motif.offset = output_motif_offset
                                        output_motif.score = n_cor2
					output_motif.uncount = output_motif_uncount
                                        strand = motif_dic[ strand_header].upper()
                                        if strand == Constants.DIRECT_STRAND:
                                            output_motif.strand = Constants.POSITIVE_STRAND
                                        elif strand == Constants.REVERSE_STRAND:
                                            output_motif.strand = Constants.NEGATIVE_STRAND
                                        else:
                                            Log.log( "MotifProcessor.analyseCompareMatricesResult : Unknown strand : '" + strand + "' for result = " + str( motif_dic))
                                        if not input_motif in analysis_result.keys():
                                            analysis_result[ input_motif] = []
                                        # add the new motif to the corresponding conserved region dictionnary
                                        analysis_result[ input_motif].append( output_motif)
                                        # add the new motif to the result list
                                        motif_list.append( output_motif)
					total_kept = total_kept + 1
                                except ParsingException, par_exce:
                                    Log.log( "MotifProcessor.analyseCompareMatricesResult : Some value of a motif correlation is not a float : " + str( motif_dic) + ". From:\n\t--->" + str( par_exce))
                    
                    # add the new motifs to the alignment motif list and remove the corresponding conserved regions
                    #Log.info( " -len( analysis_result) = " + str( len( analysis_result)))
                    if len( analysis_result) > 0:
                        
                        self.threadLock.acquire()
                        
                        for input_motif in analysis_result.keys():
                            # remove the input motif from the alignement and add all the corresponding output motifs
                            if not alignment in input_motif_to_remove.keys():
                                input_motif_to_remove[ alignment] = []
                            input_motif_to_remove[alignment].append( input_motif)
                            alignment.motifs.extend( analysis_result[ input_motif])
                            # create the motif statistics if required and
                            # increment the hit score of each output motif
                            for output_motif in analysis_result[ input_motif]:
                                Log.info( " --Found output motif = " + output_motif.name)
                                output_motif_name = output_motif.name
                                if not output_motif_name in output_commstruct.motifStatistics.keys():
                                    output_commstruct.motifStatistics[ output_motif_name] = MotifStatistics( output_motif_name)
                                motif_stats = output_commstruct.motifStatistics[ output_motif_name]
                                motif_stats.motifSize = output_motif.indexEnd - output_motif.indexStart
                                motif_stats.motifID = output_motif.id
				# Add hit score to statistics
                                if motif_stats.hasAttribute( MotifStatistics.MOTIF_HIT_SCORE):
                                    score = motif_stats.getAttributeAsint( MotifStatistics.MOTIF_HIT_SCORE)
                                else:
                                    score = 0
                                motif_stats.setAttribute( MotifStatistics.MOTIF_HIT_SCORE, score + 1)
				# Add uncount position to statistics
                                if motif_stats.hasAttribute( MotifStatistics.MOTIF_UNCOUNT):
                                    score_uncount = motif_stats.getAttributeAsint( MotifStatistics.MOTIF_UNCOUNT)
                                else:
                                    score_uncount = 0
				motif_stats.setAttribute( MotifStatistics.MOTIF_UNCOUNT, score_uncount + output_motif.uncount)
                                    
                        self.threadLock.release()
        
	Log.trace( "MotifProcessor.analyseCompareMatricesResult : Total number of motif kept = " + str(total_kept))

        return input_motif_to_remove


    # ---------------------------------------------------------------------------------------------
    # Analyse the hypergeometric compare-matrices result and store the information in the output CommStruct
    def analyseCompareMatricesHypResult(self, output_commstruct, result, correlation_limit):
        
        # define some constants for the analysis
        ncor2_header      = "Ncor2"
        id2_header        = "id2"
	uncount_header    = "uncounted"
        
        self.threadLock.acquire()
        
        for input_motif_name in result.keys():
            motif_dics = result[ input_motif_name]
            #for each entry in the result series, increment the MOTIF_HYP_HIT_SCORE stat if ncor2 > correlation_limit
            for motif_dic in motif_dics:
                try:
                    n_cor2 = self.getTokenAsfloat( motif_dic[ ncor2_header])
                    if n_cor2 >= correlation_limit:
                        # If the motif have been indentified in the identification phase, increment the MOTIF_HYP_HIT_SCORE stat
                        output_motif_name_perm = motif_dic[ id2_header]
                        output_motif_name = output_motif_name_perm[:-6] ## remove the '_perm1' suffix added by the matrix permutation tool
                        if output_motif_name in output_commstruct.motifStatistics.keys():
                            motif_stats = output_commstruct.motifStatistics[ output_motif_name]
                            # Add hit score to statistics
                            if motif_stats.hasAttribute( MotifStatistics.MOTIF_HYP_HIT_SCORE):
                                score = motif_stats.getAttributeAsint( MotifStatistics.MOTIF_HYP_HIT_SCORE)
                            else:
                                score = 0
                            motif_stats.setAttribute( MotifStatistics.MOTIF_HYP_HIT_SCORE, score + 1)
                            # Add uncount position to statistics
			    output_motif_uncount = self.getTokenAsint( motif_dic[ uncount_header])
                            if motif_stats.hasAttribute( MotifStatistics.MOTIF_HYP_UNCOUNT):
                                score_uncount = motif_stats.getAttributeAsint( MotifStatistics.MOTIF_HYP_UNCOUNT)
                            else:
                                score_uncount = 0
			    motif_stats.setAttribute( MotifStatistics.MOTIF_HYP_UNCOUNT, score_uncount + output_motif_uncount)
                                
                except ParsingException, par_exce:
                    Log.log( "MotifProcessor.analyseCompareMatricesHypResult : Some value of a motif correlation is not a float : " + str( motif_dic) + ". From:\n\t--->" + str( par_exce))
        
        self.threadLock.release()


    # #######################################
    #
    #     MEME TOMTOM FUNCTIONS
    #
    # #######################################


    # ---------------------------------------------------------------------------------------------
    # Search for the motifs using the TOMTOM method    
    def launchTOMTOM(self, arguments, sub_motif_list, output_commstruct, out_path):
        
        Log.trace( "MotifProcessor.launchTOMTOM : Sending " + str( len( sub_motif_list)) + " conserved regions to identification")

        database_file_list = arguments[ MotifProcessor.MOTIF_DATABASE_FILE_LIST_PARAM] 
        command_options = arguments[ MotifProcessor.COMMAND_OPTIONS_PARAM]
    
        # Prepare the input file and output directory for rsat compare-matrices
        file_info = self.outputMotifListToMEMEFile( sub_motif_list, out_path, threading.currentThread().getName())
        dir_path = file_info[0]
        file_path = file_info[1]
        
        # Retrieve the motif ID from transfac motif database
        motifs_details = MotifUtils.getMotifsDetailsFromTransfac()
        motif_ids = motifs_details[0]
        
        # Change directory to output dir
        os.chdir( dir_path)
        
        # Execute the comparison with all the listed databases
        MEME_PATH = self.component.getParameter( Constants.MEME_DIR_PARAM)
        for index in range( len( database_file_list)):
            # Compose the TOMTOM command line with all required options
            cmd = os.path.join( MEME_PATH , "bin/tomtom")
            cmd += " -oc " + dir_path
            #cmd += " -text"
            cmd += " " + command_options
            cmd += " " + file_path 
            #cmd += " " + database_file_path
            
            Log.info( "MotifProcessor.launchTOMTOM : Starting comparison with database : " + database_file_list[index])
            Log.info( "MotifProcessor.launchTOMTOM : Command used is : " + cmd)
            
            # Execute the command
            cmd_result = commands.getstatusoutput( cmd)
            if cmd_result[0] != 0:
                Log.log( "MotifProcessor.launchTOMTOM : Status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
                Log.log( "MotifProcessor.launchTOMTOM : Command output is = \n" + str( cmd_result[1]))
                continue
                        
            # parse the result of the compare-matrices command to get the result list
            result = self.parseTOMTOMResult( os.path.join( dir_path, "tomtom.txt"))
            if result == None:
                continue
            
            # analyze the result list and store the conserved elements in the final_result of the analysis
            identified_motifs = self.analyseTOMTOMResult( output_commstruct, result, motif_ids)
            
            # Retrieve the identified motif's PWM if required
            if arguments[ MotifProcessor.REPORT_MOTIF_PWM] == True:
                self.retrieveMotifsMatrices( identified_motifs, database_file_list[index], "meme")

        # Remove the RSAT compare-matrices results
        os.chdir( out_path)
        shutil.rmtree( dir_path, True)
        
        

    # ---------------------------------------------------------------------------------------------
    # Output the list of motif to a file in MEME format
    def outputMotifListToMEMEFile(self, motif_list, out_path, folder_prefix):

        try:
            dir_path = os.path.join( out_path, folder_prefix)
            shutil.rmtree( dir_path, True)
            FileUtils.createDirectory( dir_path, 0777)
            file_name = "motifs"
            file_path = os.path.join( dir_path, file_name + ".meme")
            meme_file = open( file_path, "w")
            meme_file.write( "MEME version 4.4\n\n")
            meme_file.write( "ALPHABET= ACGT\n\n")
            meme_file.write( "strands: + -\n\n")
            meme_file.write( "Background letter frequencies (from\n")
            for letter in sorted( Constants.HG_BACKGROUND_MODEL.keys()):
                meme_file.write( letter + " " + str( Constants.HG_BACKGROUND_MODEL[ letter]) + "\t" )
            meme_file.write( "\n\n")
            meme_file.flush()
            for motif in motif_list:
                meme_file.write( motif.pwm.convertToMEME( motif.name))
                meme_file.write( "\n")
                meme_file.flush()
            meme_file.close()
            return (dir_path, file_path)
        except IOError, io_exce:
            raise ExecutionException( "MotifProcessor.outputMotifListToMEMEFile : Unable to save motifs to tab meme_file : '" + file_path + "'. From \n\t" +str( io_exce))



    # ---------------------------------------------------------------------------------------------
    # Parse the tab file of the compare-matrices result directory
    def parseTOMTOMResult(self,  file_path):
        
        result = {}
        
        try:
            tomtom_file = open( file_path,  "r")
           
            # define some constants for the parsing
            input_motif_col = 0
            
            # load the result tab tomtom_file in a dictionnary : key = input_motif / Value = dictionnary of tomtom_file column header/value
            headers_ok = False
            headers = []
            for line in tomtom_file:
                # Read the headers of the columns
                if line[0] == "#":
                    headers = line.split( "\t")
                    headers[0] = headers[0][1:]
                    headers_ok = True
                else:
                    if headers_ok:
                        tokens = line.split()
                        input_motif = tokens[ input_motif_col]
                        if not input_motif in result.keys():
                            result[ input_motif] = []
                        output_dic = {}
                        for index in range( len( tokens)):
                            token = tokens[index]
                            header = headers[ index]
                            output_dic[ header] = token
                        result[ input_motif].append( output_dic)

            tomtom_file.close()
        except (KeyError, IndexError, IOError),  io_exce:
            Log.log( "MotifProcessor.parseTOMTOMResult : Unable to open/parse the TOMTOM result : '" + file_path + "'. From:\n\t---> " +str( io_exce))
            return None
            #raise ExecutionException( "MotifProcessor.parseTOMTOMResult : Unable to open the TOMTOM result : '" + file_path + "'. From \n\t" +str( io_exce))
 
        return result



    # ---------------------------------------------------------------------------------------------
    # Analyse the compare-matrices result and store the information in the final result
    def analyseTOMTOMResult( self, output_commstruct, result, motif_ids):
        
        
        # define some constants for the analysis
        offset_header     = "Optimal offset"
        id2_header        = "Target ID"
        consensus2_header = "Target consensus"
        
        motif_list = []
        
        # parse the bedseq and the bedseq alignment to get the motifs
        for bedseq in output_commstruct.bedToMA.keys():
            for alignment in output_commstruct.bedToMA[ bedseq]:
                if alignment.motifs != None or len( alignment.motifs) > 0:
                    analysis_result = {}
                    for input_motif in alignment.motifs:
                        # if the motif name is in the results, retrieve the list of result series
                        if input_motif.name in result.keys():
                            motif_dics = result[ input_motif.name]
                            #for each entry in the result series, create an output motif
                            for motif_dic in motif_dics:
                                try:
                                    output_motif_offset = int( motif_dic[ offset_header])
                                    output_motif_name = motif_dic[ id2_header]
                                    consensus2 = motif_dic[ consensus2_header]
                                    output_motif_length = len( consensus2)
                                    # Here the offset is substract because it seems it is the offset of the input motif respect to the database one
                                    output_motif = Motif( input_motif.indexStart - output_motif_offset, input_motif.indexStart - output_motif_offset + output_motif_length, output_motif_name, None )
                                    output_motif.consensus = consensus2
                                    output_motif.offset = - output_motif_offset
                                    if output_motif_name in  motif_ids.keys():
                                        output_motif.id = motif_ids[ output_motif_name]
                                    if not input_motif in analysis_result.keys():
                                        analysis_result[ input_motif] = []
                                    analysis_result[ input_motif].append( output_motif)
                                    motif_list.append( output_motif)
                                except ValueError, val_exce:
                                    Log.log( "MotifProcessor.analyseTOMTOMResult : The offset or length value of a motif is not an integer : " + str( motif_dic) + ". From:\n\t--->" + str( val_exce))
                                    
                    for input_motif in analysis_result.keys():
                        alignment.motifs.remove( input_motif)
                        alignment.motifs.extend( analysis_result[ input_motif])

        return motif_list



    # #######################################
    #
    #     METHODS COMMON FUNCTIONS
    #
    # #######################################
    

    # ---------------------------------------------------------------------------------------------
    # Generates the list of unrecognized motifs using first the possible conserved blocks
    # If no conserved block are found, takes the whole peaks as motif
    def getMotifList(self, output_commstruct, desired_species_list):
        
        motif_list = []
        motif_size_list = []
        bedseq_size_list = []
        min_size = 1000000
        max_size = -1
        total_size = 0
        for bedseq in output_commstruct.bedToMA.keys():
            bedseq_size_list.append( bedseq.indexEnd - bedseq.indexStart)
            for alignment in output_commstruct.bedToMA[ bedseq]:
                if alignment.motifs != None or len( alignment.motifs) > 0:
                    for motif in alignment.motifs:
                        motif_list.append( motif)
                        motif_length = motif.indexEnd - motif.indexStart
                        motif_size_list.append( motif_length)
                        if motif_length < min_size:
                            min_size = motif_length
                        if motif_length > max_size:
                            max_size = motif_length
                        total_size += motif_length
                        

        Log.trace( "MotifProcessor.getMotifList : Conserved regions found : " + str( len( motif_list)))
        
        if len( motif_list) == 0:
            for bedseq in output_commstruct.bedToMA.keys():
                for alignment in output_commstruct.bedToMA[ bedseq]:
                    pwm = PWM()
                    pwm.initFromAlignment( alignment, desired_species_list)
                    new_motif = Motif( bedseq.indexStart, bedseq.indexEnd, "", pwm)
                    new_motif.composeName( bedseq.name)
                    motif_list.append( new_motif)
                    motif_length = new_motif.indexEnd - new_motif.indexStart
                    motif_size_list.append( motif_length)
                    if motif_length < min_size:
                        min_size = motif_length
                    if motif_length > max_size:
                        max_size = motif_length
                    total_size += motif_length
                    
            Log.trace( "MotifProcessor.getMotifList : Using complete sequences as conserved regions : " + str( len( motif_list)))

        mean_size = (int) (total_size / float( len( motif_list)))

        #Memorize the statistics
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_NUMBER] = len( motif_list)
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_MIN_SIZE] = min_size
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_MAX_SIZE] = max_size
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_MEAN_SIZE] = mean_size
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_TOTAL_SIZE] = total_size
        
        # Output the histogram
        file_infos = RSATUtils.outputHistogram( motif_size_list, 10, self.out_path, "conservedRegionSize", self.component.pipelineName, "", "Conserved block size", "Number of blocks", ('5', '6'), False)
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_SIZE_PATH] = file_infos[0]
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_SIZE_GRAPH_PATH] = file_infos[1]
        
        return ( motif_list, motif_size_list, bedseq_size_list)



    # --------------------------------------------------------------------------------------
    # Create and start a new Thread to parse the next motif pack in the queue
    def startNewThread(self, count_thread, method, arguments, motif_pack_queue, output_commstruct, out_path):
        
        if not motif_pack_queue.empty():
            # get the motif pack to analyse
            motif_pack = motif_pack_queue.get()
            # create the thread according to the method used
            if method == MotifProcessor.METHOD_VALUE_RSAT:
                my_thread = threading.Thread( None, self.launchCompareMatrices, "MotifPack_" + str( count_thread), ( arguments, motif_pack, output_commstruct, out_path))
            elif method == MotifProcessor.METHOD_VALUE_TOMTOM:
                my_thread = threading.Thread( None, self.launchTOMTOM, "MotifPack_" + str( count_thread), ( arguments, motif_pack, output_commstruct, out_path))
            else:
                raise ExecutionException( "MotifProcessor.startNewThread : Unable to start Thread. Provided method is unknown : " + method)
            
            # star the thread
            my_thread.start()
            
            return my_thread
        else:
            return None
        
            


    # ---------------------------------------------------------------------------------------------
    # Retrieve the PWM of given motifs in given database file
    def retrieveMotifsMatrices(self, motif_list, database_file_path, database_format):
        
        if database_format == "transfac" or database_format == "tf":
            MotifUtils.getMotifsPWMFromJasparTF( motif_list, database_file_path)
        elif database_format == "meme":
            MotifUtils.getMotifsPWMFromMeme( motif_list, database_file_path)
        else:
            Log.trace( "MotifProcessor : Unknown database format : " + database_format)


    # ---------------------------------------------------------------------------------------------
    # Remove the dots at the extremities of the token and return also the number of dots removedirs
    # at the beginning
    def removeNoInfoChar(self, token):
        
        result = ""
        if token != None and len( token) > 0:
            left_index = 0
            for index in range( len( token)): 
                if token[ index] != ".":
                    left_index = index
                    break
            for index in range( len( token)):
                if token[ -index-1] != ".":
                    if index == 0:
                        right_index = len( token)
                    else:
                        right_index = -index
                    break
            
            result = token[ left_index:right_index]
            
        return result


    # ---------------------------------------------------------------------------------------------
    # Try to wait until the file at the given path exists and has a constant size
    def verifyResultFile(self, path):
        
        tries_nb = 0
        file_exists = False
        while tries_nb < 5 and file_exists == False:
            if not os.path.exists( path):
                Log.log("MotifProcessor.verifyResultFile : " + threading.currentThread().getName() + " : result file '" + path + "' is not accessible. Retries = " + str( tries_nb))
                time.sleep( 1 + tries_nb)
                tries_nb += 1
            else:
                file_exists = True

        if file_exists == True:
            old_file_size = FileUtils.getSize( path)
            time.sleep( 1)
            tries_nb = 0
            constant_size = False
            while tries_nb < 5 and constant_size == False:
                file_size = FileUtils.getSize( path)
                if old_file_size != file_size:
                    Log.log("MotifProcessor.verifyResultFile : " + threading.currentThread().getName() + " : result file size changed. Retries = " + str( tries_nb))
                    tries_nb += 1
                    time.sleep( 1 + tries_nb)
                else:
                    constant_size = True
                    
            if constant_size != True:
                Log.log( "MotifProcessor: verifyResultFile : Unable to get the result file at constant size after 5 retries at :" + path)
                return False
            
        else:
            Log.log( "MotifProcessor: verifyResultFile : Unable to get access to the result file after 5 retries at :" + path)
            file_list = FileUtils.getFileList( os.path.dirname( path), "tab", False)
            for file_path in file_list:
                Log.log( "MotifProcessor: verifyResultFile : tab file in :" + file_path)
            return False
        
        return True
