
import os, shutil, commands

from processor.Processor import Processor
from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct

from utils.log.Log import Log
from utils.MotifUtils import MotifUtils
from utils.Constants import Constants
from utils.exception.ExecutionException import ExecutionException

from manager.ProgressionManager import ProgressionManager

from common.Motif import Motif

class CompareIdentifiedMotifsProcessor( Processor):
    
    
    MOTIF_LIST_PARAM = "MotifList"
    MOTIF_DATABASE_FILE_PARAM = "MotifDatabaseFile"
    MOTIF_DATABASE_FORMAT_PARAM = "MotifDatabaseFormat"
    
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
        
        return "Comparison of identified motifs"
        


    #------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return ( CompareIdentifiedMotifsProcessor.MOTIF_DATABASE_FILE_PARAM, CompareIdentifiedMotifsProcessor.MOTIF_DATABASE_FORMAT_PARAM, CompareIdentifiedMotifsProcessor.MOTIF_LIST_PARAM)


    
    # ---------------------------------------------------------------------------------------------
    # Execute the processor
    def execute( self, input_commstructs):
        
        if input_commstructs == None or len( input_commstructs) == 0:
            raise ExecutionException( "CompareIdentifiedMotifsProcessor.execute : No inputs")
        
        input_commstruct = input_commstructs[0]

        # Retrieve the processor parameters
        database_file = self.getParameter( CompareIdentifiedMotifsProcessor.MOTIF_DATABASE_FILE_PARAM)

        database_format = self.getParameter( CompareIdentifiedMotifsProcessor.MOTIF_DATABASE_FORMAT_PARAM)
        
        motif_list_line = self.getParameter( CompareIdentifiedMotifsProcessor.MOTIF_LIST_PARAM)
        motif_name_list = motif_list_line.split()
        
        # Retrieve the PWM of the reference motifs
        reference_motif_list = self.getMotifMatrices( motif_name_list, database_file, database_format)
        
        # Retrieve the list of identified motifs
        identified_motifs = self.getIdentifiedMotifs( input_commstruct)

        # Compare motifs
        self.compareMotifs( reference_motif_list, identified_motifs)
        
        return input_commstruct
        

    # ---------------------------------------------------------------------------------------------
    # Retrieve the definition of given motifs in given database file
    def getMotifMatrices( self, motif_name_list, database_file_path, database_format):
        
        motif_list = []
        for name in motif_name_list:
            motif_list.append( Motif( 0, 0, name, None))
        
        if database_format == "transfac" or database_format == "tf":
           MotifUtils.getMotifsPWMFromJasparTF( motif_list, database_file_path)
        elif database_format == "meme":
           MotifUtils.getMotifsPWMFromMeme( motif_list, database_file_path)
        else:
            raise ExecutionException( "CompareIdentifiedMotifsProcessor.getMotifMatrices : Unknown database format : " + database_format)
       
        return motif_list
        
    
    
    # ---------------------------------------------------------------------------------------------
    # Retrieve the list of identified motifs avoiding repetition
    def getIdentifiedMotifs(self, input_commstruct):
        
        identified_motifs = {}
        for bedseq in input_commstruct.bedToMA.keys():
            for msa in input_commstruct.bedToMA[ bedseq]:
                for motif in msa.motifs:
                    if not motif.name in identified_motifs.keys():
                        identified_motifs[ motif.name] = []
                        identified_motifs[ motif.name].append( motif)
                    else:
                        exists = False
                        for other_motif in identified_motifs[ motif.name]:
                            if motif.pwm.equals( other_motif.pwm):
                                exists = True
                        if not exists:
                            identified_motifs[ motif.name].append( motif)    
        
        return identified_motifs

    # ---------------------------------------------------------------------------------------------
    # Retrieve the list of identified motifs avoiding repetition
    def compareMotifs(self, reference_motifs, identified_motifs):
        
        # Retrieve required parameters
        RSAT_PATH = self.component.getParameter( Constants.RSAT_DIR_PARAM)
        
        #Prepare outputdir
        dir_path = os.path.join( self.component.outputDir, self.component.getComponentPrefix())
        shutil.rmtree( dir_path, True)
        os.mkdir( dir_path)
        old_working_dir = os.getcwd()
        os.chdir( dir_path)
        
        # Establish the progression
        total_length = len( reference_motifs) * len( identified_motifs.keys())
        ProgressionManager.setTaskProgression( "Comparing motifs",  self.component, 0.0)

        progress = 0
        for reference_motif in reference_motifs:
            ref_file_info = self.outputMotifToTransfacFile( reference_motif, dir_path)
            for identified_motif_name in identified_motifs.keys():
                progress += 1
                if reference_motif.name != identified_motif_name:
                    count = 0
                    identified_motif_list = identified_motifs[ identified_motif_name]
                    for identified_motif in identified_motif_list:
                        count += 1
                        ident_file_info = self.outputMotifToTransfacFile( identified_motif, dir_path)
                    
                        # Compose the compare-matrices command line with all required options
                        cmd = os.path.join( RSAT_PATH , "perl-scripts/compare-matrices")
                        cmd += " -file1 " + ref_file_info[1] + " -format1 tf"
                        cmd += " -file2 " + ident_file_info[1] + " -format2 tf"
                        cmd += " -mode matches"
                        cmd += " -return all"
                        if len( identified_motif_list) > 1:
                            cmd += " -o " + reference_motif.name + "_" + identified_motif_name + "_" + count
                            count += 1
                        else:
                            cmd += " -o " + reference_motif.name + "_" + identified_motif_name
                        
                        # Execute the command
                        cmd_result = commands.getstatusoutput( cmd)
                        if cmd_result[0] != 0:
                            Log.log( "CompareIdentifiedMotifsProcessor.compareMotifs : status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
                            Log.log( "CompareIdentifiedMotifsProcessor.compareMotifs : command output is = \n" + str( cmd_result[1]))
                            continue
                
                if progress%10 == 0:
                    ProgressionManager.setTaskProgression( "Identifying motifs",  self.component, progress / float( total_length))

        # returns to initial working dir
        os.chdir( old_working_dir)
        
        
        

    # ---------------------------------------------------------------------------------------------
    # Output the motif to file and prepare the output directory
    def outputMotifToTransfacFile(self, motif, dir_path):
        
        try:
            file_path = os.path.join( dir_path, motif.name + ".tf")
            file = open( file_path, "w")
            file.write( "AC\t" + motif.name + "\n")
            file.write( "XX\n")
            file.write( "ID\t" + motif.name + "\n")
            file.write( "XX\n")
            file.write( motif.pwm.convertToTransfac())
            file.write( "XX\n")
            file.write("//\n")
            file.flush()
            file.close()
            return (dir_path, file_path)
        except IOError, io_exce:
            raise ExecutionException( "CompareIdentifiedMotifsProcessor.outputMotifListToTransfacFile : Unable to save motifs to tab file : '" + file_path + "'. From:\n\t---> " +str( io_exce))
