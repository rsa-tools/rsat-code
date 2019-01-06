
from utils.FileUtils import FileUtils
from utils.Constants import Constants
from utils.log.ListenerLog import ListenerLog
from utils.exception.ParsingException import ParsingException

import time

class PipelineListener:
    
    
    # --------------------------------------------------------------------------------------
    def __init__(self, pipeline_manager, dir_path):
        
        self.pipelineManager = pipeline_manager
        self.dirPath = dir_path
    
    
    # --------------------------------------------------------------------------------------
    # Create the listener instance and start it
    @staticmethod
    def start( pipeline_manager, dir_path):
        
        FileUtils.createDirectory( dir_path, 0777)
        listener = PipelineListener( pipeline_manager, dir_path)
        listener.run()
        
    
    # --------------------------------------------------------------------------------------
    # Look permanently to the reference directory to check for new files
    # When files are found, they are parsed to put in queue a new pipeline command
    # and the files are removed
    def run(self):
        
        while True:

            file_pathes = FileUtils.getFileList( self.dirPath, "ticket")
            if file_pathes != None and len( file_pathes) > 0:
                for path in file_pathes:
                    command_params = self.readCommand( path)
                    if command_params[0] != None:
                        pipelines_filepath = command_params[0]
                        verbosity = self.getIntValue( command_params[1])
                        resume = command_params[2]
                        working_dir = command_params[3]
                        ListenerLog.trace( "PipelineListener.run : Adding command to queue : " + str( command_params))
                        self.pipelineManager.addToQueue( pipelines_filepath, None, verbosity, resume, working_dir)
                    FileUtils.removeFile( path)

            time.sleep(5)

    
    # --------------------------------------------------------------------------------------
    # read the file at the given path to get the command options
    def readCommand(self, path):
        
        print" READING : " + path
        
        result = [None, 0, "True", None]
        tries = 0
        while tries <= 1:
            try:
                file = FileUtils.openFile( path)
                while True:
                    line = file.readline()
                    if len( line) == 0:
                        file.close()
                        return result
                    elif not line.isspace() and line[0] != Constants.COMMENT_CHAR:
                        tokens = line.split()
                        if len( tokens) > 0 and len( tokens) <= 4:
                            for index in range( len( tokens)):
                                result[ index] = tokens[ index]
                file.close()
            except IOError, io_exce:
                print "ioerror : " + str( io_exce)
                ListenerLog.log( "PipelineListener.readCommand : Unable to read file : " + path + ". From:\n\t---> " + str( io_exce))
                tries += 1
                if tries <= 1:
                    ListenerLog.log( "PipelineListener.readCommand : Retrying in 5 seconds")
                    time.sleep( 5)
                
            except ParsingException, par_exce:
                print "ParsingException : " + str( par_exce)
                ListenerLog.log( "PipelineListener.readCommand : Unable to read file : " + path + ". From:\n\t---> " + str( par_exce))
                tries += 1
                if tries <= 1:
                    ListenerLog.log( "PipelineListener.readCommand : Retrying in 5 seconds")
                    time.sleep( 5)

        if result[0] == None:
            ListenerLog.log( "PipelineListener.readCommand : File is not correctly formated : " + path + ". Removing it")
        
        print "Result = " + str( result)
        
        return result



    # --------------------------------------------------------------------------------------
    # Return the int value of the given token
    def getIntValue(self, token):
        
        try:
            return int( token)
        except ValueError,  val_exce:
            raise ParsingException( "PipelineListener.getIntValue : Unable to get integer value of '" + token + "'. From:\n\t---> " + str( val_exce))
