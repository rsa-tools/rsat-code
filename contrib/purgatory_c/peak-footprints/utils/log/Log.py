
import os, datetime

from utils.Constants import Constants
from utils.FileUtils import FileUtils
from utils.exception.ConfigException import ConfigException

class Log:
    
    logInstance = None

    # --------------------------------------------------------------------------------------
    def __init__(self, file_path, verbosity, option):
        
        self.verbosity = verbosity
        self.log = True
        
        if self.verbosity > 0:
            self.trace = True
        else:
            self.trace = False
            
        if self.verbosity > 1:
            self.info = True
        else:
            self.info = False
        
        try:
            if self.log:
                log_file_path = os.path.join( file_path, Constants.LOG_FILE_NAME)
                self.logFile = FileUtils.openFile( log_file_path, option, 0666)
                
            if self.trace:
                trace_file_path = os.path.join( file_path, Constants.TRACE_FILE_NAME)
                self.traceFile = FileUtils.openFile( trace_file_path, option, 0666)
        except Exception, exce:
            raise ConfigException( "Log.__init__ : Unable to create log file in directory '" + file_path + "'. From:\n\t---> " + str( exce))



    # --------------------------------------------------------------------------------------
    @staticmethod
    def initLog( file_path, verbosity = 0):
        
        if Log.logInstance != None:
            verbosity = Log.logInstance.verbosity
            
        Log.logInstance = Log( file_path, verbosity, "w")



    # --------------------------------------------------------------------------------------
    @staticmethod
    def switchFiles( file_path, verbosity = 0):
        
        if Log.isInitialized():
            Log.closeFiles()
            if os.path.isfile( os.path.join( file_path, Constants.LOG_FILE_NAME)):
                option = "a"
            else:
                option = "w"
            
            Log.logInstance = Log( file_path, Log.logInstance.verbosity, option)
        
        else:
            Log.initLog( file_path, verbosity)
        

    # --------------------------------------------------------------------------------------
    @staticmethod
    def closeFiles():
        
        if Log.logInstance.log:
            Log.logInstance.logFile.close()
        if Log.logInstance.trace:
            Log.logInstance.traceFile.close()


    # --------------------------------------------------------------------------------------
    @staticmethod
    def log( message):
        
        if Log.logInstance.log == True:
            timestamp = str( datetime.datetime.now())
            print "LOG : " + timestamp + " : " + message
            Log.logInstance.logFile.write( "LOG : " + timestamp + " : " + message + "\n")
            Log.logInstance.logFile.flush()


    # --------------------------------------------------------------------------------------
    @staticmethod
    def trace( message):
        
        if Log.logInstance.trace == True:
            timestamp = str( datetime.datetime.now())
            print "TRACE : " + timestamp + " : " + message
            Log.logInstance.traceFile.write( "TRACE : " + timestamp + " : " + message + "\n")
            Log.logInstance.traceFile.flush()


    # --------------------------------------------------------------------------------------
    @staticmethod
    def info( message):
        
        if Log.logInstance.info == True:
            timestamp = str( datetime.datetime.now())
            print "INFO : " + timestamp + " : " + message
            Log.logInstance.traceFile.write( "INFO : " + timestamp + " : " + message + "\n")
            Log.logInstance.traceFile.flush()            


    # --------------------------------------------------------------------------------------
    @staticmethod
    def isInitialized():
        
        return Log.logInstance != None

    # --------------------------------------------------------------------------------------
    @staticmethod
    def isLogEmpty():
        
        if Log.isInitialized():
            size = FileUtils.getSize( Log.logInstance.logFile.name)
            if size != 0:
                return False

        return True

# eflag: FileType = Python2
