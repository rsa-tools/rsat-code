
import os, datetime

from utils.Constants import Constants
from utils.FileUtils import FileUtils
from utils.exception.ConfigException import ConfigException

class ListenerLog:
    
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
                self.logFile = FileUtils.openFile( os.path.join( file_path, "listener_" + Constants.LOG_FILE_NAME), option)
                
            if self.trace:
                self.traceFile = FileUtils.openFile( os.path.join( file_path, "listener_" + Constants.TRACE_FILE_NAME), option)
        except Exception, exce:
            raise ConfigException( "Log.__init__ : Unable to create log file in directory '" + file_path + "'. From:\n\t---> " + str( exce))



    # --------------------------------------------------------------------------------------
    @staticmethod
    def initLog( file_path, verbosity = 0):
        
        if ListenerLog.logInstance != None:
            verbosity = ListenerLog.logInstance.verbosity
            
        ListenerLog.logInstance = ListenerLog( file_path, verbosity, "w")



    # --------------------------------------------------------------------------------------
    @staticmethod
    def switchFiles( file_path, verbosity = 0):
        
        if ListenerLog.isInitialized():
            ListenerLog.closeFiles()
            if os.path.isfile( os.path.join( file_path, Constants.LOG_FILE_NAME)):
                option = "a"
            else:
                option = "w"
            
            ListenerLog.logInstance = ListenerLog( file_path, ListenerLog.logInstance.verbosity, option)
        
        else:
            ListenerLog.initLog( file_path, verbosity)
        

    # --------------------------------------------------------------------------------------
    @staticmethod
    def closeFiles():
        
        if ListenerLog.logInstance.log:
            ListenerLog.logInstance.logFile.close()
        if ListenerLog.logInstance.trace:
            ListenerLog.logInstance.traceFile.close()


    # --------------------------------------------------------------------------------------
    @staticmethod
    def log( message):
        
        if ListenerLog.logInstance.log == True:
            timestamp = str( datetime.datetime.now())
            print "LOG : " + timestamp + " : " + message
            ListenerLog.logInstance.logFile.write( "LOG : " + timestamp + " : " + message + "\n")
            ListenerLog.logInstance.logFile.flush()


    # --------------------------------------------------------------------------------------
    @staticmethod
    def trace( message):
        
        if ListenerLog.logInstance.trace == True:
            timestamp = str( datetime.datetime.now())
            print "TRACE : " + timestamp + " : " + message
            ListenerLog.logInstance.traceFile.write( "TRACE : " + timestamp + " : " + message + "\n")
            ListenerLog.logInstance.traceFile.flush()


    # --------------------------------------------------------------------------------------
    @staticmethod
    def info( message):
        
        if ListenerLog.logInstance.info == True:
            timestamp = str( datetime.datetime.now())
            print "INFO : " + timestamp + " : " + message
            ListenerLog.logInstance.traceFile.write( "INFO : " + timestamp + " : " + message + "\n")
            ListenerLog.logInstance.traceFile.flush()            


    # --------------------------------------------------------------------------------------
    @staticmethod
    def isInitialized():
        
        return ListenerLog.logInstance != None

    # --------------------------------------------------------------------------------------
    @staticmethod
    def isLogEmpty():
        
        if ListenerLog.isInitialized():
            size = FileUtils.getSize( ListenerLog.logInstance.logFile.name)
            if size != 0:
                return False

        return True

# eflag: FileType = Python2
