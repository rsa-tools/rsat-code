
        
import os, shutil, glob

class FileUtils:
    
    # --------------------------------------------------------------------------------------
    # Return the size of the given file or directory
    @staticmethod
    def getSize( path):
        
        return os.path.getsize( path)
    
    # --------------------------------------------------------------------------------------
    @staticmethod
    # Search for files with given extensionat the given path or into a subdirectory having the species name
    def getFileList( path, extension, subdir = None):

        result = []
        # test if the path is a file path
        if os.path.isfile( path):
            # If it is a file with given extension, return only this file. If not, return None
            if path[ -len( extension)-1:].lower() == "." + extension:
                result.append( path)
            else:
                return None
        else:
            # if the path is a directory path, test if it contains files with the given extension
            result = glob.glob(  os.path.join( path,  "*." + extension))
            if result == None or len( result) == 0:
                # if directory does not contain serached files, search in given subdirectory if any
                if subdir != None:
                    sub_path = os.path.join( path, subdir)
                    result = glob.glob(  os.path.join( sub_path,   "*." + extension))
                    if result == None or len( result) == 0:
                        return None
                else:
                    return None
                

        return result
        
        
        
    # --------------------------------------------------------------------------------------
    # Return the list of directories present in the folder at the given path
    @staticmethod
    def getDirectoryList( path):
        
        result = []
        if path != None and len( path) > 0:
            list_entries = os.listdir( path)
            if list_entries != None:
                for entry in list_entries:
                    if not os.path.isfile( os.path.join( path, entry)):
                        result.append( entry)

        return result



    # --------------------------------------------------------------------------------------
    # create a directory at the given path is possible
    @staticmethod
    def createDirectory( path):
        
        if os.path.exists( path):
            if not os.path.isdir( path):
                return False
        else:
            os.mkdir( path)
        
        return True


    # --------------------------------------------------------------------------------------
    # Remove the file at the given path
    @staticmethod
    def removeFile( path):
        
        if os.path.isfile( path):
            os.remove( path)



    # --------------------------------------------------------------------------------------
    # Remove the file at the given path
    @staticmethod
    def copyFile( file_path, destination_path):
        
        if os.path.isfile( file_path):
            shutil.copy( file_path, destination_path)


    # --------------------------------------------------------------------------------------
    # Return the relative path of the given path respect to the given directory
    @staticmethod
    def getRelativePath( path, directory) :
        
        return os.path.relpath( path, directory)
        