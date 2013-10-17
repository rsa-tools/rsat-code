

from utils.exception.ExecutionException import ExecutionException
from utils.exception.ParsingException import ParsingException

from utils.log.Log import Log

class Processor:
    
    
    # --------------------------------------------------------------------------------------    
    # This def method must be called from the def method of the inherited classes
    def __init__( self):
        
        self.component = None
        self.parameters = {}
        self.pipeline = None

    # ---------------------------------------------------------------------------------------
    # Set the pipeline the processor is part of
    def setPipeline(self, pipe):
        self.pipeline = pipe

    # --------------------------------------------------------------------------------------
    # This method must be redefined in the inherited classes
    # it must return the name of the CommStruct class used as input 
    # (None if no input CommStruct)
    @staticmethod
    def getInputCommStructClass():
        
        Log.log( "The method 'getInputCommStructClass' must be implemented at the inherited class level")
        return ("Not defined", )


    # --------------------------------------------------------------------------------------
    # This method must be redefined in the inherited classes
    # it must return the name of the CommStruct class used as output
    # (None if no output CommStruct)
    @staticmethod
    def getOutputCommStructClass():
        
        Log.log( "The method 'getOutputCommStructClass' must be implemented at the inherited class level")
        return ("Not defined", )



    # --------------------------------------------------------------------------------------
    # This method must be redefined in the inherited classes
    # it must return a name that will be used as display name in the user friendly outputs
    @staticmethod
    def getDisplayName():
        
        Log.log( "The method 'getOutputCommStructClass' must be implemented at the inherited class level")
        return Processor.__class__.__name__ + " (no display name defined)"


    # --------------------------------------------------------------------------------------
    # This method must be redefined in the inherited classes
    # it must return a list of parameters name that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        Log.log( "The method 'getRequiredParameters' must be implemented at the inherited class level")
        return None



    # --------------------------------------------------------------------------------------
    # This method must be redefined in the inherited classes
    # it must execute the processor analysis part
    def execute(self, comm_struct, pipeline):
        
        Log.log( "The method 'execute' must be implemented at the inherited class level")
        return None


    # --------------------------------------------------------------------------------------
    # Return the value of the parameter if exists
    def getParameter(self, param_name, mandatory = True):
        
        try:
            param_value =  self.parameters[ param_name]
            return param_value
        except (TypeError, KeyError), key_exce:
            if mandatory:
                raise ExecutionException( "Processor.getParameter : Processor does not know the parameter :'" + param_name + "'. From:\n\t---> " + str( key_exce))
            else:
                return None


    # --------------------------------------------------------------------------------------
    # Return the value of the parameter as int if exists
    def getParameterAsint(self,  param_name, mandatory = True):
        
        try:
            param_value = int( self.getParameter( param_name, mandatory))
            return param_value
        except (TypeError, ValueError), val_exce:
            if mandatory:
                raise ExecutionException( "Processor.getParameterAsint : Unable to convert the value of parameter :'" + param_name + "'. From:\n\t---> " + str( val_exce))
            else:
                return None


    # --------------------------------------------------------------------------------------
    # Return the value of the parameter as float if exists
    def getParameterAsfloat(self,  param_name, mandatory = True):
        
        try:
            param_value = float( self.getParameter( param_name, mandatory))
            return param_value
        except (TypeError, ValueError), val_exce:
            if mandatory:
                raise ExecutionException( "Processor.getParameterAsfloat : Unable to convert the value of parameter :'" + param_name + "'. From:\n\t---> " + str( val_exce))
            else:
                return None


    # --------------------------------------------------------------------------------------
    # Return the value of the token as int if possible and raise an exception if the value is required but not converted
    def getTokenAsint( self, token, required=True):
        
        try:
            att_value = int( token)
            return att_value
        except (TypeError, ValueError), val_exce:
            if required:
                raise ParsingException( "Processor.getTokenAsint : Unable to convert the token to int :'" + token + "'. From:\n\t---> " + str( val_exce))
            else:
                return None
                

    # --------------------------------------------------------------------------------------
    # Return the value of the token as int if possible and raise an exception if the value is required but not converted
    def getTokenAsfloat( self, token, required=True):
        
        try:
            att_value = float( token)
            return att_value
        except (TypeError, ValueError), val_exce:
            if required:
                raise ParsingException( "Processor.getTokenAsfloat : Unable to convert the token to float :'" + token + "'. From:\n\t---> " + str( val_exce))
            else:
                return None



    # --------------------------------------------------------------------------------------
    # Add indentation to the ElementTree in order to have a pretty print
    # in the XML file (used by Processor subclasses)
    def indent(self, elem, level=0):
            
        i = "\n" + level*"  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for elem in elem:
                self.indent(elem, level+1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i   
                
                
# eflag: FileType = Python2
