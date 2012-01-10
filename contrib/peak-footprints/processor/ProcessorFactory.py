
import sys

from utils.exception.ExecutionException import ExecutionException

class ProcessorFactory:

    # --------------------------------------------------------------------------------------
    @staticmethod
    def getProcessorInstance( processorName, component):
        
        processor_instance = None
        
        try:
            class_info = ProcessorFactory.decomposeProcessorName( processorName)
            module_name = class_info[0]
            class_name = class_info[1]
            
            __import__( module_name)
            
            processor_instance = getattr( sys.modules[ module_name], class_name)()
            
            processor_instance.component = component
            processor_instance.parameters = component.parameters
            
        except Exception, exce:
            raise ExecutionException( "ProcessorFactory : Unable to create processor '" + processorName + "' instance. From:\n\t---> " +  str( exce))
    
        return processor_instance


    # --------------------------------------------------------------------------------------
    # Return the module name and class name of the given processor name
    @staticmethod
    def decomposeProcessorName( processor_name):
        
        paths = processor_name.split('.')
        module_name = '.'.join( paths[:-1])
        class_name = paths[-1]
        
        if len( module_name) > 1 and class_name != None and len( class_name) > 0:
            return [ module_name, class_name]
        else:
            raise Exception( "ProcessorFactory : Unable to identify the module or class of processor '" + processor_name + "' : module=" + module_name + " class=" + class_name)
        
# eflag: FileType = Python2
