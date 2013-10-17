
import sys, os, shutil, gc

from manager.ProgressionManager import ProgressionManager

from utils.exception.ExecutionException import ExecutionException
from utils.exception.ParsingException import ParsingException
from utils.log.Log import Log
from utils.FileUtils import FileUtils

from processor.ProcessorFactory import ProcessorFactory


class Component:
    
    INPUT_FILE_PARAM = "InputFile"
    
    # --------------------------------------------------------------------------------------
    def __init__( self, processor_name, rank, branch):
        
        self.pipelineName = ""
        self.rank = str( rank)
        self.branch = branch
        self.processorName = processor_name
        self.previousComponents = []
        self.nextComponents = []
        self.parameters = {}
        self.runtimeParameters = {}
        self.executed = False
        self.resumed = False
        self.resultClass = None
        self.outputDir = ""
        
        # compute processor short name
        invert_proc_name = self.processorName[::-1]
        cut_invert = invert_proc_name[:invert_proc_name.index( ".")]
        self.processorShortName =  cut_invert[::-1]
        
        # Get the processor display name
        processor_class_info = ProcessorFactory.decomposeProcessorName( self.processorName)
        __import__( processor_class_info[0])
        self.processorDisplayName = getattr( sys.modules[ processor_class_info[0]], processor_class_info[1]).getDisplayName()
    
    
    # --------------------------------------------------------------------------------------
    def addParameters(self, key, value):
        
        if key != None and value!= None:
            self.parameters[ key] = value


    # --------------------------------------------------------------------------------------
    # Add a component as a previous component in the execution order
    def addPreviousComponents(self, components):
        
        if components != None:
            for component in components:
                if not component in self.previousComponents:
                    self.previousComponents.append( component)


    # --------------------------------------------------------------------------------------
    # Add a component as a next component in the execution order
    def addNextComponent(self, component):
        
        if not component in self.nextComponents:
            self.nextComponents.append( component)


    # --------------------------------------------------------------------------------------
    # Return True if all the previous components have been executed or have been resumed
    def canStart(self):
        
        can_execute = True
        
        if len( self.previousComponents) == 0:
            return can_execute
            
        for component in self.previousComponents:
            if not component.executed and not component.resumed:
                can_execute = False
        
        return can_execute


    # --------------------------------------------------------------------------------------
    # Return True if all the previous components have been resumed (and none have been executed)
    def canResume(self):
        
        can_resume = True
        
        if len( self.previousComponents) == 0:
            return can_resume
            
        for component in self.previousComponents:
            if component.executed == True or component.resumed == False:
                can_resume = False
        
        return can_resume



    # --------------------------------------------------------------------------------------
    # If resume is True, verify if the output data of the associate processor can be retrieved from
    # previous runs. If so, the class of the output CommStruct is retrieved and the Component is declared as 'resumed'
    # If not, the processor is run and the component declared as 'executed'
    def start( self, pipeline, pipeline_out, runtime_params, resume = False):

        self.outputDir = pipeline_out
        self.runtimeParameters = runtime_params
        
        if resume == True:
            # Test if the previous component were all resumed
            if self.canResume():
                self.resumed = False
                # test if the Component parameters have changed since the previous run. If so, the processor cannot
                # be resumed and must be re-run
                if self.verifyConfig():
                    # Test if an output file of a previous run of the associated processor can be retrieved
                    # If so (or if the processor has output no files), the Component is declared as resumed and returns True
                    try:
                        output_filepath = self.getOutputFilePath()
                        if os.path.isfile( output_filepath):
                            authorized_output_classes = self.getAuthorizedOutputClasses()
                            if authorized_output_classes != None:
                                for output_class in authorized_output_classes:
                                    try:
                                        output_commstruct = output_class.fromXMLFile( output_filepath)
                                        if output_commstruct != None:
                                            self.resultClass = output_class
                                            self.resumed = True
                                            self.executed = False
                                            ProgressionManager.setComponentStatus( self, ProgressionManager.RESUMED_STATUS)
                                            Log.trace( "Component.execute : Resuming data from file : " + output_filepath)
                                            output_commstruct = None
                                            gc.collect()
                                            return True
                                    except BaseException, exce:
                                        Log.info( "Component.execute : Tried to resume output file with class '" + str( output_class) + "' : "+ str( exce))
                                        pass
                            else:
                                self.resumed = True
                                self.executed = False
                                ProgressionManager.setComponentStatus( self, ProgressionManager.RESUMED_STATUS)
                                return True
                    except IOError, io_exce:
                        Log.trace ("Component.execute : Unable to open output file to resume processor '" + self.processorName + "'. From\n\t---> " + str( io_exce))
                    
                    # Here, the processor cannot be resumed, for any reason linked to outfiles, 
                    Log.trace ("Component.execute : No output file found for processor '" + self.processorName + "': executing it")                
                    
                self.removePreviousOutputs()

            # If the processor does not have to be resumed because previous components were not resumed,
            # removes all old output files and the processor is executed
            else:
                Log.trace( "Component.execute : Processor '" + self.processorName + "' cannot be resumed since previous components have been executed.")
                self.removePreviousOutputs()

        
        return self.executeProcessor( pipeline)




    # --------------------------------------------------------------------------------------
    # Instanciate and execute the processor associated to the component, 
    # store the config and the result and declares the Component as 'executed'
    def executeProcessor( self, pipeline):
        
        ProgressionManager.setComponentStatus( self, ProgressionManager.RUNNING_STATUS)
        
        # Save the config to file
        self.outputConfig()
        
        #Retrieve the list of authorized inputs for the associated processor
        input_commstructs = self.getInputCommStructs()
            
        # Creates an instance of associated processor and execute it
        processor_instance = ProcessorFactory.getProcessorInstance( self.processorName, self)
        self.processorDisplayName = processor_instance.getDisplayName()
        processor_instance.setPipeline( pipeline)
        output_commstruct = processor_instance.execute( input_commstructs)
        
        # Manage the result and generate the XML output
        self.resultClass = output_commstruct.__class__
        authorized_output_classes = self.getAuthorizedOutputClasses()
        if output_commstruct != None:
            #print "\nRESULTS OF PROCESSOR '" + self.processorName +" ' : \n\n" + output_commstruct.toString()
            
            if authorized_output_classes != None:
                if self.resultClass in authorized_output_classes:
                    output_commstruct.setOrigin( self)
                    ProgressionManager.setTaskProgression( "Writing result to file", self, 0.0)
                    self.outputResultToXml( output_commstruct)
                    ProgressionManager.setTaskProgression( "Writing result to file", self, 1.0)
                    self.executed = True
                    self.resumed = False
                else:
                    raise ExecutionException( "Component.executeProcessor : The processor '" + self.processorName + "' returned an output of class '" + str( self.resultClass) + " instead of one of '" + str( authorized_output_classes) + "'")
            else:
                raise ExecutionException( "Component.executeProcessor : The processor '" + self.processorName + "' returned an output of class '" + str( self.resultClass) + " while it should not generate outputs")
        else:
            if authorized_output_classes != None:
                raise ExecutionException( "Component.executeProcessor : The processor '" + self.processorName + "' returned no output but should generate an output of class in '" + str( authorized_output_classes) + "'")
            else:
                self.executed = True
                self.resumed = False
        
        ProgressionManager.setComponentStatus( self, ProgressionManager.EXECUTED_STATUS)
        
        input_commstructs = None
        output_commstruct = None
        gc.collect()
        
        return self.executed
        
       


    # --------------------------------------------------------------------------------------
    # Return the list of Commstruct outputed from previous components
    def getInputCommStructs(self):
        
        authorized_input_classes = self.getAuthorizedInputClasses()
        
        input_commstructs = []
        if authorized_input_classes != None:
            input_file = self.getParameter( Component.INPUT_FILE_PARAM, False)
            if input_file == None:
                #Compares the list of authorized inputs to outputs of previous components
                for component in self.previousComponents:
                    previous_result_class = component.resultClass
                    if previous_result_class in authorized_input_classes:
                        input_commstruct = previous_result_class.fromXMLFile( component.getOutputFilePath())
                        if input_commstruct != None:
                            input_commstructs.append( input_commstruct)
                    else:
                        raise ExecutionException( "Component.getInputCommStructs : input is not of the right class. Class is '" + previous_result_class + "' but waited classes are " + str( authorized_input_classes))
            else:
                #Try to read the input file using classes authorized as input
                for input_class in authorized_input_classes:
                    try:
                        Log.trace( "Component.getInputCommStructs : Trying to load data from file : " + input_file)
                        input_commstruct = input_class.fromXMLFile( input_file)
                        if input_commstruct != None:
                            input_commstructs.append( input_commstruct)
                        Log.trace( "Component.getInputCommStructs : Data correctly loaded")
                    except Exception, exce:
                        Log.trace( "Component.getInputCommStructs : Data not loaded using class '" + str( input_class) + "' : "+ str( exce))
                        pass
                if len( input_commstructs) == 0:
                    raise ExecutionException( "Component.getInputCommStructs : The provided input file does not contain information the processor '" + self.processorName + "' can manage : " + input_file)
                

        return input_commstructs
        


    # --------------------------------------------------------------------------------------
    # Return the list of classes component processor can manage as input
    def getAuthorizedInputClasses(self ):
        
        processor_class_info = ProcessorFactory.decomposeProcessorName( self.processorName)
        __import__( processor_class_info[0])
        authorized_input_classes = getattr( sys.modules[ processor_class_info[0]], processor_class_info[1]).getInputCommStructClass()
        
        return authorized_input_classes



    # --------------------------------------------------------------------------------------
    # Return the list of classes component processor manages as output
    def getAuthorizedOutputClasses(self ):
        
        processor_class_info = ProcessorFactory.decomposeProcessorName( self.processorName)
        __import__( processor_class_info[0])
        authorized_output_classes = getattr( sys.modules[ processor_class_info[0]], processor_class_info[1]).getOutputCommStructClass()
        
        return authorized_output_classes



    # --------------------------------------------------------------------------------------
    # Saves the output CommStruct to an suitable XML file
    def outputResultToXml(self, output_commstruct):
        
        output_commstruct.toXMLFile( self.getOutputFilePath())



    # --------------------------------------------------------------------------------------
    # Check if the config from previous run and the current config math
    def verifyConfig(self):
        
        previous_config = self.getConfigFromFile()
        
        if previous_config == None:
            Log.trace("Component.verifyConfig : No config file found for processor '" + self.processorName + "'. Executing it")
            return False
        
        same = False
        if len( previous_config) == len( self.parameters):
            if len( previous_config) == 0:
                same = True
            else:
                for param_name in previous_config.keys():
                    if param_name in self.parameters.keys():
                        if previous_config[ param_name] == self.parameters[ param_name]:
                            same = True
                        else:
                            same = False
                            break
                    else:
                        same = False
                        break
        
        if not same:
            Log.trace( "Component.verifyConfig : Configuration of processor '" + self.processorName + "' has changed since previous run: Previous = " + str( previous_config) + " while current = " + str( self.parameters)+ ". Executing processor")
        
        return same


    # --------------------------------------------------------------------------------------
    # Saves the Component parameters to text file
    def outputConfig(self):
                
        try:
            output_path = self.getConfigFilePath()
            config_file = FileUtils.openFile( output_path, "w")
            for param in self.parameters.keys():
                config_file.write( param + "=" + self.parameters[ param] + "\n")
                config_file.flush()
            config_file.close()
        except IOError, io_exce:
            raise ExecutionException( "Component.outputConfig : Unable to write component config in file '" + output_path + "'. From:\n\t---> " + str( io_exce))



    # --------------------------------------------------------------------------------------
    # Reads the Component parameters from text file
    def getConfigFromFile(self):
        
        config = {}
        
        try:
            output_path = self.getConfigFilePath()
            config_file = FileUtils.openFile( output_path)
            for line in config_file:
                tokens = line.split( "=")
                if tokens != None and len( tokens) == 2:
                    if tokens[1][-1] == "\n":
                        value = tokens[1][:-1]
                    else:
                        value = tokens[1]
                    config[ tokens[0]] = value
                else:
                    raise ParsingException( "Component.getConfigFromFile : Wrongly formatted config file. Should have '<param_name> = <param_value>' instead of " + line)
            config_file.close()
        except IOError:
            return None

        return config



    # --------------------------------------------------------------------------------------
    # Provides the full path name for the file where result will be stored
    def getOutputFilePath(self):
        
        return os.path.join( self.outputDir, self.getComponentPrefix() + "_output.xml")




    # --------------------------------------------------------------------------------------
    # Provides the full path name for the file where config will be stored
    def getConfigFilePath(self):
        
        return os.path.join( self.outputDir, self.getComponentPrefix() + "_config.txt") 



    # --------------------------------------------------------------------------------------
    # Provide a prefix containing the pipeline name and the component rank, branch and processor name
    def getComponentPrefix(self):
        
        return self.pipelineName + "_" + self.rank + "_" + self.branch + self.processorShortName



    # --------------------------------------------------------------------------------------
    # Return the value of the parameter if exists
    def getParameter(self, param_name, mandatory = True):
        
        try:
            if param_name in self.parameters.keys():
                param_value =  self.parameters[ param_name]
            else:
                param_value =  self.runtimeParameters[ param_name]
            return param_value
        except KeyError, key_exce:
            if mandatory:
                raise ExecutionException( "Component.getParameter : Runtime parameter :'" + param_name + "' does not exists. From:\n\t---> " + str( key_exce))
            else:
                return None


    # --------------------------------------------------------------------------------------
    # Return the value of the parameter if exists
    def getParameterDic(self):
        
        return self.parameters

    # --------------------------------------------------------------------------------------
    # Remove the possible output files and dir remaning from previous runs
    def removePreviousOutputs(self):
        
        try:
            shutil.rmtree( os.path.join( self.outputDir, self.getComponentPrefix()), True)
            os.remove( self.getOutputFilePath())
        except ( OSError, IOError) :
            pass
 

    # --------------------------------------------------------------------------------------
    # Return a string representation of the component
    def toString( self):
        
        result = "Component " + str( self.rank) + "-" + self.branch + ":" + self.processorName + " (" + self.pipelineName + ")\n"
        result += "|--Parameters\n"
        for key in self.parameters.keys():
            result += "|--|-- " + key + " = " + self.parameters[ key] + "\n"
        result += "|--Previous Components\n"
        for component in self.previousComponents:
            result += "|--|--" + component.processorName + "\n"
        result += "|--Next Components\n"
        for component in self.nextComponents:
            result += "|--|--" + component.processorName + "\n"
        
        return result
