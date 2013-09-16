
from utils.exception.ParsingException import ParsingException

from xml.etree.ElementTree import parse
from xml.etree.ElementTree import Element, ElementTree

from Pipeline import Pipeline
from Component import Component

from utils.log.Log import Log
from utils.FileUtils import FileUtils

import os

class PipelineXMLParser:
    
    RANK = 0


    # --------------------------------------------------------------------------------------
    # Parse the given XML and create the described pipelines
    @staticmethod
    def getPipelines( pipelines_filepath):
        
        PipelineXMLParser.RANK = 0
        
        file = None
        root_element = None
        
        try:
            file = FileUtils.openFile( pipelines_filepath)
            tree = parse( file)
            root_element = tree.getroot()
            file.close()
        except IOError,  io_exce:
            raise ParsingException( "PipelineXMLParser.getPipelines : unable to open/close XML file '" + pipelines_filepath + "'. From:\n\t---> " + str( io_exce))
        
        pipelines = []
        
        for node_pipeline in root_element:
            if node_pipeline.tag.lower() == PipelineXMLParser.PIPELINE_TAG:
                pipeline = PipelineXMLParser.getPipeline( node_pipeline)
                if pipeline != None:
                    pipelines.append( pipeline)
                            
        for pipeline in pipelines:
            for component in pipeline.componentList:
                if len( component.previousComponents) == 0:
                    pipeline.addFirstComponent( component)
                if len( component.nextComponents) == 0:
                    pipeline.addLastComponent( component)
        
        return pipelines


    # --------------------------------------------------------------------------------------
    # Verify the pipeline is correctly defined
    # if so, create a Pipeline and add it to the given Pipelines
    @staticmethod
    def getPipeline( node_pipeline):

        name = PipelineXMLParser.getAttribute( node_pipeline, PipelineXMLParser.PIPELINE_NAME_ATT)
        if name != None and len( name) > 0:
            pipeline = Pipeline()
            pipeline.name = name
            previous_components = []
            PipelineXMLParser.analyseNode( node_pipeline, previous_components, pipeline, "")
            return pipeline
        else:
            raise ParsingException ( "PipelineXMLParser.getPipeline : Malformed pipeline - unable to retrieve pipeline name")


    # --------------------------------------------------------------------------------------
    # Analyse the content of the given node
    @staticmethod
    def analyseNode( parent_node, initial_previous_components, pipeline, initial_prefix):
        
        previous_components = initial_previous_components
        final_components = []
        
        branch_counter = 0
        for node in parent_node:
            if node.tag.lower() == PipelineXMLParser.COMPONENT_TAG:
                prefix = initial_prefix
                final_components = ( PipelineXMLParser.analyseComponent( node, previous_components, pipeline, prefix), )
                previous_components = final_components
            if node.tag.lower() == PipelineXMLParser.BRANCH_TAG:
                branch_counter += 1
                branch_name = PipelineXMLParser.getAttribute( node, PipelineXMLParser.BRANCH_NAME_ATT, False) 
                if branch_name != None:
                    prefix = initial_prefix + branch_name + ":"
                else:
                    prefix = initial_prefix + "Branch" + str( branch_counter) + ":"
                final_components.extend( PipelineXMLParser.analyseNode( node, previous_components, pipeline, prefix))
            elif node.tag.lower() == PipelineXMLParser.NODE_TAG:
                final_components = PipelineXMLParser.analyseNode( node, previous_components, pipeline, prefix)
                previous_components = []
                for component in final_components:
                    previous_components.append( component)
        
        return final_components


    # --------------------------------------------------------------------------------------
    # Parse a component tag
    @staticmethod
    def analyseComponent( component_node, last_components, pipeline, prefix):
        
        current_component = PipelineXMLParser.getComponent( component_node, prefix)
        current_component.pipelineName = pipeline.name
        current_component.addPreviousComponents( last_components)
        for component in last_components:
            component.addNextComponent( current_component)
        pipeline.addComponent( current_component)
        
        return current_component




    # --------------------------------------------------------------------------------------
    # Verify the component is correctly defined
    # if so, create a Component and return it
    @staticmethod
    def getComponent( node_component, prefix):
        
        processor_name = PipelineXMLParser.getAttribute( node_component, PipelineXMLParser.COMPONENT_PROCESSOR_ATT)
        if processor_name != None and len( processor_name) > 0:
            PipelineXMLParser.RANK += 1
            component = Component( processor_name, str(PipelineXMLParser.RANK), prefix)
            if component != None:
                for node in node_component:
                    if node.tag.lower() == PipelineXMLParser.PARAM_TAG:
                        PipelineXMLParser.getParam( node, component)
                return component
            else:
                raise ParsingException( "PipelineXMLParser.getComponent : Unable to create Component '" + processor_name)            
        else:
            raise ParsingException( "PipelineXMLParser.getComponent : Malformed component - unable to retrieve processor name")            


    # --------------------------------------------------------------------------------------
    # Verify the param is correctly defined
    # if so, add a parameter to the component
    @staticmethod
    def getParam( node_param, component):
        
        param_name = PipelineXMLParser.getAttribute( node_param, PipelineXMLParser.PARAM_NAME_ATT,)
        param_value = PipelineXMLParser.getAttribute( node_param, PipelineXMLParser.PARAM_VALUE_ATT)
        if param_name != None and len( param_name) > 0:
           if param_value != None and len( param_value) > 0:
               component.addParameters( param_name, param_value)
           else:
               raise ParsingException( "PipelineXMLParser.getParam : Malformed parameter - unable to retrieve parameter value in component '" +  component.processorName + "'")
        else:
            raise ParsingException( "PipelineXMLParser.getParam : Malformed parameter - unable to retrieve parameter name in component '" +  component.processorName + "'")


    # --------------------------------------------------------------------------------------
    # Returns the attribute if exists and raise an exception if the attribute is required but not found
    @staticmethod
    def getAttribute( node, att_name, required = True):
        
        try:
            att_value =  node.get( att_name)
            return att_value
        except Exception, exce:
            if required:
                raise ParsingException( "PipelineXMLParser.getAttribute : Node '" + node.tag + "' does not know the attribute :'" + att_name + "'. From:\n\t---> " + str( exce))
            else:
                return None


    # --------------------------------------------------------------------------------------
    # Export the pipelines to an XML file
    @staticmethod
    def toXMLFile( outpath, pipelines):
        
        pipelines_element = Element( PipelineXMLParser.PIPELINES_TAG)
        
        for pipeline in pipelines:
            pipeline_element = Element( PipelineXMLParser.PIPELINE_TAG)
            pipelines_element.append( pipeline_element)
            pipeline_element.attrib[ PipelineXMLParser.PIPELINE_NAME_ATT] = pipeline.name
            
            for component in pipeline.componentList:
                component_element = Element( PipelineXMLParser.COMPONENT_TAG)
                pipeline_element.append( component_element)
                component_element.attrib[ PipelineXMLParser.COMPONENT_PROCESSOR_ATT] = component.processorName
                for param_name, param_value in component.parameters.iteritems():
                    param_element = Element( PipelineXMLParser.PARAM_TAG)
                    component_element.append( param_element)
                    param_element.attrib[ PipelineXMLParser.PARAM_NAME_ATT] = str( param_name)
                    param_element.attrib[ PipelineXMLParser.PARAM_VALUE_ATT] = str( param_value)

        try:
            PipelineXMLParser.indent( pipelines_element, 0)
	    outfile = os.path.join( outpath, pipeline.name + ".xml")
            ElementTree( pipelines_element).write( outfile)
        except IOError, exce:
            Log.log( "PipelineXMLParser.toXMLFile : Unable to write Pipelines to XML file. From:\n\t---> " + str( exce))
        except ParsingException, par_exce:
            Log.log( "PipelineXMLParser.toXMLFile : Unable to save Pipelines to XML file. From:\n\t---> " + str( par_exce))


    # --------------------------------------------------------------------------------------
    # Add indentation to the ElementTree in order to have a pretty print
    # in the XML file (used by subclasses)
    @staticmethod
    def indent( elem, level=0):
            
        i = "\n" + level*"  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for elem in elem:
                PipelineXMLParser.indent(elem, level+1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i   

    # --------------------------------------------------------------------------------------
    # XML Tag used in Pipeline declaration
    
    PIPELINES_TAG = "pipelines"
    PIPELINE_TAG = "pipeline"
    PIPELINE_NAME_ATT = "name"
    NODE_TAG = "node"
    BRANCH_TAG = "branch"
    BRANCH_NAME_ATT = "name"
    COMPONENT_TAG = "component"
    COMPONENT_PROCESSOR_ATT = "processor"
    COMPONENT_INPUT_FILE_ATT = "inputFile"
    PARAM_TAG = "param"
    PARAM_NAME_ATT = "name"
    PARAM_VALUE_ATT = "value"
    
    
# eflag: FileType = Python2
