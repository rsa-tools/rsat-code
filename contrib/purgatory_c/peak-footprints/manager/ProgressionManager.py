
import os, math, time

from xml.etree.ElementTree import Element
from xml.etree import cElementTree as ET

from utils.Constants import Constants
from utils.log.Log import Log
from utils.FileUtils import FileUtils

class ProgressionManager:
    
    instance = None
    
    # Constants describing possible statuses
    NOT_STARTED_STATUS = "Not Started"
    RUNNING_STATUS = "Running"
    FINISHED_STATUS = "Finished"
    FINISHED_ERROR_STATUS = "Finished with errors"
    EXECUTED_STATUS = "Executed"
    RESUMED_STATUS = "Resumed"
    FAILED_STATUS = "Failed"
    

    # --------------------------------------------------------------------------------------
    # XML tags for output
    
    PIPELINES_TAG = "pipelines"
    PIPELINE_TAG = "pipeline"
    COMPONENT_TAG = "component"
    TASK_TAG = "task"
    STATUS_ATT = "status"
    START_TIME_ATT = "startTime"
    END_TIME_ATT = "endTime"
    NAME_ATT = "name"
    DISPLAY_NAME_ATT = "displayName"
    BRANCH_ATT = "branch"
    RANK_ATT = "rank"
    PROGRESSION_VALUE_ATT = "progression"
    ELAPSED_TIME_ATT = "elapsedTime"
    RESULT_ATT = "result"

    
    
    # --------------------------------------------------------------------------------------
    def __init__( self, output_path, config_path):

        self.pipelinesProgressions = {}
        self.outputPath = output_path
        #self.stylesheetPath = os.path.join( config_path, "resources/xsl/progression/progression.xsl")
        self.stylesheetPath = Constants.PROGRESSION_XSL_FILE



    # --------------------------------------------------------------------------------------
    # Initialize the manager, create the progression object struture relative to the provided
    # pipelines components structure
    @staticmethod
    def initialize( pipelines, output_path, config_path):
        
        ProgressionManager.instance = ProgressionManager( output_path, config_path)
        
        for pipeline in pipelines:
            ProgressionManager.instance.pipelinesProgressions[ pipeline] = PipelineProgression()
            for component in pipeline.componentList:
                ProgressionManager.instance.addNewComponentProgression( component, pipeline)



    # --------------------------------------------------------------------------------------
    # Add a progression object relative to the given pipeline component
    def addNewComponentProgression( self, component, pipeline):
        
        self.pipelinesProgressions[ pipeline].componentProgressions.append( ComponentProgression( component))



    # --------------------------------------------------------------------------------------
    # Set the progression status of the given pipeline to the given value
    @staticmethod
    def setPipelineStatus( pipeline, status):

        if pipeline in ProgressionManager.instance.pipelinesProgressions.keys():
            ProgressionManager.instance.pipelinesProgressions[ pipeline].setStatus( status)
            ProgressionManager.outputProgression( pipeline)



    # --------------------------------------------------------------------------------------
    # Set the component status of the given component to the given value
    @staticmethod
    def setComponentStatus( component, status):
        
        for pipeline in ProgressionManager.instance.pipelinesProgressions.keys():
            pipeline_prog = ProgressionManager.instance.pipelinesProgressions[ pipeline]
            for component_prog in pipeline_prog.componentProgressions:
                if component_prog.component == component:
                    component_prog.setStatus( status)
                    ProgressionManager.outputProgression( pipeline)
                    return
        


    # --------------------------------------------------------------------------------------
    # Set the progression value of the given component to the given value
    @staticmethod
    def setComponentProgression( component, value):
        
        for pipeline in  ProgressionManager.instance.pipelinesProgressions.keys():
            pipeline_prog = ProgressionManager.instance.pipelinesProgressions[ pipeline]
            for component_prog in pipeline_prog.componentProgressions:
                if component_prog.component == component:
                    component_prog.setProgression( value)
                    ProgressionManager.outputProgression( pipeline)
                    return
        
        
    # --------------------------------------------------------------------------------------
    # Set the progression value of the given compone task to the given value
    @staticmethod
    def setTaskProgression( task, component, value):
        
        for pipeline in  ProgressionManager.instance.pipelinesProgressions.keys():
            pipeline_prog = ProgressionManager.instance.pipelinesProgressions[ pipeline]
            for component_prog in pipeline_prog.componentProgressions:
                if component_prog.component == component:
                    component_prog.setTaskProgression( task, value)
                    ProgressionManager.outputProgression( pipeline)
                    return



    # --------------------------------------------------------------------------------------
    # Output the progression status to XML file
    @staticmethod
    def outputProgression( pipeline):

        try:
            # create the pipeline element and set its attributes
            pipeline_element = Element( ProgressionManager.PIPELINE_TAG)
            pipeline_element.attrib[ ProgressionManager.NAME_ATT] = pipeline.name
            pipeline_prog = ProgressionManager.instance.pipelinesProgressions[pipeline]
            pipeline_element.attrib[ ProgressionManager.STATUS_ATT] = pipeline_prog.status
            pipeline_element.attrib[ ProgressionManager.START_TIME_ATT] = time.strftime("%b %d %Y %H:%M:%S", time.localtime( pipeline_prog.startTime))
            if pipeline_prog.status == ProgressionManager.RUNNING_STATUS or pipeline_prog.status == ProgressionManager.NOT_STARTED_STATUS:
                pipeline_element.attrib[ ProgressionManager.END_TIME_ATT] = "0"
            else:
                pipeline_element.attrib[ ProgressionManager.END_TIME_ATT] = time.strftime("%b %d %Y %H:%M:%S", time.localtime( pipeline_prog.endTime))
            pipeline_elapsed_time = pipeline_prog.getElapsedTime()
            if pipeline_elapsed_time > 0:
                pipeline_element.attrib[ ProgressionManager.ELAPSED_TIME_ATT] = str( pipeline_elapsed_time)
            
            # Parse the component list to create the component element and set their attributes
            for component_prog in pipeline_prog.componentProgressions:
                component_element = Element( ProgressionManager.COMPONENT_TAG)
                pipeline_element.append( component_element)
                component_element.attrib[ ProgressionManager.NAME_ATT] = component_prog.component.processorShortName
                component_element.attrib[ ProgressionManager.DISPLAY_NAME_ATT] = component_prog.component.processorDisplayName
                component_element.attrib[ ProgressionManager.BRANCH_ATT] = component_prog.component.branch
                component_element.attrib[ ProgressionManager.RANK_ATT] = component_prog.component.rank 
                component_element.attrib[ ProgressionManager.STATUS_ATT] = component_prog.status
                component_elapsed_time = component_prog.getElapsedTime()
                if component_elapsed_time >= 0:
                    component_element.attrib[ ProgressionManager.ELAPSED_TIME_ATT] = ProgressionManager.convertTime( component_elapsed_time)
                if component_prog.status == ProgressionManager.RUNNING_STATUS:
                    # If component is running, look for tasks to create corresponding elements and attributes
                    if len( component_prog.tasks) > 0:
                        for task in component_prog.tasks:
                            task_element = Element( ProgressionManager.TASK_TAG)
                            component_element.append( task_element)
                            task_element.attrib[ ProgressionManager.NAME_ATT] = task
                            task_element.attrib[ ProgressionManager.PROGRESSION_VALUE_ATT] = str( int( math.ceil( component_prog.taskProgression[task]*1000.0) )/float(10)) + "%"
                    else:
                        # If no task exists, set the progrssion attribute at the component level
                        component_element.attrib[ ProgressionManager.PROGRESSION_VALUE_ATT] = str( int( math.ceil( component_prog.getProgression()*1000.0) )/float(10)) + "%"
                # If component is not running and is not 'not started', set the output put result
                elif component_prog.status != ProgressionManager.NOT_STARTED_STATUS:
                    component_element.attrib[ ProgressionManager.RESULT_ATT] = component_prog.component.getOutputFilePath()

            ProgressionManager.indent( pipeline_element, 0)
            doc = ET.ElementTree( pipeline_element)
            pipeline_output_dir = os.path.join( ProgressionManager.instance.outputPath, pipeline.name)
            progression_file = os.path.join( pipeline_output_dir, Constants.PROGRESSION_XML_FILE)
            #ElementTree( pipeline_element).write( progression_file)
            outfile = FileUtils.openFile( progression_file, 'w')
            outfile.write('<?xml version="1.0" encoding="utf-8"?>\n')
            outfile.write('<?xml-stylesheet type="text/xsl" href="' + ProgressionManager.instance.stylesheetPath + '"?>\n')
            doc.write( outfile)
            outfile.close()

        except IOError, exce:
            Log.log( "ProgressionManager.outputProgression : Unable to write progresssion to XML file. From:\n\t---> " + str( exce))


    # return a h:m:s representation of given duration
    @staticmethod
    def convertTime( duration):
        
        hours = int( duration / float(3600))
        minutes = int( (duration - hours*3600) / float(60))
        seconds = int( duration - hours*3600 - minutes*60)
        
        return "%.2d" %hours + ":" + "%.2d" %minutes + ":" + "%.2d" %seconds
        

    
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
                ProgressionManager.indent(elem, level+1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i   
    

            
# ################################################

class PipelineProgression:
    

    # --------------------------------------------------------------------------------------
    def __init__( self):
        
        self.status = ProgressionManager.NOT_STARTED_STATUS
        self.componentProgressions = []
        self.startTime = 0.0
        self.endTime = 0.0
    
    
    # --------------------------------------------------------------------------------------
    # Set the pipeline progression status to the given value. If this value mean a failure,
    # the status of running component is also set to failure
    def setStatus(self, status):
        
        self.status = status
        
        if self.status == ProgressionManager.RUNNING_STATUS:
            self.startTime = time.time()
        elif self.status == ProgressionManager.FAILED_STATUS:
            self.endTime = time.time()
            for component_prog in self.componentProgressions:
                if component_prog.status == ProgressionManager.RUNNING_STATUS:
                    component_prog.setStatus( ProgressionManager.FAILED_STATUS)
        else:
            self.endTime = time.time()
    
    
    # --------------------------------------------------------------------------------------
    # Add a progression linked to a component to the list
    def addComponentProgression(self, component_progression):
        
        if component_progression != None:
            self.componentProgressions.append( component_progression)
    
    
    # --------------------------------------------------------------------------------------
    # Return the percentage of executed components
    def getProgression( self):
        
        executed = 0
        for component in self.componentProgressions:
            if component.status == ProgressionManager.FINISHED_STATUS:
                executed += 1
                
        return (executed / float( len( self.components))) * 100.0


    # --------------------------------------------------------------------------------------
    # Return the number of seconds elapsed since component start
    def getElapsedTime(self):
        
        return int((self.endTime - self.startTime) *100.0) / 100.0

# ################################################

class ComponentProgression:

    # --------------------------------------------------------------------------------------
    def __init__( self, component):
        
        self.component = component
        self.status = ProgressionManager.NOT_STARTED_STATUS
        self.currentProgress = 0
        self.tasks = []
        self.taskProgression = {}
        self.startTime = 0.0
        self.currentTime = 0.0


    # --------------------------------------------------------------------------------------
    # Set the component progression status to the given value. If this value mean a failure,
    # the status of running tasks is also set to failure
    def setStatus(self, status):
        
        self.status = status
        
        if self.status == ProgressionManager.RUNNING_STATUS:
            self.startTime = time.time()
        elif self.status == ProgressionManager.RESUMED_STATUS:
            self.startTime = time.time()
            self.endTime = self.startTime
        else:
            self.currentTime = time.time()


    # --------------------------------------------------------------------------------------
    # Set the progression value to the given task (creates the task if useful)
    def setTaskProgression( self, task, value):
        
        if not task in self.tasks:
            self.tasks.append( task)
            self.taskProgression[ task] = 0.0
            
        self.taskProgression[ task] = value
        self.currentTime = time.time()


    # --------------------------------------------------------------------------------------
    # Set the progression value
    def setProgression(self, value):
        
        self.currentProgress = value
        self.currentTime = time.time()


    # --------------------------------------------------------------------------------------
    # Return the percentage representing the progression of the task
    def getProgression( self):
        
        if self.status == ProgressionManager.FINISHED_STATUS:
            return 1.0
        else:
            return self.currentProgress


    # --------------------------------------------------------------------------------------
    # Return the number of seconds elapsed since component start
    def getElapsedTime(self):
        
        return int((self.currentTime - self.startTime) *100.0) / 100.0


