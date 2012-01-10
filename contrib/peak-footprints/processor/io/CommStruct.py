
from utils.exception.ParsingException import ParsingException

class CommStruct:
    
    PROCESSOR_NAME_ATT = "processorName"
    PROCESSOR_PIPELINE_ATT = "pipeline"
    PROCESSOR_RANK_ATT = "processorRank"
    PROCESSOR_BRANCH_ATT = "processorBranch"

    # --------------------------------------------------------------------------------------
    def __init__( self):
        
        self.origin = CommStructOrigin()


    # --------------------------------------------------------------------------------------
    def setOrigin(self, component):

        self.origin.processorPipeline = component.pipelineName
        self.origin.processorName = component.processorName
        self.origin.processorRank = component.rank
        self.origin.processorBranch = component.branch[0:-1]

    # ###############################
    # FUNCTIONS TO READ THE XLM FILE
    # ###############################


    # --------------------------------------------------------------------------------------
    # Get the attribute stored at top level tag in CommStruct
    def getCommStructData(self, node):
        
        self.origin.processorName = self.getAttribute( node, CommStruct.PROCESSOR_NAME_ATT)
        self.origin.processorPipeline = self.getAttribute( node, CommStruct.PROCESSOR_PIPELINE_ATT)
        self.origin.processorRank = self.getAttributeAsint( node, CommStruct.PROCESSOR_RANK_ATT)
        self.origin.processorBranch = self.getAttribute( node, CommStruct.PROCESSOR_BRANCH_ATT)



    # --------------------------------------------------------------------------------------
    @staticmethod
    def fromXMLFile( self,  input_filepath):

        raise ParsingException( "The method fromXMLFile must be implemented at the inherited class level")
        return None



    # ###############################
    # FUNCTIONS TO WRITE THE XLM FILE
    # ###############################


    # --------------------------------------------------------------------------------------
    # Get the attribute stored at top level tag in CommStruct
    def addCommStructDataToElement(self, element):
        
        element.attrib[ CommStruct.PROCESSOR_NAME_ATT] = self.origin.processorName
        element.attrib[ CommStruct.PROCESSOR_PIPELINE_ATT] = self.origin.processorPipeline
        element.attrib[ CommStruct.PROCESSOR_RANK_ATT] = str( self.origin.processorRank)
        element.attrib[ CommStruct.PROCESSOR_BRANCH_ATT] = self.origin.processorBranch


    # --------------------------------------------------------------------------------------
    def toXMLFile( self, output_filepath):
        
        raise ParsingException( "The method toXMLFile must be implemented at the inherited class level")



    # #############################
    # COMMON FUNCTIONS
    # #############################


    # --------------------------------------------------------------------------------------
    # Returns the attribute if exists and raise an exception if the attribute is required but not found
    @staticmethod
    def getAttribute( node, att_name, required = True):
        
        try:
            att_value =  node.get( att_name)
            return att_value
        except Exception, exce:
            if required:
                raise ParsingException( "CommStruct.getAttribute : Node '" + node.tag + "' does not know the attribute :'" + att_name + "'. From:\n\t---> " + str( exce))
            else:
                return None
    
    
    
    # --------------------------------------------------------------------------------------
    # Return the value of the attribute as int if exists and raise an exception if the attribute is required but not found
    @staticmethod
    def getAttributeAsint( node, att_name, required = True):
        
        try:
            att_value = int( float( CommStruct.getAttribute( node, att_name, required)))
            return att_value
        except (TypeError, ValueError), val_exce:
            if required:
                raise ParsingException( "CommStruct.getAttributeAsint : Unable to convert the value of attribute :'" + att_name + "'. From:\n\t---> " + str( val_exce))
            else:
                return None



    # --------------------------------------------------------------------------------------
    # Return the value of the attribute as float if exists and raise an exception if the attribute is required but not found
    @staticmethod
    def getAttributeAsfloat( node, att_name, required = True):
        
        try:
            att_value = float( CommStruct.getAttribute( node, att_name, required))
            return att_value
        except (TypeError, ValueError), val_exce:
            if required:
                raise ParsingException( "CommStruct.getAttributeAsint : Unable to convert the value of parameter :'" + att_name + "'. From:\n\t---> " + str( val_exce))
            else:
                return None
    
    
    
    # --------------------------------------------------------------------------------------
    # Add indentation to the ElementTree in order to have a pretty print
    # in the XML file (used by subclasses)
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
    
    
    # --------------------------------------------------------------------------------------
    # Return a string representation of the CommStruct
    def toString( self):
        
        return "CommStruct toString() method not implemented"   



# ##############################################
#  CLASS COMMSTRUCTORIGIN
# ##############################################


class CommStructOrigin:
    
    def __init__( self):
        
        self.processorPipeline = ""
        self.processorName = ""
        self.processorBranch = ""
        self.processorRank = 0
    
    
    # --------------------------------------------------------------------------------------
    # Return a string representation of the CommStructOrigin
    def toString( self):
        
        result = "CommStructOrigin " + str( self.processorRank) + "-" + self.processorBranch + ":" + self.processorName + " (" + self.processorPipeline + ")"

        return result
        
        
