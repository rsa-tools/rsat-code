
from processor.io.CommStruct import CommStruct

from xml.etree.ElementTree import parse
from xml.etree.ElementTree import Element, ElementTree

from utils.log.Log import Log
from utils.exception.ParsingException import ParsingException

from common.Motif import Motif

class MotifStatisticsCommStruct( CommStruct):
    
    # --------------------------------------------------------------------------------------
    def __init__( self):
        
        CommStruct.__init__( self)
        self.motifList = []
        self.motifToStatistics = {}
        


    # --------------------------------------------------------------------------------------
    def addMotif( self, motif):
        
        if motif != None:
            self.motifList.append( motif)
    
    # --------------------------------------------------------------------------------------
    def addMotifStatistics(self, motif, statistics):
        
        if motif != None and statistics != None:
            self.motifToStatistics[motif] = statistics
            
    
    # --------------------------------------------------------------------------------------
    def toXMLFile( self, output_filepath):
        
        try:
            root_element = self.convertCommStructToElementTree()
            self.indent( root_element,  0)
            ElementTree( root_element).write( output_filepath)
        except IOError, exce:
            Log.log( "MotifStatisticsCommStruct.toXMLFile : Unable to write CommStruct to XML file. From:\n\t---> " + str( exce))
        except ParsingException, par_exce:
            Log.log( "MotifStatisticsCommStruct.toXMLFile : Unable to save CommStruct to XML file. From:\n\t---> " + str( par_exce))


    # --------------------------------------------------------------------------------------
    @staticmethod
    def fromXMLFile( input_filepath):

        try:
            return MotifStatisticsCommStruct.getCommStructFromXML( input_filepath)
        except ParsingException, par_exce:
            Log.log( "MotifStatisticsCommStruct.fromXMLFile : Unable to get CommStruct from XML file '" + input_filepath + "'. From:\n\t---> " + str( par_exce))
            return None



    # #############################
    # METHODS TO READ THE XLM FILE
    # #############################

    # --------------------------------------------------------------------------------------
    @staticmethod
    def getCommStructFromXML( commstruct_filepath):
        
        file = None
        root_element = None
        
        try:
            file = open( commstruct_filepath, "r")
            tree = parse( file)
            root_element = tree.getroot()
            file.close()
        except IOError, io_exce:
            raise ParsingException( "MotifStatisticsCommStruct.getCommStructFromXML : Unable to open/close XML file '" + commstruct_filepath, "' : " + str( io_exce))
        
        comm_struct = MotifStatisticsCommStruct()
        
        comm_struct.getCommStructData( root_element)
        
        for node_motif in root_element:
            if node_motif.tag.lower() == MotifStatisticsCommStruct.MOTIF_TAG:
                MotifStatisticsCommStruct.getMotif( node_motif, comm_struct)
        
        return comm_struct
        
        
    # --------------------------------------------------------------------------------------
    # Retrieve the list of motifs of the alignment
    # and store it into the given SequenceAlignment
    @staticmethod
    def getMotif( node_motif, comm_struct):
        
            name = CommStruct.getAttribute( node_motif, MotifStatisticsCommStruct.MOTIF_NAME_ATT)
            consensus = CommStruct.getAttribute( node_motif, MotifStatisticsCommStruct.MOTIF_CONSENSUS_ATT, False)
            
            if name != None:
                motif = Motif( 0, 0, name, None)
                if consensus != None:
                    motif.consensus = consensus
                comm_struct.addMotif( motif)
                statistics = MotifStatisticsCommStruct.getMotifStatistics( node_motif, motif)
                comm_struct.addMotifStatistics( motif, statistics)
            else:
                raise ParsingException( "MotifStatisticsCommStruct.getAlignmentMotifs : The motif is missing required attribute 'name'")

    # --------------------------------------------------------------------------------------
    # Verify the param is correctly defined
    # if so, add a parameter to the component
    @staticmethod
    def getMotifStatistics( node_motif, motif):
        
        statistics = MotifStatistics()
        
        for node_param in node_motif:
            if node_param.tag.lower() == MotifStatisticsCommStruct.PARAM_TAG:
                param_name = MotifStatisticsCommStruct.getAttribute( node_param, MotifStatisticsCommStruct.PARAM_NAME_ATT, False)
                param_value = MotifStatisticsCommStruct.getAttribute( node_param, MotifStatisticsCommStruct.PARAM_VALUE_ATT, False)
                if param_name != None and len( param_name) > 0:
                    if param_value != None and len( param_value) > 0:
                        if param_name == MotifStatisticsCommStruct.CHI2_PARAM_NAME:
                            statistics.chi2 = MotifStatisticsCommStruct.getTokenAsfloat( param_value, False)
                        elif param_name == MotifStatisticsCommStruct.HISTOGRAM_GRAPH_PATH_PARAM_NAME:
                            statistics.histogramGraphPath = param_value
                        elif param_name == MotifStatisticsCommStruct.HISTOGRAM_PARAM_NAME:
                            statistics.histogram = param_value.split( MotifStatisticsCommStruct.HISTOGRAM_ENTRY_SEPARATOR_CHAR)
                        elif param_name == MotifStatisticsCommStruct.NULL_HISTOGRAM_PARAM_NAME:
                            statistics.nullHistogram = param_value.split( MotifStatisticsCommStruct.HISTOGRAM_ENTRY_SEPARATOR_CHAR)
                        else:
                            Log.log( "MotifStatisticsCommStruct.getMotifAttributes : Unknown attribute name : " + param_name)
                    else:
                        raise ParsingException( "MotifStatisticsCommStruct.getMotifAttributes : Malformed parameter - unable to retrieve parameter value in motif '" +  motif.name + "'")
                else:
                    raise ParsingException( "MotifStatisticsCommStruct.getMotifAttributes : Malformed parameter - unable to retrieve parameter name in motif '" +  motif.name + "'")

        return statistics

    # --------------------------------------------------------------------------------------
    # Returns the attribute if exists and raise an exception if the attribute is required but not found
    @staticmethod
    def getAttribute( node, att_name, required = True):
        
        try:
            att_value =  node.get( att_name)
            return att_value
        except Exception, exce:
            if required:
                raise ParsingException( "MotifStatisticsCommStruct.getAttribute : Node '" + node.tag + "' does not know the attribute :'" + att_name + "'. From:\n\t---> " + str( exce))
            else:
                return None

    # --------------------------------------------------------------------------------------
    # Return the value of the token as int if possible and raise an exception if the value is required but not converted
    @staticmethod
    def getTokenAsfloat( token, required=True):
        
        try:
            att_value = float( token)
            return att_value
        except (TypeError, ValueError), val_exce:
            if required:
                raise ParsingException( "MotifStatisticsCommStruct.getTokenAsfloat : Unable to convert the token to float :'" + token + "'. From:\n\t---> " + str( val_exce))
            else:
                return None


    # #############################
    # METHODS TO WRITE THE XLM FILE
    # #############################


    # --------------------------------------------------------------------------------------
    # Convert the CommStruct to an ElementTree containing all the required
    # suitably organized information
    def convertCommStructToElementTree( self):
        
        data_element = Element( MotifStatisticsCommStruct.DATA_TAG)
        
        self.addCommStructDataToElement( data_element)
        
        for motif in self.motifList:
            motif_element = Element( MotifStatisticsCommStruct.MOTIF_TAG)
            data_element.append( motif_element)
            motif_element.attrib[ MotifStatisticsCommStruct.MOTIF_NAME_ATT] = motif.name
            if motif.consensus != None and len( motif.consensus) > 0:
                motif_element.attrib[ MotifStatisticsCommStruct.MOTIF_NAME_ATT] = motif.consensus
                
            # Add the statistics
            motif_stats = self.motifToStatistics[ motif]
            if motif_stats.chi2 != None:
                param_chi2_element = Element( MotifStatisticsCommStruct.PARAM_TAG)
                motif_element.append( param_chi2_element)
                param_chi2_element.attrib[ MotifStatisticsCommStruct.PARAM_NAME_ATT] = MotifStatisticsCommStruct.CHI2_PARAM_NAME
                param_chi2_element.attrib[ MotifStatisticsCommStruct.PARAM_VALUE_ATT] = str( motif_stats.chi2)
            if  motif_stats.histogramGraphPath != None:
                param_histopath_element = Element( MotifStatisticsCommStruct.PARAM_TAG)
                motif_element.append( param_histopath_element)
                param_histopath_element.attrib[ MotifStatisticsCommStruct.PARAM_NAME_ATT] = MotifStatisticsCommStruct.HISTOGRAM_GRAPH_PATH_PARAM_NAME
                param_histopath_element.attrib[ MotifStatisticsCommStruct.PARAM_VALUE_ATT] = str( motif_stats.histogramGraphPath)
            if motif_stats.histogram != None:
                line = ""
                for token in motif_stats.histogram:
                    line += token + MotifStatisticsCommStruct.HISTOGRAM_ENTRY_SEPARATOR_CHAR
                line = line[:-1]
                param_histo_element = Element( MotifStatisticsCommStruct.PARAM_TAG)
                motif_element.append( param_histo_element)
                param_histo_element.attrib[ MotifStatisticsCommStruct.PARAM_NAME_ATT] = MotifStatisticsCommStruct.HISTOGRAM_PARAM_NAME
                param_histo_element.attrib[ MotifStatisticsCommStruct.PARAM_VALUE_ATT] = line
            if motif_stats.nullHistogram != None:
                line = ""
                for token in motif_stats.nullHistogram:
                    line += token + MotifStatisticsCommStruct.HISTOGRAM_ENTRY_SEPARATOR_CHAR
                line = line[:-1]
                param_nullhisto_element = Element( MotifStatisticsCommStruct.PARAM_TAG)
                motif_element.append( param_nullhisto_element)
                param_nullhisto_element.attrib[ MotifStatisticsCommStruct.PARAM_NAME_ATT] = MotifStatisticsCommStruct.NULL_HISTOGRAM_PARAM_NAME
                param_nullhisto_element.attrib[ MotifStatisticsCommStruct.PARAM_VALUE_ATT] = line            
        
        return data_element



    # --------------------------------------------------------------------------------------
    # Definition of the tags and attributes names used in the XML format
    
    DATA_TAG = "data"
    MOTIF_TAG = "motif"
    MOTIF_NAME_ATT = "name"
    MOTIF_CONSENSUS_ATT = "consensus"
    PARAM_TAG = "param"
    PARAM_NAME_ATT = "name"
    PARAM_VALUE_ATT = "value"
    CHI2_PARAM_NAME = "Chi2"
    HISTOGRAM_PARAM_NAME = "Histogram"
    NULL_HISTOGRAM_PARAM_NAME = "NullHistogram"
    HISTOGRAM_GRAPH_PATH_PARAM_NAME = "HistogramGraphPath"
    HISTOGRAM_ENTRY_SEPARATOR_CHAR = ";"
    HISTOGRAM_VALUE_SEPARATOR_CHAR = ":"


class MotifStatistics:
    
    def __init__(self):
        
        self.chi2 = None
        self.histogram = None
        self.histogramGraphPath = None
        self.nullHistogram = None
        
