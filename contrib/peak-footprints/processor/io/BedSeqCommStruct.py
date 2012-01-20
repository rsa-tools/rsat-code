
from processor.io.CommStruct import CommStruct

from utils.exception.ParsingException import ParsingException
from utils.log.Log import Log

from xml.etree.ElementTree import parse
from xml.etree.ElementTree import Element, ElementTree

from common.BEDSequence import BEDSequence

class BedSeqCommStruct( CommStruct):

    BEDSEQUENCES_NUMBER = "BEDSequencesNumber"
    BEDSEQUENCES_MIN_SIZE = "BEDSequencesMinSize"
    BEDSEQUENCES_MAX_SIZE = "BEDSequencesMaxSize"
    BEDSEQUENCES_MEAN_SIZE = "BEDSequencesMeanSize"
    BEDSEQUENCES_TOTAL_SIZE = "BEDSequencesTotalSize"
    BED_SEQUENCES_SIZE_PATH = "BEDSequencesSizeHistogram"
    BED_SEQUENCES_SIZE_GRAPH_PATH = "BEDSequencesSizeHistogramGraph"

    # --------------------------------------------------------------------------------------
    def __init__( self):
        
        CommStruct.__init__( self)
        self.baseSpecies = ""
        # dicitonnary of the BED Sequences : key = <species>.<chromosom> / value = list of BEDSequence
        self.bedSequencesDict = {}
        # dictionnary of statistics parameters: key = parameter name / Value = parameter value
        self.paramStatistics = {}        


    # --------------------------------------------------------------------------------------
    # Add a sequence in BED format to the CommStruct
    def addBEDSequence(self, bed_sequence):
        
        if bed_sequence != None:
            key = bed_sequence.getKey()
            if not self.bedSequencesDict.has_key( key):
                self.bedSequencesDict[ key] = []
            self.bedSequencesDict[ key].append( bed_sequence)


    # --------------------------------------------------------------------------------------
    # Return a string representation of the CommStruct
    def toString( self):
        
        result = "BedSeqCommStruct\n"

        for chrom in self.bedSequencesDict.keys():
            for sequence in self.bedSequencesDict[chrom]:
                result += "|--BED Sequence " + sequence.toString() + "\n"
            
        return result    


    # --------------------------------------------------------------------------------------
    # Read the CommStruct from the given XML file
    @staticmethod
    def fromXMLFile( input_filepath):
        
        try:
            return BedSeqCommStruct.getCommStructFromXML( input_filepath)
        except ParsingException, par_exce:
            Log.log( "BedSeqCommStruct.fromXMLFile : unable to get CommStruct from XML file '" + input_filepath + "'. From:\n\t--->" + str( par_exce))
            return None


    # --------------------------------------------------------------------------------------
    # Write the CommStruct to the given XML file
    def toXMLFile( self, output_filepath):
        
        try:
            root_element = self.convertCommStructToElementTree()
            self.indent( root_element, 0)
            ElementTree( root_element).write( output_filepath)
        except IOError, exce:
            Log.log( "BedSeqCommStruct.toXMLFile : Unable to write CommStruct to XML file. From:\n\t---> " + str( exce))
        except ParsingException, par_exce:
            Log.log( "BedSeqCommStruct.toXMLFile : Unable to save CommStruct to XML file. From:\n\t---> " + str( par_exce))


    # #############################
    # METHODS TO READ THE XLM FILE
    # #############################

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
            raise ParsingException( "BedSeqCommStruct.getCommStructFromXML : Unable to open/close XML file '" + commstruct_filepath, "'. From:\n\t---> " + str( io_exce))
        
        comm_struct = BedSeqCommStruct()
        
        comm_struct.getCommStructData( root_element)
        
        for root_son in root_element:
            if root_son.tag.lower() == BedSeqCommStruct.BEDSEQ_TAG:
                bed_sequence = BedSeqCommStruct.getBEDSequence( root_son, comm_struct)
                base_species = bed_sequence.species
            elif root_son.tag.lower() == BedSeqCommStruct.STATISTICS_TAG:
                BedSeqCommStruct.getStatistics( root_son, comm_struct)
            else:
                raise ParsingException( "BedSeqCommStruct.getCommStructFromXML : The data contains an unauthorized element : '" + root_son.tag.lower() +  "'")                    
        
        comm_struct.baseSpecies = base_species
        
        return comm_struct
   
    # --------------------------------------------------------------------------------------
    # Retrieve the motif statistics
    @staticmethod
    def getStatistics( statistics_node, comm_struct):
        
        for son_node in statistics_node:
            if son_node.tag.lower() == BedSeqCommStruct.PARAM_TAG:
                att_name = CommStruct.getAttribute( son_node, BedSeqCommStruct.PARAM_NAME_ATT)
                att_value = CommStruct.getAttribute( son_node, BedSeqCommStruct.PARAM_VALUE_ATT)
                comm_struct.paramStatistics[ att_name] = att_value
            else:
                raise ParsingException( "BedSeqCommStruct.getStatistics : The statistics contains an unauthorized element : '" + son_node.tag.lower() +  "'")                    



    # --------------------------------------------------------------------------------------
    # Retrieve the BED sequence information, create a BEDSequence
    # and store it in the CommStruct
    @staticmethod
    def getBEDSequence( node_bedseq, comm_struct):
        
        species = CommStruct.getAttribute( node_bedseq, BedSeqCommStruct.BEDSEQ_SPECIES_ATT)
        chrom = CommStruct.getAttribute( node_bedseq, BedSeqCommStruct.BEDSEQ_CHROM_ATT)
        start = CommStruct.getAttributeAsint( node_bedseq, BedSeqCommStruct.BEDSEQ_START_ATT)
        end = CommStruct.getAttributeAsint( node_bedseq, BedSeqCommStruct.BEDSEQ_END_ATT)
        score = CommStruct.getAttributeAsint( node_bedseq, BedSeqCommStruct.BEDSEQ_SCORE_ATT, False)
        max = CommStruct.getAttributeAsint( node_bedseq, BedSeqCommStruct.BEDSEQ_MAX_PEAK_ATT, False)
        id = CommStruct.getAttribute( node_bedseq, BedSeqCommStruct.BEDSEQ_ID_ATT, False)
        
        if species != None and chrom != None and start != None and end != None:
            bed_sequence = BEDSequence( species, chrom, start, end)
            if score != None:
                bed_sequence.score = score
            if max != None:
                bed_sequence.referenceIndex = max
            if id != None:
                bed_sequence.id = id
                
            comm_struct.addBEDSequence( bed_sequence)
            return bed_sequence
        else:
            raise ParsingException ( "BedSeqCommStruct.getBEDSequence : Malformed BED Sequence - unable to retrieve sequence information")


    # #############################
    # METHODS TO WRITE THE XLM FILE
    # #############################

    # Convert the CommStruct to an ElementTree containing all the required
    # suitably organized information
    def convertCommStructToElementTree( self):
        
        data_element = Element( BedSeqCommStruct.DATA_TAG)
        
        self.addCommStructDataToElement( data_element)

        # Output the statistics data
        statistics_element = Element( BedSeqCommStruct.STATISTICS_TAG)
        data_element.append( statistics_element)
        for param_name in sorted( self.paramStatistics.keys()):
            paramstats_element = Element( BedSeqCommStruct.PARAM_TAG)
            statistics_element.append( paramstats_element)
            paramstats_element.attrib[ BedSeqCommStruct.PARAM_NAME_ATT] = param_name
            paramstats_element.attrib[ BedSeqCommStruct.PARAM_VALUE_ATT] = str( self.paramStatistics[ param_name])

        # output the bedseq data
        for chrom in self.bedSequencesDict.keys():
            for bedseq in self.bedSequencesDict[chrom]:
                bedseq_element = Element( BedSeqCommStruct.BEDSEQ_TAG)
                data_element.append( bedseq_element)
                bedseq_element.attrib[ BedSeqCommStruct.BEDSEQ_SPECIES_ATT] = bedseq.species
                bedseq_element.attrib[ BedSeqCommStruct.BEDSEQ_CHROM_ATT] = bedseq.chromosom
                bedseq_element.attrib[ BedSeqCommStruct.BEDSEQ_START_ATT] = str( bedseq.indexStart)
                bedseq_element.attrib[ BedSeqCommStruct.BEDSEQ_END_ATT] = str( bedseq.indexEnd)
                bedseq_element.attrib[ BedSeqCommStruct.BEDSEQ_SCORE_ATT] = str( bedseq.score)
                if bedseq.referenceIndex != 0:
                    bedseq_element.attrib[ BedSeqCommStruct.BEDSEQ_MAX_PEAK_ATT] = str( bedseq.referenceIndex)
                if bedseq.id != "":
                    bedseq_element.attrib[ BedSeqCommStruct.BEDSEQ_ID_ATT] = str( bedseq.id)

        return data_element
        
        

    # Definition of the tags and attributes names used in the XML format

    DATA_TAG = "data"
    DATA_SPECIES_ATT = "species"
    DATA_GENE_ATT = "gene"
    DATA_TF_ATT = "tf"
    
    STATISTICS_TAG = "statistics"
    PARAM_TAG = "param"
    PARAM_NAME_ATT = "name"
    PARAM_VALUE_ATT = "value"    
    
    BEDSEQ_TAG = "bedseq"
    BEDSEQ_SPECIES_ATT = "species"
    BEDSEQ_CHROM_ATT = "chr"
    BEDSEQ_START_ATT = "start"
    BEDSEQ_END_ATT = "end"
    BEDSEQ_SCORE_ATT = "score"
    BEDSEQ_ID_ATT = "id"
    BEDSEQ_MAX_PEAK_ATT = "max"


# eflag: FileType = Python2
