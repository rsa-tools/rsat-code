    
from processor.io.CommStruct import CommStruct

from utils.exception.ParsingException import ParsingException
from utils.log.Log import Log
from utils.Constants import Constants
from utils.exception.ExecutionException import ExecutionException

from xml.etree.ElementTree import parse
from xml.etree.ElementTree import Element, ElementTree

from common.BEDSequence import BEDSequence
from common.SequenceAlignment import SequenceAlignment
from common.Motif import Motif
from common.PWM import PWM

import math

class BedSeqAlignmentStatsCommStruct( CommStruct):

    BEDSEQUENCES_NUMBER = "BEDSequencesNumber"
    BEDSEQUENCES_MIN_SIZE = "BEDSequencesMinSize"
    BEDSEQUENCES_MAX_SIZE = "BEDSequencesMaxSize"
    BEDSEQUENCES_MEAN_SIZE = "BEDSequencesMeanSize"
    BEDSEQUENCES_TOTAL_SIZE = "BEDSequencesTotalSize"
    BED_SEQUENCES_SIZE_PATH = "BEDSequencesSizeHistogram"
    BED_SEQUENCES_SIZE_GRAPH_PATH = "BEDSequencesSizeHistogramGraph"

    REFERENCE_SPECIES = "ReferenceSpecies"
    ALIGNED_SPECIES = "AlignedSpecies"
    
    MSA_NUMBER = "MSANumber" #Number of BEDSequence with MSA
    MSA_MIN_SIZE = "MSAMinSize"
    MSA_MAX_SIZE = "MSAMaxSize"
    MSA_MEAN_SIZE = "MSAMeanSize"
    MSA_TOTAL_SIZE = "MSATotalSize"
    MSA_SIZE_PATH = "MSASizeHistogram"
    MSA_SIZE_GRAPH_PATH = "MSASizeHistogramGraph"
    
    CONSERVED_BLOCKS_NUMBER = "ConservedBlocksNumber"
    CONSERVED_BLOCKS_MIN_SIZE = "ConservedBlocksMinSize"
    CONSERVED_BLOCKS_MAX_SIZE = "ConservedBlocksMaxSize"
    CONSERVED_BLOCKS_MEAN_SIZE = "ConservedBlocksMeanSize"
    CONSERVED_BLOCKS_TOTAL_SIZE = "ConservedBlocksTotalSize"
    CONSERVED_BLOCKS_SIZE_PATH = "ConservedBlocksSizeHistogram"
    CONSERVED_BLOCKS_SIZE_GRAPH_PATH = "ConservedBlocksSizeHistogramGraph"

    
    REFERENCE_MOTIF = "ReferenceMotif"
    CLASSIFICATION_SCORE = "ClassificationScore"
    BED_OUTPUT_PATH = "BEDOutput"
    BIGBED_OUTPUT_PATH = "BigBEDOutput"

    # --------------------------------------------------------------------------------------
    def __init__( self):
        
        CommStruct.__init__( self)
        self.baseSpecies = ""
        # dicitonnary of the BED Sequences : key = chromosom / value = list of BED Sequence
        self.bedSequencesDict = {}
        # dictionnary of the Multiple Alignments : key = BED Sequence / value = list of multiple alignments
        self.bedToMA = {}
        # dictionnary of statistics parameters: key = parameter name / Value = parameter value
        self.paramStatistics = {}
        # dictionnary of identified motifs statistics : key = motif name / Value = MotifStatistics instance
        self.motifStatistics = {}



    # --------------------------------------------------------------------------------------
    # Add a sequence in BED format to the CommStruct
    def addBEDSequence(self, bed_sequence):
        
        if bed_sequence != None:
            key = bed_sequence.getKey()
            if not self.bedSequencesDict.has_key( key):
                self.bedSequencesDict[ key] = []
            self.bedSequencesDict[ key].append( bed_sequence)


    # --------------------------------------------------------------------------------------
    # Add a Sequence Alignement associated to the given sequence in BED format
    def addSequenceAlignment( self, bed_sequence, seq_align):
        
        if bed_sequence != None and seq_align != None:
            if not self.bedToMA.has_key( bed_sequence):
                self.bedToMA[ bed_sequence] = []
            self.bedToMA[ bed_sequence].append( seq_align)


    # --------------------------------------------------------------------------------------
    # Return a string representation of the CommStruct
    def toString( self):
        
        result = "BedSeqAlignmentStatsCommStruct\n"
        for chrom in self.bedSequencesDict.keys():
            for sequence in self.bedSequencesDict[chrom]:
                result += "|--" + sequence.toString() + "\n"
                if self.bedToMA.has_key( sequence):
                    for sa in self.bedToMA[sequence]:
                        result+= sa.toString()
                else:
                    result += "| No Sequence Alignment"
                result += "\n"
            
        return result    


    # --------------------------------------------------------------------------------------
    # Read the CommStruct from the given XML file
    @staticmethod
    def fromXMLFile( input_filepath):
        
        try:
            return BedSeqAlignmentStatsCommStruct.getCommStructFromXML( input_filepath)
        except ParsingException, par_exce:
            Log.log( "BedSeqAlignmentStatsCommStruct.fromXMLFile : Unable to get CommStruct from XML file '" + input_filepath + "'. From:\n\t---> " + str( par_exce))
            return None


    # --------------------------------------------------------------------------------------
    # Write the CommStruct to the given XML file
    def toXMLFile( self, output_filepath):
        
        try:
            root_element = self.convertCommStructToElementTree()
            self.indent( root_element,  0)
            ElementTree( root_element).write( output_filepath)
        except IOError, exce:
            Log.log( "BedSeqAlignmentStatsCommStruct.toXMLFile : Unable to write CommStruct to XML file. From:\n\t---> " + str( exce))
        except ParsingException, par_exce:
            Log.log( "BedSeqAlignmentStatsCommStruct.toXMLFile : Unable to save CommStruct to XML file. From:\n\t---> " + str( par_exce))
    



    # #############################
    # METHODS TO READ THE XLM FILE
    # #############################

    # --------------------------------------------------------------------------------------
    @staticmethod
    def getCommStructFromXML( commstruct_filepath):
        
        commstruct_file = None
        root_element = None
        
        try:
            commstruct_file = open( commstruct_filepath, "r")
            tree = parse( commstruct_file)
            root_element = tree.getroot()
            commstruct_file.close()
        except IOError, io_exce:
            raise ParsingException( "BedSeqAlignmentStatsCommStruct.getCommStructFromXML : Unable to open/close XML commstruct_file '" + commstruct_filepath, "' : " + str( io_exce))
        
        comm_struct = BedSeqAlignmentStatsCommStruct()
        
        comm_struct.getCommStructData( root_element)
        
        base_species = ""
        for root_son in root_element:
            if root_son.tag.lower() == BedSeqAlignmentStatsCommStruct.BEDSEQ_TAG:
                bedseq = BedSeqAlignmentStatsCommStruct.getBEDSequence( root_son, comm_struct)
                base_species = bedseq.species
                for align_node in root_son:
                    if align_node.tag.lower() == BedSeqAlignmentStatsCommStruct.ALIGNMENT_TAG:
                        seqalign = BedSeqAlignmentStatsCommStruct.getAlignment( align_node, bedseq, comm_struct)
                        for sub_node in align_node:
                            if sub_node.tag.lower() == BedSeqAlignmentStatsCommStruct.SEQUENCES_TAG:   
                                BedSeqAlignmentStatsCommStruct.getAlignmentSequences( sub_node, bedseq, seqalign)
                            elif sub_node.tag.lower() == BedSeqAlignmentStatsCommStruct.MOTIFS_TAG:
                                BedSeqAlignmentStatsCommStruct.getAlignmentMotifs( sub_node, bedseq, seqalign)
                            else:
                                raise ParsingException( "BedSeqAlignmentStatsCommStruct.getCommStructFromXML : The alignment of '" + bedseq.toString() + "' contains an unauthorized element : '" + sub_node.tag.lower() +  "'")
                    else:
                        raise ParsingException( "BedSeqAlignmentStatsCommStruct.getCommStructFromXML : The BED Sequence '" + bedseq.toString() + "' contains an unauthorized element : '" + align_node.tag.lower() +  "'")
            elif root_son.tag.lower() == BedSeqAlignmentStatsCommStruct.STATISTICS_TAG:
                BedSeqAlignmentStatsCommStruct.getStatistics( root_son, comm_struct)
            else:
                raise ParsingException( "BedSeqAlignmentStatsCommStruct.getCommStructFromXML : The data contains an unauthorized element : '" + root_son.tag.lower() +  "'")                    
        
        comm_struct.baseSpecies = base_species
        
        return comm_struct


    # --------------------------------------------------------------------------------------
    # Retrieve the motif statistics
    @staticmethod
    def getStatistics( statistics_node, comm_struct):
        
        for son_node in statistics_node:
            if son_node.tag.lower() == BedSeqAlignmentStatsCommStruct.MOTIF_STATS_TAG:
                name = CommStruct.getAttribute( son_node, BedSeqAlignmentStatsCommStruct.MOTIF_STATS_NAME_ATT)
                if name != None and len( name) > 0:
                    motif_id = CommStruct.getAttribute( son_node, BedSeqAlignmentStatsCommStruct.MOTIF_STATS_ID_ATT)
                    family = CommStruct.getAttribute( son_node, BedSeqAlignmentStatsCommStruct.MOTIF_STATS_FAMILY_ATT)
                    classe = CommStruct.getAttribute( son_node, BedSeqAlignmentStatsCommStruct.MOTIF_STATS_CLASS_ATT)
                    motif_type = CommStruct.getAttribute( son_node, BedSeqAlignmentStatsCommStruct.MOTIF_STATS_TYPE_ATT)
                    size = CommStruct.getAttributeAsint( son_node, BedSeqAlignmentStatsCommStruct.MOTIF_STATS_SIZE_ATT)
                    
                    motif_stats = MotifStatistics( name)
                    motif_stats.motifID = motif_id
                    motif_stats.motifFamily = family
                    motif_stats.motifClass = classe
                    motif_stats.motifType = motif_type
                    motif_stats.motifSize = size
                    
                    for param_node in son_node:
                        if param_node.tag.lower() == BedSeqAlignmentStatsCommStruct.PARAM_TAG:
                            att_name = CommStruct.getAttribute( param_node, BedSeqAlignmentStatsCommStruct.PARAM_NAME_ATT)
                            att_value = CommStruct.getAttribute( param_node, BedSeqAlignmentStatsCommStruct.PARAM_VALUE_ATT)
                            motif_stats.setAttribute( att_name, att_value)
                    
                    comm_struct.motifStatistics[ name] = motif_stats
                    
            elif son_node.tag.lower() == BedSeqAlignmentStatsCommStruct.PARAM_TAG:
                att_name = CommStruct.getAttribute( son_node, BedSeqAlignmentStatsCommStruct.PARAM_NAME_ATT)
                att_value = CommStruct.getAttribute( son_node, BedSeqAlignmentStatsCommStruct.PARAM_VALUE_ATT)
                comm_struct.paramStatistics[ att_name] = att_value
            else:
                raise ParsingException( "BedSeqAlignmentStatsCommStruct.getStatistics : The statistics contains an unauthorized element : '" + son_node.tag.lower() +  "'")                    
                


    # --------------------------------------------------------------------------------------
    # Retrieve the BED sequence information, create a BEDSequence
    # and store it in the CommStruct
    @staticmethod
    def getBEDSequence( node_bedseq, comm_struct):
        
        try:
            species = CommStruct.getAttribute( node_bedseq, BedSeqAlignmentStatsCommStruct.BEDSEQ_SPECIES_ATT)
            chrom = CommStruct.getAttribute( node_bedseq, BedSeqAlignmentStatsCommStruct.BEDSEQ_CHROM_ATT)
            start = CommStruct.getAttributeAsint( node_bedseq, BedSeqAlignmentStatsCommStruct.BEDSEQ_START_ATT)
            end = CommStruct.getAttributeAsint( node_bedseq, BedSeqAlignmentStatsCommStruct.BEDSEQ_END_ATT)
            score = CommStruct.getAttributeAsint( node_bedseq, BedSeqAlignmentStatsCommStruct.BEDSEQ_SCORE_ATT, False)
            peak_max = CommStruct.getAttributeAsint( node_bedseq, BedSeqAlignmentStatsCommStruct.BEDSEQ_MAX_PEAK_ATT, False)
            peak_id = CommStruct.getAttribute( node_bedseq, BedSeqAlignmentStatsCommStruct.BEDSEQ_ID_ATT, False)
        except ParsingException,  par_exce:
            raise ParsingException ( "BedSeqAlignmentStatsCommStruct.getBEDSequence : Malformed BED Sequence - some attributes are not numbers. From:\n\t---> " + str( par_exce))
 
        if species != None and chrom != None and start != None and end != None:
            bed_sequence = BEDSequence( species, chrom, start, end)
            if score != None:
                bed_sequence.score = score
            if peak_max != None:
                bed_sequence.referenceIndex = peak_max
            if peak_id != None:
                bed_sequence.id = peak_id
                
            comm_struct.addBEDSequence( bed_sequence)
            return bed_sequence
        else:
            raise ParsingException ( "BedSeqAlignmentStatsCommStruct.getBEDSequence : Malformed BED Sequence - unable to retrieve sequence required attributes")



    # --------------------------------------------------------------------------------------
    # create the SequenceAlignment object
    @staticmethod
    def getAlignment( align_node, bedseq, comm_struct):
        
        name = CommStruct.getAttribute( align_node, BedSeqAlignmentStatsCommStruct.ALIGNMENT_NAME_ATT)
        seq_align = SequenceAlignment()
        seq_align.name = name
        seq_align.referenceSpecies = bedseq.species
        comm_struct.addSequenceAlignment( bedseq, seq_align)
        
        return seq_align


    # --------------------------------------------------------------------------------------
    # Retrieve the list of sequence of the alignment
    # and store it into the given SequenceAlignment
    @staticmethod
    def getAlignmentSequences( sub_node, bedseq, seqalign):
        
        for node_sequence in sub_node:
            if node_sequence.tag.lower() == BedSeqAlignmentStatsCommStruct.SEQUENCE_TAG:
                species = CommStruct.getAttribute( node_sequence, BedSeqAlignmentStatsCommStruct.SEQUENCE_SPECIES_ATT)
                text = list( CommStruct.getAttribute( node_sequence, BedSeqAlignmentStatsCommStruct.SEQUENCE_TEXT_ATT))
                
                if species != None and text != None:
                    seqalign.addSequence( species, text)
                else:
                    raise ParsingException( "BedSeqAlignmentStatsCommStruct.getAlignmentSequences : A sequence of the alignment of '" + bedseq.toString() + "' is missing required attributes")
                    
            else:
                raise ParsingException( "BedSeqAlignmentStatsCommStruct.getAlignmentSequences : The sequences of the alignment of '" + bedseq.toString() + "' contains an unauthorized element : '" + node_sequence.tag.lower() +  "'")


    # --------------------------------------------------------------------------------------
    # Retrieve the list of motifs of the alignment
    # and store it into the given SequenceAlignment
    @staticmethod
    def getAlignmentMotifs( sub_node, bedseq, seqalign):
        
        for node_motif in sub_node:
            if node_motif.tag.lower() == BedSeqAlignmentStatsCommStruct.MOTIF_TAG:
                
                start = CommStruct.getAttributeAsint( node_motif, BedSeqAlignmentStatsCommStruct.MOTIF_START_ATT)
                end = CommStruct.getAttributeAsint( node_motif, BedSeqAlignmentStatsCommStruct.MOTIF_END_ATT)
                name = CommStruct.getAttribute( node_motif, BedSeqAlignmentStatsCommStruct.MOTIF_NAME_ATT)
                motif_id = CommStruct.getAttribute( node_motif, BedSeqAlignmentStatsCommStruct.MOTIF_ID_ATT, False)
                consensus = CommStruct.getAttribute( node_motif, BedSeqAlignmentStatsCommStruct.MOTIF_CONSENSUS_ATT, False)
                nb_species = CommStruct.getAttributeAsint( node_motif, BedSeqAlignmentStatsCommStruct.MOTIF_NBSPECIES_ATT, False)
                strand = CommStruct.getAttribute( node_motif, BedSeqAlignmentStatsCommStruct.MOTIF_STRAND_ATT, False)
                offset = CommStruct.getAttributeAsint( node_motif, BedSeqAlignmentStatsCommStruct.MOTIF_OFFSET_ATT)
                score = CommStruct.getAttributeAsfloat( node_motif, BedSeqAlignmentStatsCommStruct.MOTIF_SCORE_ATT)
                
                # Retrieve the PWM of the motif

                pwm_s1 = CommStruct.getAttribute( node_motif, BedSeqAlignmentStatsCommStruct.MOTIF_PWM_ATT)
                if pwm_s1 != None and len( pwm_s1) > 0:
                    pwm_matrix = {}
                    pwm_s2 = pwm_s1.split(";")
                    for line in pwm_s2:
                        pwm_s3 = line.split(":")
                        if len( pwm_s3) > 1:
                            pwm_s4 = pwm_s3[1].split()
                            try:
                                length = 0
                                for value in pwm_s4:
                                    length += 1
                                    if not pwm_matrix.has_key( pwm_s3[0]):
                                        pwm_matrix[ pwm_s3[0]]=[]
                                    pwm_matrix[ pwm_s3[0]].append( int( value))
                            except ValueError, val_exce:
                                raise ParsingException( "BedSeqAlignmentStatsCommStruct.getAlignmentMotifs : Unable to get integer value for Motif '" + BedSeqAlignmentStatsCommStruct.MOTIF_PWM_ATT + "' attributes. From:\n\t---> " + str( val_exce))
                     
                    pwm = PWM()
                    pwm.matrix = pwm_matrix
                    pwm.totalLength = length
                    pwm.nbSequences = nb_species
                else:
                    pwm = None
                
                if start != None and end != None and name != None:
                    motif = Motif( start, end, name, pwm)
                    motif.offset = offset
                    motif.score = score
                    if consensus != None:
                        motif.consensus = consensus
                    if motif_id != None:
                        motif.id = motif_id
                    if strand != None:
                        motif.strand = strand
                    seqalign.addMotif( motif)
                else:
                    raise ParsingException( "BedSeqAlignmentStatsCommStruct.getAlignmentMotifs : The motifs of the alignment of '" + bedseq.toString() + "' is missing required attributes")
                    
            else:
                raise ParsingException( "BedSeqAlignmentStatsCommStruct.getAlignmentMotifs : The motifs of the alignment of '" + bedseq.toString() + "' contains an unauthorized element : '" + node_motif.tag.lower() +  "'")





    # #############################
    # METHODS TO WRITE THE XLM FILE
    # #############################


    # --------------------------------------------------------------------------------------
    # Convert the CommStruct to an ElementTree containing all the required
    # suitably organized information
    def convertCommStructToElementTree( self):
        
        data_element = Element( BedSeqAlignmentStatsCommStruct.DATA_TAG)
        
        self.addCommStructDataToElement( data_element)
        
        # Output the statistics data
        statistics_element = Element( BedSeqAlignmentStatsCommStruct.STATISTICS_TAG)
        data_element.append( statistics_element)
        for param_name in sorted( self.paramStatistics.keys()):
            paramstats_element = Element( BedSeqAlignmentStatsCommStruct.PARAM_TAG)
            statistics_element.append( paramstats_element)
            paramstats_element.attrib[ BedSeqAlignmentStatsCommStruct.PARAM_NAME_ATT] = param_name
            paramstats_element.attrib[ BedSeqAlignmentStatsCommStruct.PARAM_VALUE_ATT] = str( self.paramStatistics[ param_name])
            
        for motif_name in sorted( self.motifStatistics.keys()):
            motifstats_element = Element( BedSeqAlignmentStatsCommStruct.MOTIF_STATS_TAG)
            statistics_element.append( motifstats_element)
            
            motif_statistics = self.motifStatistics[ motif_name]
            motifstats_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_STATS_NAME_ATT] = motif_statistics.motifName
            motifstats_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_STATS_ID_ATT] = motif_statistics.motifID
            motifstats_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_STATS_FAMILY_ATT] = motif_statistics.motifFamily
            motifstats_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_STATS_CLASS_ATT] = motif_statistics.motifClass
            motifstats_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_STATS_TYPE_ATT] = motif_statistics.motifType
            motifstats_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_STATS_SIZE_ATT] = str( motif_statistics.motifSize)

            for att_name in sorted( motif_statistics.attributes.keys()):
                attribute_element = Element( BedSeqAlignmentStatsCommStruct.PARAM_TAG)
                motifstats_element.append( attribute_element)
                attribute_element.attrib[ BedSeqAlignmentStatsCommStruct.PARAM_NAME_ATT] = att_name
                attribute_element.attrib[ BedSeqAlignmentStatsCommStruct.PARAM_VALUE_ATT] = str( motif_statistics.getAttribute( att_name))
                
            
        # output the bedseq data
        for chrom in self.bedSequencesDict.keys():
            for bedseq in self.bedSequencesDict[chrom]:
                bedseq_element = Element( BedSeqAlignmentStatsCommStruct.BEDSEQ_TAG)
                data_element.append( bedseq_element)
                bedseq_element.attrib[ BedSeqAlignmentStatsCommStruct.BEDSEQ_SPECIES_ATT] = bedseq.species
                bedseq_element.attrib[ BedSeqAlignmentStatsCommStruct.BEDSEQ_CHROM_ATT] = bedseq.chromosom
                bedseq_element.attrib[ BedSeqAlignmentStatsCommStruct.BEDSEQ_START_ATT] = str( bedseq.indexStart)
                bedseq_element.attrib[ BedSeqAlignmentStatsCommStruct.BEDSEQ_END_ATT] = str( bedseq.indexEnd)
                bedseq_element.attrib[ BedSeqAlignmentStatsCommStruct.BEDSEQ_SCORE_ATT] = str( bedseq.score)
                if bedseq.referenceIndex != 0:
                    bedseq_element.attrib[ BedSeqAlignmentStatsCommStruct.BEDSEQ_MAX_PEAK_ATT] = str( bedseq.referenceIndex)
                if bedseq.id != "":
                    bedseq_element.attrib[ BedSeqAlignmentStatsCommStruct.BEDSEQ_ID_ATT] = str( bedseq.id)
                    
                if self.bedToMA.has_key( bedseq):
                    for alignment in self.bedToMA[ bedseq]:
                        alignment_element = Element( BedSeqAlignmentStatsCommStruct.ALIGNMENT_TAG)
                        alignment_element.attrib[ BedSeqAlignmentStatsCommStruct.ALIGNMENT_NAME_ATT] = alignment.name
                        bedseq_element.append( alignment_element)
                        sequences_element = Element( BedSeqAlignmentStatsCommStruct.SEQUENCES_TAG)
                        alignment_element.append( sequences_element)
                        for species in alignment.sequences.keys():
                            sequence_element = Element( BedSeqAlignmentStatsCommStruct.SEQUENCE_TAG)
                            sequences_element.append( sequence_element)
                            sequence_element.attrib[ BedSeqAlignmentStatsCommStruct.SEQUENCE_SPECIES_ATT] = species
                            sequence_element.attrib[ BedSeqAlignmentStatsCommStruct.SEQUENCE_TEXT_ATT] = "".join( alignment.sequences[species])
                        motifs_element = Element( BedSeqAlignmentStatsCommStruct.MOTIFS_TAG)
                        alignment_element.append( motifs_element)
                        alignment.motifs.sort( Motif.compare)
                        for motif in alignment.motifs:
                            motif_element = Element( BedSeqAlignmentStatsCommStruct.MOTIF_TAG)
                            motifs_element.append( motif_element)
                            motif_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_START_ATT] = str( motif.indexStart)
                            motif_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_END_ATT] = str( motif.indexEnd)
                            motif_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_NAME_ATT] = motif.name
                            motif_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_STRAND_ATT] = motif.strand
                            motif_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_OFFSET_ATT] = str( motif.offset)
                            motif_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_SCORE_ATT] = str( motif.score)
                            if motif.id != None and len( motif.id) > 0:
                                motif_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_ID_ATT] = motif.id
                            if motif.consensus != None and len( motif.consensus) > 0:
                                motif_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_CONSENSUS_ATT] = str( motif.consensus)
                            if motif.pwm != None:
                                pwm = ""
                                for letter in Constants.DNA_ALPHABET:
                                    pwm += letter +":"
                                    for value in motif.pwm.matrix[letter]:
                                        pwm += str( value) + " "
                                    pwm += ";"
                                motif_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_PWM_ATT] = pwm
                                motif_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_NBSPECIES_ATT] = str( motif.pwm.nbSequences)
                            else:
                                motif_element.attrib[ BedSeqAlignmentStatsCommStruct.MOTIF_PWM_ATT] = ""
                                
                       
        return data_element


    # --------------------------------------------------------------------------------------
    # Definition of the tags and attributes names used in the XML format
    
    DATA_TAG = "data"
    DATA_SPECIES_ATT = "species"
    DATA_GENE_ATT = "gene"
    DATA_TF_ATT = "tf"
    
    STATISTICS_TAG = "statistics"
    MOTIF_STATS_TAG = "motifstats"
    MOTIF_STATS_NAME_ATT = "name"
    MOTIF_STATS_ID_ATT = "id"
    MOTIF_STATS_FAMILY_ATT = "family"
    MOTIF_STATS_CLASS_ATT = "class"
    MOTIF_STATS_TYPE_ATT = "type"
    MOTIF_STATS_SIZE_ATT = "size"
    
    BEDSEQ_TAG = "bedseq"
    BEDSEQ_SPECIES_ATT = "species"
    BEDSEQ_CHROM_ATT = "chr"
    BEDSEQ_START_ATT = "start"
    BEDSEQ_END_ATT = "end"
    BEDSEQ_SCORE_ATT = "score"
    BEDSEQ_ID_ATT = "id"
    BEDSEQ_MAX_PEAK_ATT = "max"
    
    ALIGNMENT_TAG = "alignment"
    ALIGNMENT_NAME_ATT = "name"
    
    SEQUENCES_TAG = "sequences"
    SEQUENCE_TAG = "sequence"
    SEQUENCE_SPECIES_ATT = "species"
    SEQUENCE_TEXT_ATT = "text"

    MOTIFS_TAG = "motifs"
    MOTIF_TAG = "motif"
    MOTIF_START_ATT = "start"
    MOTIF_END_ATT = "end"
    MOTIF_NAME_ATT = "name"
    MOTIF_ID_ATT = "id"
    MOTIF_CONSENSUS_ATT = "consensus"
    MOTIF_NBSPECIES_ATT = "nbSpecies"
    MOTIF_PWM_ATT = "text"
    MOTIF_STRAND_ATT = "strand"
    MOTIF_OFFSET_ATT = "offset"
    MOTIF_SCORE_ATT = "score"
    
    PARAM_TAG = "param"
    PARAM_NAME_ATT = "name"
    PARAM_VALUE_ATT = "value"



    # #############################
    # #############################
    #     CLASS MotifStatistics
    # #############################
    # #############################


class MotifStatistics:
    
    MOTIF_HIT_SCORE = "HitScore"
    MOTIF_UNCOUNT = "Uncount"
    MOTIF_HYP_HIT_SCORE = "HypergeometricTestHitScore"
    MOTIF_HYP_UNCOUNT = "HypergometricUncount"
    MOTIF_HYP_PVALUE = "HypergeometricPValue"
    MOTIF_HYP_EVALUE = "HypergeometricEValue"
    MOTIF_CHI2 = "Chi2"
    MOTIF_CHI2_PVALUE = "Chi2PValue"
    
    MOTIF_DISTANCE_HISTOGRAM = "DistanceHistogram"
    MOTIF_DISTANCE_HISTOGRAM_GRAPH = "DistanceHistogramGraph"
    
    MOTIF_PEAK_SCORE_HISTOGRAM = "PeakScoreHistogram"
    MOTIF_PEAK_SCORE_HISTOGRAM_GRAPH ="PeakScoreHistogramGraph"
    
    MOTIF_COLOCATION_HISTOGRAM = "CoLocationHistogram"
    MOTIF_COLOCATION_HISTOGRAM_GRAPH = "CoLocationHistogramGraph"
    
    MOTIF_RANK = "Rank"
    MOTIF_FAMILY_RANK = "FamilyRank"
    MOTIF_RANK_IN_FAMILY = "RankInFamily"
    
    MOTIF_RATIO_HOMOLOCATION = "RatioHomoLocation"


    # --------------------------------------------------------------------------------------
    # 
    def __init__( self, name):
        
        self.motifName = name
        self.motifID = ""
        self.motifFamily = ""
        self.motifClass = ""
        self.motifType = ""
        self.motifSize = 0
        self.attributes = {}


    # --------------------------------------------------------------------------------------
    # Set the attribute to the given value
    def setAttribute(self, att_name, att_value):
        
        if att_name != None:
            self.attributes[ att_name] = att_value


    # --------------------------------------------------------------------------------------
    # Return the attribute value
    def getAttribute(self, att_name):
        
        if att_name in self.attributes.keys():
            return self.attributes[ att_name]
        else:
            return None


    # --------------------------------------------------------------------------------------
    # Return the attribute value as int
    def getAttributeAsint(self, att_name, mandatory = False):
        
        try:
            att_value = int( self.getAttribute( att_name))
            return att_value
        except (TypeError, ValueError), val_exce:
            if mandatory:
                raise ExecutionException( "MotifStatistics.getAttributeAsint : Unable to convert the value of attribute :'" + att_name + "'. From:\n\t---> " + str( val_exce))
            else:
                return 0
 

    # --------------------------------------------------------------------------------------
    # Return the attribute value as float
    def getAttributeAsfloat(self, att_name, mandatory = False):
        
        try:
            att_value = float( self.getAttribute( att_name))
            return att_value
        except (TypeError, ValueError), val_exce:
            if mandatory:
                raise ExecutionException( "MotifStatistics.getAttributeAsint : Unable to convert the value of attribute :'" + att_name + "'. From:\n\t---> " + str( val_exce))
            else:
                return 0
         

    # --------------------------------------------------------------------------------------
    # Return True if the attribute name is in the attribute dict
    def hasAttribute(self, att_name):
        
        return att_name in self.attributes.keys()
        
    # --------------------------------------------------------------------------------------
    # Compare two MotifStatistics (used for sorting) 
    @staticmethod
    def compare( stats1, stats2):
        
        # Compare the names
        if stats1.motifName == stats2.motifName:
            return 0

         # Get the values of hypergeometric and chi2 p-values
        hyp_pvalue_1 = stats1.getAttributeAsfloat( MotifStatistics.MOTIF_HYP_PVALUE)
        hyp_pvalue_2 = stats2.getAttributeAsfloat( MotifStatistics.MOTIF_HYP_PVALUE)    
        chi2_pvalue_1 = stats1.getAttributeAsfloat( MotifStatistics.MOTIF_CHI2_PVALUE)
        chi2_pvalue_2 = stats2.getAttributeAsfloat( MotifStatistics.MOTIF_CHI2_PVALUE)
        
        # Compare the products of log of hypergeometric and chi2 p-values to have the greatest first
        product_1 = 1.0
        if hyp_pvalue_1 > 0:
            product_1 = product_1 * math.log10(hyp_pvalue_1)
        else:
            product_1 = -400.0 * product_1
        if chi2_pvalue_1 > 0:
            product_1 = product_1 * math.log10(chi2_pvalue_1)
        else:
            product_1 = -400.0 * product_1
    
        product_2 = 1.0
        if hyp_pvalue_2 > 0:
            product_2 = product_2 * math.log10(hyp_pvalue_2)
        else:
            product_2 = -400.0 * product_2
        if chi2_pvalue_2 > 0:
            product_2 = product_2 * math.log10(chi2_pvalue_2)
        else:
            product_2 = -400.0 * product_2
    
        if product_1 != product_2:
            if product_1 < product_2:
                return 1
            else:
                return -1
    
        # Compare the hypergeometric p-value to have the littlest first
        if hyp_pvalue_1 != hyp_pvalue_2:
            if hyp_pvalue_1 < hyp_pvalue_2:
                return -1
            else:
                return 1
        
        # Compare the chi2 p-value to have the littlest first
        
        
        if chi2_pvalue_1 != chi2_pvalue_2:
            if chi2_pvalue_1 < chi2_pvalue_2:
                return -1
            else:
                return 1
    
        # Compare the chi2 to have the greatest first
        chi2_1 = stats1.getAttributeAsfloat( MotifStatistics.MOTIF_CHI2)
        chi2_2 = stats2.getAttributeAsfloat( MotifStatistics.MOTIF_CHI2)
        
        if chi2_1 == None:
            chi2_1 = 0
        if chi2_2 == None:
            chi2_2 = 0
        
        if chi2_1 != chi2_2:
            if chi2_1 > chi2_2:
                return -1
            else:
                return 1
            
        return 0


# eflag: FileType = Python2
