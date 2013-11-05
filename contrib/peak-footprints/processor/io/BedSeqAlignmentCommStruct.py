
from processor.io.CommStruct import CommStruct

from utils.exception.ParsingException import ParsingException
from utils.log.Log import Log
from utils.Constants import Constants

from xml.etree.ElementTree import parse
from xml.etree.ElementTree import Element, ElementTree

from common.BEDSequence import BEDSequence
from common.SequenceAlignment import SequenceAlignment
from common.Motif import Motif
from common.PWM import PWM

class BedSeqAlignmentCommStruct( CommStruct):


    # --------------------------------------------------------------------------------------
    def __init__( self):
        
        CommStruct.__init__( self)
        self.baseSpecies = ""
        # dicitonnary of the BED Sequences : key = chromosom / value = list of BED Sequence
        self.bedSequencesDict = {}
        # dicitonnary of the Multiple Alignments : key = BED Sequence / value = list of multiple alignments
        self.bedToMA = {}


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
        
        result = "BedSeqAlignmentCommStruct\n"
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
            return BedSeqAlignmentCommStruct.getCommStructFromXML( input_filepath)
        except ParsingException, par_exce:
            Log.log( "BedSeqAlignmentCommStruct.fromXMLFile : Unable to get CommStruct from XML file '" + input_filepath + "'. From:\n\t---> " + str( par_exce))
            return None


    # --------------------------------------------------------------------------------------
    # Write the CommStruct to the given XML file
    def toXMLFile( self, output_filepath):
        
        try:
            root_element = self.convertCommStructToElementTree()
            self.indent( root_element,  0)
            ElementTree( root_element).write( output_filepath)
        except IOError, exce:
            Log.log( "BedSeqAlignmentCommStruct.toXMLFile : Unable to write CommStruct to XML file. From:\n\t---> " + str( exce))
        except ParsingException, par_exce:
            Log.log( "BedSeqAlignmentCommStruct.toXMLFile : Unable to save CommStruct to XML file. From:\n\t---> " + str( par_exce))
    



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
            raise ParsingException( "BedSeqAlignmentCommStruct.getCommStructFromXML : Unable to open/close XML file '" + commstruct_filepath, "' : " + str( io_exce))
        
        comm_struct = BedSeqAlignmentCommStruct()
        
        comm_struct.getCommStructData( root_element)
        
        base_species = ""
        for node_bedseq in root_element:
            if node_bedseq.tag.lower() == BedSeqAlignmentCommStruct.BEDSEQ_TAG:
                bedseq = BedSeqAlignmentCommStruct.getBEDSequence( node_bedseq, comm_struct)
                base_species = bedseq.species
                for align_node in node_bedseq:
                    if align_node.tag.lower() == BedSeqAlignmentCommStruct.ALIGNMENT_TAG:
                        seqalign = BedSeqAlignmentCommStruct.getAlignment( align_node, bedseq, comm_struct)
                        for sub_node in align_node:
                            if sub_node.tag.lower() == BedSeqAlignmentCommStruct.SEQUENCES_TAG:   
                                BedSeqAlignmentCommStruct.getAlignmentSequences( sub_node, bedseq, seqalign)
                            elif sub_node.tag.lower() == BedSeqAlignmentCommStruct.MOTIFS_TAG:
                                BedSeqAlignmentCommStruct.getAlignmentMotifs( sub_node, bedseq, seqalign)
                            else:
                                raise ParsingException( "BedSeqAlignmentCommStruct.getCommStructFromXML : The alignment of '" + bedseq.toString() + "' contains an unauthorized element : '" + sub_node.tag.lower() +  "'")
                    else:
                        raise ParsingException( "BedSeqAlignmentCommStruct.getCommStructFromXML : The BED Sequence '" + bedseq.toString() + "' contains an unauthorized element : '" + align_node.tag.lower() +  "'")
            else:
                raise ParsingException( "BedSeqAlignmentCommStruct.getCommStructFromXML : The data contains an unauthorized element : '" + node_bedseq.tag.lower() +  "'")                    
        
        comm_struct.baseSpecies = base_species
        
        return comm_struct


    # --------------------------------------------------------------------------------------
    # Retrieve the BED sequence information, create a BEDSequence
    # and store it in the CommStruct
    @staticmethod
    def getBEDSequence( node_bedseq, comm_struct):
        
        try:
            species = CommStruct.getAttribute( node_bedseq, BedSeqAlignmentCommStruct.BEDSEQ_SPECIES_ATT)
            chrom = CommStruct.getAttribute( node_bedseq, BedSeqAlignmentCommStruct.BEDSEQ_CHROM_ATT)
            start = CommStruct.getAttributeAsint( node_bedseq, BedSeqAlignmentCommStruct.BEDSEQ_START_ATT)
            end = CommStruct.getAttributeAsint( node_bedseq, BedSeqAlignmentCommStruct.BEDSEQ_END_ATT)
            score = CommStruct.getAttributeAsint( node_bedseq, BedSeqAlignmentCommStruct.BEDSEQ_SCORE_ATT, False)
            max = CommStruct.getAttributeAsint( node_bedseq, BedSeqAlignmentCommStruct.BEDSEQ_MAX_PEAK_ATT, False)
            id = CommStruct.getAttribute( node_bedseq, BedSeqAlignmentCommStruct.BEDSEQ_ID_ATT, False)
        except ParsingException,  par_exce:
            raise ParsingException ( "BedSeqAlignmentCommStruct.getBEDSequence : Malformed BED Sequence - some attributes are not numbers. From:\n\t---> " + str( par_exce))
 
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
            raise ParsingException ( "BedSeqAlignmentCommStruct.getBEDSequence : Malformed BED Sequence - unable to retrieve sequence required attributes")


    # --------------------------------------------------------------------------------------
    # create the SequenceAlignment object
    @staticmethod
    def getAlignment( align_node, bedseq, comm_struct):
        
        name = CommStruct.getAttribute( align_node, BedSeqAlignmentCommStruct.ALIGNMENT_NAME_ATT)
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
            if node_sequence.tag.lower() == BedSeqAlignmentCommStruct.SEQUENCE_TAG:
                species = CommStruct.getAttribute( node_sequence, BedSeqAlignmentCommStruct.SEQUENCE_SPECIES_ATT)
                text = list( CommStruct.getAttribute( node_sequence, BedSeqAlignmentCommStruct.SEQUENCE_TEXT_ATT))
                
                if species != None and text != None:
                    seqalign.addSequence( species, text)
                else:
                    raise ParsingException( "BedSeqAlignmentCommStruct.getAlignmentSequences : A sequence of the alignment of '" + bedseq.toString() + "' is missing required attributes")
                    
            else:
                raise ParsingException( "BedSeqAlignmentCommStruct.getAlignmentSequences : The sequences of the alignment of '" + bedseq.toString() + "' contains an unauthorized element : '" + node_sequence.tag.lower() +  "'")


    # --------------------------------------------------------------------------------------
    # Retrieve the list of motifs of the alignment
    # and store it into the given SequenceAlignment
    @staticmethod
    def getAlignmentMotifs( sub_node, bedseq, seqalign):
        
        for node_motif in sub_node:
            if node_motif.tag.lower() == BedSeqAlignmentCommStruct.MOTIF_TAG:
                
                start = CommStruct.getAttributeAsint( node_motif, BedSeqAlignmentCommStruct.MOTIF_START_ATT)
                end = CommStruct.getAttributeAsint( node_motif, BedSeqAlignmentCommStruct.MOTIF_END_ATT)
                name = CommStruct.getAttribute( node_motif, BedSeqAlignmentCommStruct.MOTIF_NAME_ATT)
                id = CommStruct.getAttribute( node_motif, BedSeqAlignmentCommStruct.MOTIF_ID_ATT, False)
                consensus = CommStruct.getAttribute( node_motif, BedSeqAlignmentCommStruct.MOTIF_CONSENSUS_ATT, False)
                nb_species = CommStruct.getAttributeAsint( node_motif, BedSeqAlignmentCommStruct.MOTIF_NBSPECIES_ATT, False)
                strand = CommStruct.getAttribute( node_motif, BedSeqAlignmentCommStruct.MOTIF_STRAND_ATT, False)
                offset = CommStruct.getAttributeAsint( node_motif, BedSeqAlignmentCommStruct.MOTIF_OFFSET_ATT)
                
                # Retrieve the PWM of the motif

                pwm_s1 = CommStruct.getAttribute( node_motif, BedSeqAlignmentCommStruct.MOTIF_PWM_ATT)
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
                                raise ParsingException( "BedSeqAlignmentCommStruct.getAlignmentMotifs : Unable to get integer value for Motif '" + BedSeqAlignmentCommStruct.MOTIF_PWM_ATT + "' attributes. From:\n\t---> " + str( val_exce))
                     
                    pwm = PWM()
                    pwm.matrix = pwm_matrix
                    pwm.totalLength = length
                    pwm.nbSequences = nb_species
                else:
                    pwm = None
                
                if start != None and end != None and name != None:
                    motif = Motif( start, end, name, pwm)
                    motif.offset = offset
                    if consensus != None:
                        motif.consensus = consensus
                    if id != None:
                        motif.id = id
                    if strand != None:
                        motif.strand = strand
                    seqalign.addMotif( motif)
                else:
                    raise ParsingException( "BedSeqAlignmentCommStruct.getAlignmentMotifs : The motifs of the alignment of '" + bedseq.toString() + "' is missing required attributes")
                    
            else:
                raise ParsingException( "BedSeqAlignmentCommStruct.getAlignmentMotifs : The motifs of the alignment of '" + bedseq.toString() + "' contains an unauthorized element : '" + node_motif.tag.lower() +  "'")





    # #############################
    # METHODS TO WRITE THE XLM FILE
    # #############################


    # --------------------------------------------------------------------------------------
    # Convert the CommStruct to an ElementTree containing all the required
    # suitably organized information
    def convertCommStructToElementTree( self):
        
        data_element = Element( BedSeqAlignmentCommStruct.DATA_TAG)
        
        self.addCommStructDataToElement( data_element)
        
        for chrom in sorted( self.bedSequencesDict.keys()):
            for bedseq in self.bedSequencesDict[chrom]:
                bedseq_element = Element( BedSeqAlignmentCommStruct.BEDSEQ_TAG)
                data_element.append( bedseq_element)
                bedseq_element.attrib[ BedSeqAlignmentCommStruct.BEDSEQ_SPECIES_ATT] = bedseq.species
                bedseq_element.attrib[ BedSeqAlignmentCommStruct.BEDSEQ_CHROM_ATT] = bedseq.chromosom
                bedseq_element.attrib[ BedSeqAlignmentCommStruct.BEDSEQ_START_ATT] = str( bedseq.indexStart)
                bedseq_element.attrib[ BedSeqAlignmentCommStruct.BEDSEQ_END_ATT] = str( bedseq.indexEnd)
                bedseq_element.attrib[ BedSeqAlignmentCommStruct.BEDSEQ_SCORE_ATT] = str( bedseq.score)
                if bedseq.referenceIndex != 0:
                    bedseq_element.attrib[ BedSeqAlignmentCommStruct.BEDSEQ_MAX_PEAK_ATT] = str( bedseq.referenceIndex)
                if bedseq.id != "":
                    bedseq_element.attrib[ BedSeqAlignmentCommStruct.BEDSEQ_ID_ATT] = str( bedseq.id)
                    
                if self.bedToMA.has_key( bedseq):
                    for alignment in self.bedToMA[ bedseq]:
                        alignment_element = Element( BedSeqAlignmentCommStruct.ALIGNMENT_TAG)
                        alignment_element.attrib[ BedSeqAlignmentCommStruct.ALIGNMENT_NAME_ATT] = alignment.name
                        bedseq_element.append( alignment_element)
                        sequences_element = Element( BedSeqAlignmentCommStruct.SEQUENCES_TAG)
                        alignment_element.append( sequences_element)
                        for species in alignment.sequences.keys():
                            sequence_element = Element( BedSeqAlignmentCommStruct.SEQUENCE_TAG)
                            sequences_element.append( sequence_element)
                            sequence_element.attrib[ BedSeqAlignmentCommStruct.SEQUENCE_SPECIES_ATT] = species
                            sequence_element.attrib[ BedSeqAlignmentCommStruct.SEQUENCE_TEXT_ATT] = "".join( alignment.sequences[species])
                        motifs_element = Element( BedSeqAlignmentCommStruct.MOTIFS_TAG)
                        alignment_element.append( motifs_element)
                        alignment.motifs.sort( Motif.compare)
                        for motif in alignment.motifs:
                            motif_element = Element( BedSeqAlignmentCommStruct.MOTIF_TAG)
                            motifs_element.append( motif_element)
                            motif_element.attrib[ BedSeqAlignmentCommStruct.MOTIF_START_ATT] = str( motif.indexStart)
                            motif_element.attrib[ BedSeqAlignmentCommStruct.MOTIF_END_ATT] = str( motif.indexEnd)
                            motif_element.attrib[ BedSeqAlignmentCommStruct.MOTIF_NAME_ATT] = motif.name
                            motif_element.attrib[ BedSeqAlignmentCommStruct.MOTIF_STRAND_ATT] = motif.strand
                            motif_element.attrib[ BedSeqAlignmentCommStruct.MOTIF_OFFSET_ATT] = str( motif.offset)
                            if motif.id != None and len( motif.id) > 0:
                                motif_element.attrib[ BedSeqAlignmentCommStruct.MOTIF_ID_ATT] = motif.id
                            if motif.consensus != None and len( motif.consensus) > 0:
                                motif_element.attrib[ BedSeqAlignmentCommStruct.MOTIF_CONSENSUS_ATT] = str( motif.consensus)
                            if motif.pwm != None:
                                pwm = ""
                                for letter in Constants.DNA_ALPHABET:
                                    pwm += letter +":"
                                    for value in motif.pwm.matrix[letter]:
                                        pwm += str( value) + " "
                                    pwm += ";"
                                motif_element.attrib[ BedSeqAlignmentCommStruct.MOTIF_PWM_ATT] = pwm
                                motif_element.attrib[ BedSeqAlignmentCommStruct.MOTIF_NBSPECIES_ATT] = str( motif.pwm.nbSequences)
                            else:
                                motif_element.attrib[ BedSeqAlignmentCommStruct.MOTIF_PWM_ATT] = ""
                       
        return data_element


    # --------------------------------------------------------------------------------------
    # Definition of the tags and attributes names used in the XML format
    
    DATA_TAG = "data"
    DATA_SPECIES_ATT = "species"
    DATA_GENE_ATT = "gene"
    DATA_TF_ATT = "tf"
    
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
    
    PARAM_TAG = "param"
    PARAM_NAME_ATT = "name"
    PARAM_VALUE_ATT = "value"

# eflag: FileType = Python2
