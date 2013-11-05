
from processor.Processor import Processor
from processor.io.BedSeqCommStruct import BedSeqCommStruct

from utils.log.Log import Log
from utils.parser.BEDParser import BEDParser

class MultipleBEDProcessor( Processor):
    
    INPUT_BED_FILE_PARAM = "BEDFiles"
    SPECIES_PARAM = "Species"
    
    
    # --------------------------------------------------------------------------------------
    def __init__(self):
        Processor.__init__( self)
    
    # --------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as input 
    # (None if no input CommStruct)
    @staticmethod
    def getInputCommStructClass():
        
        return None


    # --------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as output
    # (None if no output CommStruct)
    @staticmethod
    def getOutputCommStructClass():
        
        return ( BedSeqCommStruct,  )
        

    # --------------------------------------------------------------------------------------
    # Returns a name that will be used as display name in the user friendly outputs
    @staticmethod
    def getDisplayName():
        
        return "Parsing of Multiple BED format input file"



    # --------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return ( MultipleBEDProcessor.INPUT_BED_FILES_PARAM, MultipleBEDProcessor.SPECIES_PARAM)
        
        
        
    # --------------------------------------------------------------------------------------
    # Execute the processor
    def execute(self, input_comm_structs):
        
        # Retrieve the processor parameters
        bed_filepaths = self.getParameter( MultipleBEDProcessor.INPUT_BED_FILES_PARAM)
        species = self.getParameter( MultipleBEDProcessor.SPECIES_PARAM)
        
        bed_files_list = bed_filepaths.split( ";")
        
        bedseq_dictionnaries = []
        for bed_filepath in bed_files_list:
            self.bedseq_dictionnaries.add( BEDParser.getBEDSequenceDictionnary( species, bed_filepath));
        
        # Build the list of all chromosom present in only one list
        chrom_list = []
        for bedseq_dic in bedseq_dictionnaries:
            for chrom in bedseq_dic.keys():
                if not chrom in chrom_list:
                    chrom_list.append( chrom)
        
        # Launch the peak classification for each chromosom
        for chrom in chrom_list:
            bedseq_lists = []
            for bedseq_dic in bedseq_dictionnaries:
                bedseq_lists.append( bedseq_dic.get( chrom))
            self.classifySequences( bedseq_lists)
    
    
    
    # --------------------------------------------------------------------------------------
    # Search for peaks multi-species intersections and for isolated peaks
    def classifySequences(self, bedseq_lists):
        
        global_seq_list = []
        for index in len(bedseq_lists):
            for bedseq in bedseq_lists[index]:
                global_seq_list.append( IndexedSequence( bedseq, index))
        
        global_seq_list.sort( IndexedSequence.compare)
        
        
        
        
            

    # #############################
    # #############################
    #     CLASS IndexedSequence
    # #############################
    # #############################


class IndexedSequence:
    
    def __init__( self, sequence,  index):
        
        self.sequence = sequence
        self.index = index
    
    # --------------------------------------------------------------------------------------
    # Compare two IndexedSequence (used for sorting) 
    @staticmethod
    def compare( iseq1, iseq2):
        
        start1 = iseq1.sequence.indexStart
        start2 = iseq2.sequence.indexStart
        
        if start1 < start2:
            return 1
        elif start1 > start2:
            return -1
        else:
            end1 = iseq1.sequence.indexEnd
            end2 = iseq2.sequence.indexEnd
            if end1 < end2:
                return 1
            elif end1 > end2:
                return -1
        
        return 0
            


    # #############################
    # #############################
    #     CLASS IndexedSequence
    # #############################
    # #############################


class SequenceIntersection:
    
    def __init__( self, category_nb):
        
        self.sequences = [None] * category_nb
        self.categoryNumber = category_nb
        self.indexStart = -1
        self.indexEnd = -1
            
    # --------------------------------------------------------------------------------------
    # Set the given sequence at the given category
    def setSequence(self, category,  sequence):
        
        result = None
        
        self.sequences[ category] = sequence
        
        if self.isGlobalIntersection():
            if self.indexStart < 0:
                self.indexStart = self.getCurrentStartIndex()
            self.indexEnd = self.getCurrentEndIndex()
        else:
            if self.indexStart > 0:
                result = (self.indexStart, self.indexEnd)
                self.indexStart = -1;
                self.indexEnd = -1;
                self.sequences = [None] * category_nb
                self.sequences[ category] = sequence
            else:
                
        
        return result


    # --------------------------------------------------------------------------------------
    # Remove the sequence at the given category
    def removeSequence(self, category):
        
        self.sequences[ category] = None


    # --------------------------------------------------------------------------------------
    # Return the lowest start index through all the listed sequences
    def getCurrentStartIndex(self):
        
        start = -1
        for sequence in self.sequences:
            if start < 0 or sequence.indexStart < start:
                start = sequence.indexStart
        
        return start


    # --------------------------------------------------------------------------------------
    # Return the highest end index through all the listed sequences        
    def getCurrentEndIndex(self):
        
        end = -1
        for sequence in self.sequences:
            if end < 0 or sequence.indexEnd > end:
                end = sequence.indexEnd
        
        return end

    # --------------------------------------------------------------------------------------
    # Return true if the intersection contains the maximal number of sequences  
    def isGlobalIntersection(self):
        
        count = 0
        min_end = -1
        max_start = -1
        for sequence in self.sequences:
            if sequence != None:
                count += 1
                if min_end < 0 or sequence.indexEnd < min_end:
                    min_end = sequence.indexEnd
                if max_start < 0 or sequence.indexStart > max_start:
                    max_start = sequence.indexStart
                
        if count == self.categoryNumber and max_start < min_end:
            return True
        
        return False
        


    # --------------------------------------------------------------------------------------
    # If the intersection has only one sequence returns its category, else return -1     
    def isUnique(self):
        
        count = 0
        category = -1
        for cat in len( self.sequences):
            if self.sequences[ cat] != None:
                count += 1
                category = cat
        
        if count == 1:
            return category
        else:
            return -1
        
        
