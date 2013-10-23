
from processor.Processor import Processor
from processor.io.BedSeqCommStruct import BedSeqCommStruct
from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct

from manager.ProgressionManager import ProgressionManager

from utils.FileUtils import FileUtils
from utils.Constants import Constants
from utils.RSATUtils import RSATUtils
from utils.log.Log import Log
from utils.exception.ParsingException import ParsingException
from utils.exception.ExecutionException import ExecutionException
from utils.SequenceUtils import SequenceUtils


from common.SequenceAlignment import SequenceAlignment

import Queue, threading, time, os, shutil

# This processor takes as input a list of sequence in BED format and as sources one or several
# multiple alignment files in MAF format. The species of the BED sequences must be the
# same as the one used by the MAF format as reference species.
# For each BED Sequence, the processor searches through the blocks of the MAF files the ones
# intersecting it. One all the corresponding block are found, the text of the BED sequence is
# recomposed using the text in the blocks and the corresponding alignment is so created.

# Parameters:
#   MAFFile : the path to the MAF file(s)
#   SpecializedFile : defines if the MAF files are specialized or not by chromosom
#   ReferenceSpecies : the species taken as reference in the MAF file
#   DesiredSpeciesList (optional): the species to take into account in the block of the MAF file. if not set, all the
#                                  species are used.
#   ThreadNumber (optional) : Option to activate the multi-threading during parsing. If a number greater than one is
#                             given, the parsing will launch as many thread, each thread will parse one MAF file.
#                             By default, this number is 1, meaning no multi-threading.
#   Keep gaps (optional) : Option to activate the conservation of position in the MSA that correspond to a gap 
#                          in the reference species sequence. By default, this option is set to False.


class MAFProcessor( Processor):
    
    # Parameter names for processor
    INPUT_MAF_FILE_PARAM = "MAFFile"
    SPECIALIZED_MAF_FILE_PARAM  = "SpecializedFile"
    REFERENCE_SPECIES_PARAM = "ReferenceSpecies"
    DESIRED_SPECIES_LIST_PARAM = "DesiredSpeciesList"
    THREAD_NUMBER_PARAM = "ThreadNumber"
    KEEP_GAPS = "KeepGaps"
    
    # Parameter for multi-threading
    THREAD_CHECK_DELAY = 2
    
    # Parameters for MAF file parsing
    _lineType_col = 0
    _speciesChrom_col = 1
    _startindex_col = 2
    _textlength_col = 3
    _strand_col = 4
    _source_size_col = 5
    _text_col = 6


    # --------------------------------------------------------------------------------------
    def __init__( self):
        Processor.__init__( self)
        self.desiredSpeciesList = []
        self.referenceSpecies = None
        self.bedSequencesDict = None
        self.mafBlockDic = {}
        self.parsedSpeciesList = []
        self.threadLock = None


    # --------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as input 
    # (None if no input CommStruct)
    @staticmethod
    def getInputCommStructClass():
        
        return ( BedSeqCommStruct, )


    # --------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as output
    # (None if no output CommStruct)
    @staticmethod
    def getOutputCommStructClass():
        
        return ( BedSeqAlignmentStatsCommStruct, )



    # --------------------------------------------------------------------------------------
    # Returns a name that will be used as display name in the user friendly outputs
    @staticmethod
    def getDisplayName():
        
        return "Research of matching multiple sequence alignments"



    #------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return ( MAFProcessor.INPUT_MAF_FILE_PARAM, MAFProcessor.REFERENCE_SPECIES_PARAM)
        


    # --------------------------------------------------------------------------------------
    # Execute the processor
    def execute( self, input_commstructs):
        
        if input_commstructs == None or len( input_commstructs) == 0:
            raise ExecutionException( "MAFProcessor.execute : No inputs")
        
        input_commstruct = input_commstructs[0]

        # retrieve processor parameters
        source_maffile = self.getParameter( MAFProcessor.INPUT_MAF_FILE_PARAM)
        
        specialized_file_line = self.getParameter( MAFProcessor.SPECIALIZED_MAF_FILE_PARAM, False)
        if specialized_file_line == None:
            specialized_file = False
        else:
            specialized_file = ( specialized_file_line.lower() == "true")
            
        desired_species_line = self.getParameter( MAFProcessor.DESIRED_SPECIES_LIST_PARAM, False)
        if desired_species_line != None:
            self.desiredSpeciesList = desired_species_line.split()
            
        self.referenceSpecies = self.getParameter( MAFProcessor.REFERENCE_SPECIES_PARAM)
        
        thread_number = self.getParameterAsint( MAFProcessor.THREAD_NUMBER_PARAM, False)
        if thread_number == None or thread_number < 0:
            thread_number = 1
            
        keep_gaps = self.getParameter( MAFProcessor.KEEP_GAPS, False)
        if keep_gaps == None:
            keep_gaps = False
            
        
        # Retrieve BED sequences from the input CommStruct
        self.bedSequencesDict = input_commstruct.bedSequencesDict
        
        if self.bedSequencesDict == None or len( self.bedSequencesDict) == 0:
            raise ExecutionException( "MAFProcessor.execute : No BEDSequence provided as input")
        
        # Look for MAF files to parse
        maf_file_list = FileUtils.getFileList( source_maffile, "maf", self.referenceSpecies)
        if maf_file_list == None:
            raise ExecutionException( "MAFProcessor.execute : The path '" + source_maffile + "' does not point to a MAF file or a directory containing MAF files and does not contain a subdirectory '" + self.referenceSpecies + "' containing MAF files.")
        
        # Retrieve the alignement blocks of the MAF files corresponding to each BED sequence,
        # managing the parsing of MAF file list according to the chosen number of threads
        count_parsed_files = 0
        ProgressionManager.setTaskProgression( "Parsing MAF Files", self.component, 0.0)
        if thread_number > 1:
            # store the list of MAF file in a queue
            file_queue = Queue.Queue( len( maf_file_list))
            for file in maf_file_list:
                file_queue.put( file)
                
            # Create a Thread lock used in several methods to avoid possible Thread conflicts
            self.threadLock = threading.Lock()
                
            # Launch the firsts Threads
            thread_list =[]
            while not file_queue.empty() and len( thread_list) < thread_number:
                self.startNewThread( file_queue, specialized_file, thread_list)
     
            # Manage the Threads list in order to empty the file queue
            while not file_queue.empty() or len( thread_list) > 0:
                for threads in thread_list:
                    if not threads.is_alive():
                        thread_list.remove( threads)
                        count_parsed_files += 1
                        ProgressionManager.setTaskProgression( "Parsing MAF Files", self.component, count_parsed_files/float( len( maf_file_list)))
                        self.startNewThread( file_queue, specialized_file, thread_list)
                time.sleep( MAFProcessor.THREAD_CHECK_DELAY)
        else:
            # Parse the files in the MAF file list
            for file in maf_file_list:
                self.parseFile(file, specialized_file)
                count_parsed_files += 1
                ProgressionManager.setTaskProgression( "Parsing MAF Files", self.component, count_parsed_files/float( len( maf_file_list)))
            
        # Assign the whole list of parsed species if a desired list was not set
        if len( self.desiredSpeciesList) == 0:
            self.desiredSpeciesList = self.parsedSpeciesList
        
        # create the output CommStruct
        output_commstruct = BedSeqAlignmentStatsCommStruct()
        output_commstruct.processorName = self.component.processorName
        output_commstruct.baseSpecies = input_commstruct.baseSpecies
        output_commstruct.bedSequencesDict = input_commstruct.bedSequencesDict
        output_commstruct.paramStatistics = input_commstruct.paramStatistics
        
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.REFERENCE_SPECIES] = self.referenceSpecies
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.ALIGNED_SPECIES] = ", ".join( self.desiredSpeciesList)

        # Compose the MSA corresponding to each BED sequence from the MAF blocks
        ProgressionManager.setTaskProgression( "Building MSA", self.component, 0.0)
        count = 0
        total_number_bed = len( self.mafBlockDic.keys())
        min_size = 100000000
        max_size = -1
        total_size = 0
        msa_lenghts = []
        for bed_sequence in self.mafBlockDic.keys():
            count += 1
            if count % 100 == 0:
                ProgressionManager.setTaskProgression( "Building MSA", self.component, count / float( total_number_bed))
            alignment = self.composeSequenceAlignment( bed_sequence, keep_gaps)
            align_length = alignment.totalLength
            msa_lenghts.append( align_length)
            if align_length < min_size:
                min_size = align_length
            if align_length > max_size:
                max_size = align_length
            total_size += align_length
            
            output_commstruct.addSequenceAlignment( bed_sequence, alignment)
        
        ProgressionManager.setTaskProgression( "Building MSA", self.component, 1.0)

        mean_size = (int) (total_size / float( total_number_bed))

        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.MSA_NUMBER] = count
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.MSA_MIN_SIZE] = min_size
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.MSA_MAX_SIZE] = max_size
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.MSA_MEAN_SIZE] = mean_size
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.MSA_TOTAL_SIZE] = total_size
        Log.trace( "MAFProcessor.execute : Total number of BEDsequence with associated MSA = " + str( count))
        Log.trace( "MAFProcessor.execute : Minimal size of BEDsequence with associated MSA = " + str( min_size))
        Log.trace( "MAFProcessor.execute : Maximal size of BEDsequence with associated MSA = " + str( max_size))
        Log.trace( "MAFProcessor.execute : Mean size of BEDsequence with associated MSA = " + str( mean_size))
        Log.trace( "MAFProcessor.execute : Total size of BEDsequence with associated MSA = " + str( total_size))
        
        # Output the MSA lengths histogram and graph
        self.outputMSALenghtHistogram( msa_lenghts, output_commstruct)
        
        
        return output_commstruct



    # --------------------------------------------------------------------------------------
    # Create and start a new Thread to parse the next file in the queue
    def startNewThread(self, file_queue, specialized_file, thread_list):
        
        if not file_queue.empty():
            file = file_queue.get()
            my_thread = threading.Thread( None, self.parseFile, file, ( file, specialized_file, ))
            thread_list.append( my_thread)
            Log.trace( "MAFProcessor.startNewThread : Starting new thread to parse file : '" + file + "'. Number of active Thread = " + str( len( thread_list)))
            my_thread.start()



    # --------------------------------------------------------------------------------------
    # For each BEDSequence, try to recompose the corresponding DNA sequence from the
    # reference sequences in the associated MAFblocks and build at the same time the
    # multiple alignment
    def composeSequenceAlignment( self, bed_sequence, keep_gaps):
        
        final_seq_align = SequenceAlignment()
        final_seq_align.name =  bed_sequence.name + "_1"
        final_seq_align.initializeWithDots( bed_sequence.getLength(), self.referenceSpecies, self.desiredSpeciesList)
        for maf_block in self.mafBlockDic[ bed_sequence]:
            ref_dna_seq = maf_block.sequences[ self.referenceSpecies]
            # Determine the start and end indexes (related to the string) of the intersection between
            # the BED sequence and the reference DNA sequence
            dna_index_start = max( bed_sequence.indexStart - ref_dna_seq.indexStart,  0)
            dna_index_end = min( bed_sequence.indexEnd - ref_dna_seq.indexStart, ref_dna_seq.textLength)
            
            # Compute the start and end index the block must be placed to the final alignment text
            bed_index_start = max( ref_dna_seq.indexStart - bed_sequence.indexStart,  0)
            bed_index_end = bed_index_start + (dna_index_end - dna_index_start)
            
            # Modify the start and end indexes according to the number of "-" present
            # in the sequence text to catch all the required letters
            current_index = 0
            while current_index <= dna_index_start:
                if ref_dna_seq.text[current_index] == Constants.SEQUENCE_INSERTION_CHAR:
                    dna_index_start += 1
                    dna_index_end += 1
                current_index += 1
            
            while current_index < dna_index_end:
                if ref_dna_seq.text[current_index] == Constants.SEQUENCE_INSERTION_CHAR:
                    dna_index_end += 1
                current_index += 1
            
            # Modify the bed_index_start and bed_index_end according to 
            # the number of "-" present in the sequence text
            # of the reference species currently present in the alignment

            count = 0
            indice = 0
            current = 0
            while indice <= bed_index_start:
                if final_seq_align.sequences[ self.referenceSpecies][ current] == Constants.SEQUENCE_INSERTION_CHAR:
                    count += 1
                else:
                    indice += 1
                current += 1
            bed_index_start += count
            bed_index_end += count

            # Insert the block of sequence of each species in the final alignement 
            for species in maf_block.sequences.keys():
                final_seq_align.insertSequenceBlock( species, bed_index_start, bed_index_end, maf_block.sequences[species].text[ dna_index_start : dna_index_end]   )
        
        # Cross algorithm to test the composition of the sequence
        #self.testFinalMSA( final_seq_align, self.mafBlockDic[ bed_sequence], bed_sequence)
        
        final_seq_align.finalizeSequences( keep_gaps)
            
        return final_seq_align


    
    # --------------------------------------------------------------------------------------
    # Parse the file at the given path if it is a MAF file.
    # If the file contains only information for a specific chromosom,
    # the 'is_chrom_file' boolean must be set to true in order to bypass the whole file
    # when no BED Sequence corresponds to this chromosom
    # 
    # Returns a dictionnary with the BEDSequence as key and the list of corresponding MAFBlock as value
    def parseFile(self, file_name, is_chrom_file):
        
        try:
            input_file = open( file_name, 'r')

            # Verify if the token '##maf' indicating a MAF file is found in the first lines
            is_maf_file = False
            while 1:
                line = input_file.readline()
                if len( line) == 0:
                    break
                elif not line.isspace():
                    tokens = line.split()
                    if tokens != None and len( tokens) > 0 and tokens[0] == "##maf":
                        is_maf_file = True
                        break
            
            # if it is a maf file, verify if an index file exists
            if is_maf_file == True:
                indexed = False
                try:
                    index_path =  file_name + "index"
                    input_index_file = open( index_path,  "r")
                    indexed = True
                    
                except IOError:
                    pass
                
                if indexed == True:
                    Log.trace( "MAFProcessor.parseFile : parsing file '" + file_name + "' using index '" + index_path + "'")
                    self.parseBlockListWithIndex( input_index_file, input_file)
                    self.closeFile( input_index_file)
                else:
                    Log.trace( "MAFProcessor.parseFile : parsing file '" + file_name + "'")
                    self.parseBlockListWithoutIndex( input_file, is_chrom_file)
                    
                self.closeFile( input_file)
                return
                
            else:
                self.closeFile( input_file)
                raise ParsingException( "MAFProcessor.parseFile : The file '" + file_name + "' is not a MAF file")
        except IOError,  io_exec:
            raise ParsingException( "MAFProcessor.parseFile : Enable to open file '" + file_name + "'. From:\n\t---> " + str(io_exec))


    # --------------------------------------------------------------------------------------
    # Parse the alignment block of the file linearly
    def parseBlockListWithoutIndex(self, input_file, is_chrom_file):
        
        # search for the next line starting with 'a' (meaning new alignment lbock)
        counter = 0
        while 1:
            line = input_file.readline()
            if len( line) == 0:
                break
            elif not line.isspace():
                tokens = line.split()
                if tokens != None and len( tokens) > 0 and tokens[ MAFProcessor._lineType_col] == "a":
                    counter += 1
                    if counter % 100000 == 0:
                        Log.trace( "MAFIndexerProcessor.execute : Number of MSA already parsed : " + str( counter))
                    parsed = self.parseBlock( input_file)
                    if not parsed and is_chrom_file:
                        return


    # --------------------------------------------------------------------------------------
    # Search for the right alignment blocks of the file using the index file
    def parseBlockListWithIndex( self, index_file, input_file):
        
        is_chrom_file = False
        ordered = False
        spec_chrom = None
        
        # Read the index file header to know if the file is chromosom specialized and ordered
        while 1:
            line = index_file.readline()
            if len( line) == 0:
                Log.log( "MAFProcessor.parseBlockListWithIndex : index file '" + index_file.name + "' has no header line : skipping it")
                return
            else:
                tokens = line.split()
                if tokens != None and tokens[0] == Constants.COMMENT_CHAR:
                    if len( tokens) > 1 and tokens[1] != Constants.MIXED:
                        is_chrom_file = True
                        spec_chrom = tokens[1]
                        if len( tokens) > 2 and tokens[2] == Constants.ORDERED:
                            ordered = True
                break

        # If the file is specialized by chromosom, get once for all the bed sequences concerned 
        # by this species and chromosom
        if is_chrom_file == True:
            bed_sequences = self.getAssociatedBEDSequences( spec_chrom)
            if bed_sequences == None or len( bed_sequences) == 0:
                Log.info( "MAFProcessor.parseBlockListWithIndex : No BED sequences matching for file :" + index_file.name)
                return
        else:
            bed_sequences = None

        # If file is ordered, compute the peaks extremum in order to optimize the parsing
        if ordered == True:
            min_start = 1000000000
            max_end = 0
            for bed_sequence in bed_sequences:
                if bed_sequence.indexStart < min_start:
                    min_start = bed_sequence.indexStart
                if bed_sequence.indexEnd > max_end:
                    max_end = bed_sequence.indexEnd
        
        # Parse the index file
        while 1:
            line = index_file.readline()
            if len( line) == 0:
                break
            else:
                tokens = line.split()
                if tokens != None and len( tokens) == 4:
                    # retrieve the index information
                    spec_chrom = tokens[ 0]
                    start = self.getIntValue( tokens[ 1])
                    end = self.getIntValue( tokens[ 2])
                    position = self.getIntValue( tokens[ 3])
                    
                    if ordered == True:
                        # If the file is ordered and the indexes are less than the BED indexes, skip the line
                        if end <= min_start:
                            continue 
                        # If the file is ordered and the indexes are greater than the BED indexes, skip the file
                        elif start >= max_end:
                            break
                        # If the indexes are at least in one of the BED sequences index range,
                        # the corresponding MSA block is parsed
                        else:
                            for bed_sequence in bed_sequences:
                                if end > bed_sequence.indexStart and start < bed_sequence.indexEnd:
                                    input_file.seek( position, 0)
                                    result = self.parseBlock( input_file, True)
                                    if result == False:
                                        raise ExecutionException( "MAFFile.parseBlockListWithIndex : Indexed MSA block seems not correct. You should have not updated indexes. Please see logs for more information")
                                    break
                    else:
                        # If the file is not chromosom specialized, the bed sequence list must be
                        # retrieve for each new index
                        if is_chrom_file == False:
                            bed_sequences = self.getAssociatedBEDSequences( spec_chrom)
                            if (bed_sequences == None or len( bed_sequences) == 0):
                                continue
                        
                        # If the indexes are at least in one of the BED sequences index range,
                        # the corresponding MSA block is parsed
                        for bed_sequence in bed_sequences:
                            if end > bed_sequence.indexStart and start < bed_sequence.indexEnd:
                                input_file.seek( position, 0)
                                result = self.parseBlock( input_file, True)
                                if result == False:
                                    Log.log( "MAFFile.parseBlockListWithIndex : Indexed MSA block seems not correct. You should have not updated indexes")
                                    raise ExecutionException( "MAFFile.parseBlockListWithIndex : Indexed MSA block seems not correct. You should have not updated indexes. Please, see logs for more information")
                                break



    # --------------------------------------------------------------------------------------
    # Parse a alignment block
    def parseBlock( self, input_file,  indexed = False):
                
        new_block = None
        
        #Search for the first sequence line of the block and verify if the block match
        # with at least one of the BED sequence
        while 1:
            line = input_file.readline()
            if len( line) == 0:
                break
            # Check if the line is not void
            elif not line.isspace():
                tokens = line.split()
                # Check if the line contains enough tokens
                if tokens != None and len( tokens) > MAFProcessor._text_col:
                    if tokens[ MAFProcessor._lineType_col] == 's':
                        # Verify if current sequence species match with reference species
                        spec_chrom = SequenceUtils.getSpeciesAndChrom( tokens[ MAFProcessor._speciesChrom_col])
                        species = spec_chrom[0]
                        chromosom = spec_chrom[1]
                        if species == self.referenceSpecies:
                            # Search for BED Sequences having the same <species>.<chromosom>
                            bed_sequences = self.getAssociatedBEDSequences( species + "." + chromosom)
                            if bed_sequences != None and len( bed_sequences) > 0:
                                strand = tokens[ MAFProcessor._strand_col]
                                bp_start = self.computeStartIndex( tokens, strand)
                                text_length = self.getIntValue( tokens[ MAFProcessor._textlength_col])
                                # Search for BEDSequences intersection the current sequence
                                new_block = self.findMatchingBEDSequences( bed_sequences, bp_start, text_length, strand)
                                if new_block != None:
                                    text = tokens[ MAFProcessor._text_col]
                                    new_block.addSequence( species, chromosom,  bp_start, text_length, text)
                                    break
                                else:
                                    # This block does not intersect any BEDSequence. If indexation is used, something is wrong
                                    if indexed == True:
                                        return False
                                    else:
                                        return True
                            else:
                               # This block does not match the chromosom of any BED Sequences
                               # Parsing must be stopped if the file contain only information of one chromosom
                               # Alert is raised in case of index file is used
                               if indexed == True:
                                   Log.log( "MAFProcessor.parseBlock : No BED sequences corresponds to this MSA Block")
                               return False
                        else:
                            Log.log( "MAFProcessor.parseBlock : The first sequence of the parsed block does not correspond to the reference species : " + line)
                            return False
            else:
                # This block is void but we have to continue the parsing
                return True
        
        # If the block matches with at least one of the BED sequence, 
        # parses the rest of the block and store the information
        if new_block != None:
            while 1:
                line = input_file.readline()
                if len( line) == 0:
                    break
                elif not line.isspace():
                    tokens = line.split() 
                    if tokens != None and len( tokens) > MAFProcessor._text_col:
                        if tokens[ MAFProcessor._lineType_col] == 's':
                            spec_chrom = SequenceUtils.getSpeciesAndChrom( tokens[ MAFProcessor._speciesChrom_col])
                            species = spec_chrom[0]
                            chromosom = spec_chrom[1]
                            if len( self.desiredSpeciesList) == 0 or ( len( self.desiredSpeciesList) > 0 and species in self.desiredSpeciesList):
                                bp_start = self.getIntValue( tokens[ MAFProcessor._startindex_col])
                                text_length = self.getIntValue( tokens[ MAFProcessor._textlength_col])
                                text = tokens[ MAFProcessor._text_col]
                                new_block.addSequence( species, chromosom, bp_start, text_length, text)
                            if not species in self.parsedSpeciesList:
                                self.parsedSpeciesList.append( species)
                #Block ends at the first empty line
                else:
                    break
        else:
            Log.log( "MAFProcessor.parseBlock : The parsed block does not contains any sequence (line starting with 's')")
            return False
            
        return True


    # --------------------------------------------------------------------------------------
    # Return the list of BEDsequences having the same species and chromosom
    # than the given ones
    def getAssociatedBEDSequences(self, spec_chrom):
        
        if( self.bedSequencesDict.has_key( spec_chrom)):
            return self.bedSequencesDict[spec_chrom]
        else:
            return None


    # --------------------------------------------------------------------------------------
    # Detect if the given sequence intersect with some of the given BED sequence
    # if so, a new instance of MAFblock is created and attached to the suitable BED Sequence
    def findMatchingBEDSequences(self, bed_sequences, dna_start, text_length, strand):
        
        new_block = None

        dna_end = dna_start + text_length
        for bed_sequence in bed_sequences:
            if dna_end > bed_sequence.indexStart and dna_start < bed_sequence.indexEnd:
                
                if new_block == None:
                    new_block = MAFBlock()
                    if strand == Constants.NEGATIVE_STRAND:
                        new_block.inverseSequences = True
                
                if self.threadLock != None:
                    self.threadLock.acquire()
                if not self.mafBlockDic.has_key( bed_sequence):
                    self.mafBlockDic[ bed_sequence] = []
                    
                self.mafBlockDic[ bed_sequence].append( new_block)
                if self.threadLock != None:
                    self.threadLock.release()
        
        return new_block


    # --------------------------------------------------------------------------------------
    # Compute the start index of the sequence according to the strand
    # if the strand is "-", the coordinates must be inversed
    def computeStartIndex(self, tokens,  strand):
        
        # If current sequence direction is inserved, coordinates must be transformed
        if strand == Constants.NEGATIVE_STRAND:
            source_size = self.getIntValue( tokens[ MAFProcessor._source_size_col])
            text_length = self.getIntValue( tokens[ MAFProcessor._textlength_col])
            rev_start = self.getIntValue( tokens[ MAFProcessor._startindex_col])
            bp_start = source_size + 1 - (text_length + rev_start)
        else:
            bp_start = self.getIntValue( tokens[ MAFProcessor._startindex_col])

        return bp_start


    # --------------------------------------------------------------------------------------
    # Close the given File
    def closeFile(self, file):
        
        try:
            file.close()
        except IOError,  exce:
            Log.log( "MAFProcessor.closeFile : Enable to close file '" + file + "'. From:\n\t--> " +  str(exce))


    # --------------------------------------------------------------------------------------
    # Return the int value of the given token
    def getIntValue(self, token):
        
        try:
            return int( token)
        except ValueError,  val_exce:
            raise ParsingException( "MAFProcessor : Unable to get integer value of '" + token + "'. From:\n\t---> " + str( val_exce))



    # --------------------------------------------------------------------------------------
    # Export the MSA lengths histogram and graph it
    def outputMSALenghtHistogram( self, msa_lenghts, output_commstruct):
        
        # Prepare the processor output dir
        out_path = os.path.join( self.component.outputDir, self.component.getComponentPrefix())
        shutil.rmtree( out_path, True)
        FileUtils.createDirectory( out_path, 0777)
        
        # Output the histogram
        file_infos = RSATUtils.outputHistogram( msa_lenghts, 10, out_path, "MSASize", self.component.pipelineName, "", "Conserved region size", "Number of regions", ('5', '6'), False)

        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.MSA_SIZE_PATH] = file_infos[0]
        output_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.MSA_SIZE_GRAPH_PATH] = file_infos[1]
    




    # --------------------------------------------------------------------------------------
    # Test the result of the MSA composition using an other composition method and comparing the results
    def testFinalMSA(self, final_seq_align, maf_blocks, bed_sequence):
        
        # Note: all the indexes in this first part are coordinates in the genome
        
        # initialize the list that represents the succession of sequence from the MAF blocks
        long_seq = []
        
        # # if necessary, add dots at the beginning of the long_seq representing missing information
        seq_start = maf_blocks[0].sequences[ self.referenceSpecies].indexStart
        if seq_start > bed_sequence.indexStart:
            long_seq.extend( [ Constants.SEQUENCE_INIT_CHAR]*( seq_start -  bed_sequence.indexStart))
        previous_end = seq_start
 
        # Compose the long_seq by the succession of sequence from MAF blocks 
        for maf_block in maf_blocks:
            current_start = maf_block.sequences[ self.referenceSpecies].indexStart
            # inserts dots between two sequences if they are not successive 
            long_seq.extend( ['.']*( current_start - previous_end))
            # create a list from the MAF sequence text
            text_list = list( maf_block.sequences[ self.referenceSpecies].text)
            # localize the insertion characters at the beginning of the list 
            begin = 0
            for i in range( len( text_list)):
                if text_list[i] == Constants.SEQUENCE_INSERTION_CHAR:
                    begin += 1
                else:
                    break
            # localize the insertion characters at the end of the list 
            end = len( text_list)
            for i in range( len( text_list) -1):
                if text_list[-i-1] == Constants.SEQUENCE_INSERTION_CHAR:
                    end = end - 1
                else:
                    break
            # insert the MAF sequence in the long_seq ignoring insertion character at the beginning and at the end
            long_seq.extend( text_list[ begin:end])
            previous_end = current_start + maf_block.sequences[ self.referenceSpecies].textLength

        # if necessary, add dots at the end of the long_seq representing missing information
        if previous_end < bed_sequence.indexEnd:
            long_seq.extend( [ Constants.SEQUENCE_INIT_CHAR]*( bed_sequence.indexEnd -  previous_end))
            
        # compute the index at which the BED sequence may start in the long_seq index coordinates
        if bed_sequence.indexStart < maf_blocks[0].sequences[ self.referenceSpecies].indexStart:
            index_seq_start = 0
        else:
            index_seq_start = bed_sequence.indexStart - maf_blocks[0].sequences[ self.referenceSpecies].indexStart
        
        # compute the index at which the BED sequence may end in the long_seq index coordinates
        index_seq_end = index_seq_start + bed_sequence.indexEnd - bed_sequence.indexStart
        
        # Note : in this second part, we have to consider that insertion characters exists in the sequence text
        # to compute the true star and end index of the BED sequence in the long_seq
        
        # modify the end index according to the number of insertion characters
        indice = 0
        index = 0
        count = 0
        while indice <= index_seq_start:
            if long_seq[index] != Constants.SEQUENCE_INSERTION_CHAR:
                indice += 1
            else:
                count += 1
            index +=1
            
        seq_index_start = index_seq_start + count
        
        # modify the end index according to the number of insertion characters
        indice = 0
        index = 0
        count = 0
        while indice < index_seq_end:
            if long_seq[index] != Constants.SEQUENCE_INSERTION_CHAR:
                indice += 1
            else:
                count += 1
            index +=1
            
        seq_index_end = index_seq_end + count
        
        # retrieve the sub-string of logn_seq that should represent the BED sequence
        result = long_seq[ seq_index_start:seq_index_end]
        
        # compare the string obtained above with the one if the MSA composed in the previous method
        str_result = "".join(result)
        str_final = "".join( final_seq_align.sequences[ self.referenceSpecies])
        
        # If the result are not equals, some thing is wrong
        if str_result != str_final:
            Log.log( "MAFProcessor.testFinalMSA: an error has been detected on the recomposed sequence")
            Log.log( "Composed MSA sequence = " + str_final)
            Log.log( "Test MSA sequence = " + str_result)
            str_long = "".join( long_seq)
            # try to find the MSA sequence in the long_seq directly
            index_test = str_long.find( str_final)
            # if the MSA sequence is not found, we are facing a true issue
            if index_test < 0:
                Log.log( "  The error is confirmed since the composed MSA sequence does not appear in the string composed by the succession of sequences from MAF file : ")                
                Log.log( "  Succession sequence = " + str_long)
                Log.log( "  Bed seq start = " + str( bed_sequence.indexStart))
                Log.log( "  Bed seq end = " + str( bed_sequence.indexEnd))
                Log.log( "  Associated MAF blocks : ")
                for maf_block in maf_blocks:
                   Log.log( maf_block.toString())
            # if the MSA sequence is found, the error come from an index computation issue
            else:
                Log.log( "The error is not confirmed since the composed MSA sequence appears at index " + str( index_test) + " in the string composed by the succession of sequences from MAF file:")
                Log.log( str_long)
                #print "first maf start = " + str( maf_blocks[0].sequences[ self.referenceSpecies].indexStart)
                #print "last maf end = " + str( previous_end)
                #print "bed seq start = " + str( bed_sequence.indexStart)
                #print "bed seq end = " + str( bed_sequence.indexEnd)
                #for maf_block in maf_blocks:
                #   print maf_block.toString()





# ############################################################

#               MAFBlock Class

# ############################################################

from common.DNASequence import DNASequence

class MAFBlock:
    
    # --------------------------------------------------------------------------------------
    def __init__(self):
        self.sequences = {}
        self.inverseSequences = False


    # --------------------------------------------------------------------------------------
    def addSequence(self, species,  chromosom, bp_start, text_lenght,  text):
        
        if self.inverseSequences:
            text = self.inverseText( text)
        sequence = DNASequence( species, chromosom, bp_start, text_lenght, text)
        self.sequences[ species] = sequence


    # --------------------------------------------------------------------------------------
    def inverseText(self,  text):
        
        result = ""
        length = len( text)
        for index in range( length):
            result += text[ length - index -1]
            
        return result
    
    
    # --------------------------------------------------------------------------------------
    def toString(self,  reference_species=""):
        
        result = "MAF Block : \n"
        if reference_species == "":
            for species in sorted( self.sequences.keys()):
                sequence = self.sequences[ species]
                result += "|--" + sequence.toString() + "\n"
        else:
            sequence = self.sequences[ reference_species]
            result += "|--" + sequence.toString() + "\n\t\t" + sequence.text + "\n"
            for species in sorted( self.sequences.keys()):
                if species != reference_species:
                    sequence = self.sequences[ species]
                    result += "|--" + sequence.toString() + "\n"
        return result





# eflag: FileType = Python2
