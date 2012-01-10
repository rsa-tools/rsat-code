
from utils.log.Log import Log
from utils.exception.ExecutionException import ExecutionException
from utils.Constants import Constants

class SequenceAlignment:


    # --------------------------------------------------------------------------------------
    def __init__( self):
        
        self.sequences = {}
        self.motifs = []
        self.totalLength = 0
        self.name = ""
        self.referenceSpecies = ""


    # --------------------------------------------------------------------------------------
    # Insert in the sequences dictionnary an entry for each species in the given list
    # with key = species and value = a list of '-' of length equals to given one
    def initializeWithDots(self, length, reference_species, species_list):
        
        self.referenceSpecies = reference_species
        
        if length > 0 and len( species_list) > 0:
            self.sequences[reference_species] = []
            for i in range( length):
                self.sequences[reference_species].append( Constants.SEQUENCE_INIT_CHAR)
            for species in species_list:
                self.sequences[species] = []
                for i in range( length):
                     self.sequences[species].append( Constants.SEQUENCE_INIT_CHAR)


    # --------------------------------------------------------------------------------------
    # Insert a block of text in the sequence corresponding to the given species
    def insertSequenceBlock(self, species, index_start, index_end, text):
        
        if self.sequences.has_key( species):
            self.sequences[ species][ index_start : index_end] = text


    # --------------------------------------------------------------------------------------
    # For all sequences of the alignment, replace the character used 
    # to initialize the sequences by the bp insertion character "-"
    # and remove the gaps in the sequences
    def finalizeSequences(self, keep_gaps = False):
        
        # Compute the length of the longest sequence in the MSA
        for species in self.sequences.keys():
            sequence = self.sequences[species]
            if len( sequence) > self.totalLength:
                self.totalLength = len( sequence)

        # Add insertion characters at the end of sequence that does not have the right length
        for species in self.sequences.keys():
            sequence = self.sequences[species]
            if len( sequence) != self.totalLength:
                Log.info( "SequenceAlignment.finalizeSequences : Sequence does not have the right lenght for this alignment : Alignement length = " + str( self.totalLength) + " DNA sequence length = " + str( len( sequence)) + " for species= " + species + ". Completing sequence")
                for fix_index in range( self.totalLength - len( sequence)):
                    sequence.append( Constants.SEQUENCE_INSERTION_CHAR)

        # Analyse each position of the MSA and decide if the column should be kept or removed
        species_list = self.sequences.keys()
        removed_char = 0
        for index in range( self.totalLength):
            all_insertion = True
            # if the residus on the reference species sequence indicates a missing information, it means that
            # the MAF files provide no information for this column. So the column is kept
            char = self.sequences[ self.referenceSpecies][ index - removed_char]
            if char == Constants.SEQUENCE_INIT_CHAR:
                continue
                
            # if keep_gaps is False, remove the column if the residu on the reference species sequence
            # is an insertion character
            if keep_gaps == False:
                if char == Constants.SEQUENCE_INSERTION_CHAR:
                    for index_sub_species in range( len( species_list)):
                        species = species_list[index_sub_species]
                        self.sequences[ species].pop( index - removed_char)
                    removed_char += 1
                continue
                
            # remove the columns that contains only insertion characters or initialization characters
            for index_species in range( len(species_list)):
                species = species_list[index_species]
                char = self.sequences[species][index - removed_char]
                if char != Constants.SEQUENCE_INSERTION_CHAR and char != Constants.SEQUENCE_INIT_CHAR:
                    all_insertion = False
                    break
            if all_insertion == True:
                for index_sub_species in range( len(species_list)):
                    species = species_list[index_sub_species]
                    self.sequences[species].pop( index - removed_char)
                removed_char += 1
            
        self.totalLength = self.totalLength - removed_char     


    # --------------------------------------------------------------------------------------
    # Insert directly a sequence in the sequences dictionnary
    def addSequence( self, species, sequence):
        
        if species != None and sequence != None:
            self.sequences[ species] = sequence
            seq_length = len( sequence)
            if self.totalLength == 0:
                self.totalLength = seq_length
            else:
                if seq_length != self.totalLength:
                    Log.log( "SequenceAlignment.addSequence : Added sequence does not have the right lenght for this alignment : Alignement length = " + str( self.totalLength) + " DNA sequence length = " + str( seq_length))
                    for fix_index in range( self.totalLength - len( sequence)):
                        sequence.append( Constants.SEQUENCE_INSERTION_CHAR)
                
    
    
    # --------------------------------------------------------------------------------------
    # Add a motif to the sequence motif list. If mergeContiguous is True, before inserting the new motif,
    # the list is parsed to search motifs contiguous to the new one. If it exist contiguous motif, they are
    # merged.
    def addMotif( self, new_motif, mergeContiguous = False):
        
        if new_motif != None:
            if mergeContiguous and len( self.motifs) > 0:
                merged = False
                for current_motif in self.motifs:
                    if current_motif.indexStart == new_motif.indexEnd or current_motif.indexEnd == new_motif.indexStart:
                        merged = True
                        new_motif = self.mergeMotifs( current_motif, new_motif)
                if not merged:
                    self.motifs.append( new_motif)
            else:
                self.motifs.append( new_motif)
    
    
    # --------------------------------------------------------------------------------------
    # Merge the second given motif into the first one
    def mergeMotifs(self, master_motif, added_motif):
        
        if master_motif.indexStart == added_motif.indexEnd:
            master_motif.indexStart = added_motif.indexStart
            after = False
        elif master_motif.indexEnd == added_motif.indexStart:
            master_motif.indexEnd = added_motif.indexEnd
            after = True
        else:
            raise ExecutionException( "SequenceAlignement.mergeMotifs : The two given motifs are not contiguous")
        
        master_motif.composeName( self.name)
        master_motif.pwm.mergeMatrix( added_motif.pwm, after)
    


    # --------------------------------------------------------------------------------------
    # Write the Sequence Alignment in FASTA format
    def convertToFASTA(self, species_list = None):
        
        result = ""
        for species in self.sequences.keys():
            if species_list == None or len( species_list) == 0 or species in species_list: 
                result +="> " + species + "\n"
                for residu in self.sequences[species]:
                    result += residu
                result += "\n"
        
        return result
    
    
    
    # --------------------------------------------------------------------------------------
    # modify the value of the given index considering this index is relative to the text of
    # the reference species sequence (including the insertions) and should be modified to corresponding
    # to an index relative to the genome (without the insertions)
    def fixIndex(self, text_index):
        
        if text_index < 0:
            return text_index
        
        if self.referenceSpecies in self.sequences.keys():
            limit = min( text_index + 1, len( self.sequences[ self.referenceSpecies]))
            count = 0
            for index in range( limit):
                if self.sequences[ self.referenceSpecies][index] == Constants.SEQUENCE_INSERTION_CHAR:
                    count += 1
            return text_index - count
        else:
            Log.log( "SequenceAlignement.fixIndex : Reference species is not set for Sequence Alignement : " + self.name)
            return text_index
    
    
    
    # --------------------------------------------------------------------------------------
    # Return a tring representation of the SequenceAlignement
    def toString(self):
        
        result = "|Sequence Alignement : name=" + self.name + "\n"
        result += "|--Sequences\n"
        for species in self.sequences.keys():
            result += "|--|-" + species + (" "* (10-len( species))) + "".join( self.sequences[ species]) + "\n"
        result += "|--Motifs\n"
        for motif in self.motifs:
            result += "|--|-" + motif.toString() + "\n"
            
        return result

# eflag: FileType = Python2
