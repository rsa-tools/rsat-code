
class DNASequence:


    # --------------------------------------------------------------------------------------
    def __init__( self, species, chrom, bp_start, text_length, text):
        
        self.species = species
        self.chromosom = chrom
        self.indexStart = bp_start
        self.textLength = text_length
        self.text = text.upper()
    

    # --------------------------------------------------------------------------------------
    def getKey( self):
        
        return self.species + "." + self.chromosom


    # --------------------------------------------------------------------------------------
    def toString( self):
        
        result = "DNASequence : " + self.getKey() + " start=" + str( self.indexStart) + " text length=" + str( self.textLength)
        result += "\n" + self.text
        return result
