
class BEDSequence:
    
    # --------------------------------------------------------------------------------------
    def __init__( self, species, chrom, bp_start, bp_end):
        
        self.species = species
        self.chromosom = chrom
        self.indexStart = bp_start
        self.indexEnd = bp_end
        self.name = self.getKey() + "_" + str( self.indexStart) + "." + str( self.indexEnd)
        self.id = ""
        self.referenceIndex =  self.indexStart + int( (self.indexEnd - self.indexStart) / float(2))
        self.score = 0.0


    # --------------------------------------------------------------------------------------
    def getKey( self):
        
        return self.species + "." + self.chromosom



    # --------------------------------------------------------------------------------------
    def getLength( self):
    
        return self.indexEnd - self.indexStart

    # --------------------------------------------------------------------------------------
    @staticmethod
    def compare( bedseq1, bedseq2):
        
        if bedseq1.chromosom == bedseq2.chromosom:
            return bedseq1.indexStart - bedseq2.indexStart
        elif bedseq1.chromosom < bedseq2.chromosom:
            return -1
        elif bedseq1.chromosom > bedseq2.chromosom:
            return 1

    # --------------------------------------------------------------------------------------
    def toString( self):
        
        result = "BEDSequence : name=" + self.name + " start=" + str( self.indexStart) + " end=" + str( self.indexEnd)
        
        return result

