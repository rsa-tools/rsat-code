
from utils.Constants import Constants

class Motif:
    
    # --------------------------------------------------------------------------------------
    def __init__( self, start, end, name, pwm):
        
        self.indexStart = start
        self.indexEnd = end
        self.name = name
        self.id = ""
        self.pwm = pwm
        self.consensus = ""
        self.strand = Constants.POSITIVE_STRAND
        self.offset = 0
        self.score = 0
	self.uncount = 0


    # --------------------------------------------------------------------------------------
    # Create the name on the basis of provided information
    def composeName(self, base_name):

       self.name = base_name + "_" + str( self.indexStart) + "." + str( self.indexEnd)
    
    
    # --------------------------------------------------------------------------------------
    # Compares the two given motifs
    @staticmethod
    def compare( motif1, motif2):
        
        result = 0
        if motif1.indexStart != motif2.indexStart:
            result = motif1.indexStart - motif2.indexStart
        else:
            if motif1.name < motif2.name:
                result = -1
            elif motif1.name > motif2.name:
                result = 1

        return result
    
    
    # --------------------------------------------------------------------------------------
    # Return a string representation of the motif
    def toString(self):
        
        result = "Motif :" + " name=" + self.name + " start=" + str( self.indexStart) + " end=" + str( self.indexEnd) + " strand=" + self.strand
        if self.consensus != None and len( self.consensus) > 0:
            result += " consensus=" + self.consensus
        if self.pwm != None:
           result += "\nPWM=" + self.pwm.toString( False)

        return result
