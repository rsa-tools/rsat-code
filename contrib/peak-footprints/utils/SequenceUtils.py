
import re

from utils.Constants import Constants

class SequenceUtils:

    # --------------------------------------------------------------------------------------
    # Extract the species name and the chromosom name (if any) from a string token
    # with format '<species>.<chromosom>' or '<species>'
    @staticmethod
    def getSpeciesAndChrom( spec_chrom):
        
        matcher = re.compile('(.*)\.(.*)')
        matches = matcher.search( spec_chrom)
        
        if matches != None and len( matches.groups()) >= 2:
            species_and__chrom = [ matches.group(1),  matches.group(2)]
        else:
            species_and__chrom = [ spec_chrom,  ""]

        return species_and__chrom    
    
    
    
    # --------------------------------------------------------------------------------------
    # Build a new sequence text removing the insertion characters from the original sequence text
    @staticmethod
    def cleanSequenceText( text):
        
        cleaned = ""
        for index in range( len( text)):
            if text[index] != Constants.SEQUENCE_INSERTION_CHAR:
                cleaned += text[index]
        
        return cleaned

# eflag: FileType = Python2
