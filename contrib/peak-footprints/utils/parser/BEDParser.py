
from utils.log.Log import Log
from utils.exception.ParsingException import ParsingException

from common.BEDSequence import BEDSequence

class BEDParser:
    
    _chrom_col = 0
    _startindex_col = 1
    _endindex_col = 2
    _id_col = 3
    _strand_col = 5


    # --------------------------------------------------------------------------------------
    # Parse the given BED file and return a dictionnary of the BED Sequences
    # grouped by sequence keys ('species'.'chromosom')
    @staticmethod
    def getBEDSequenceDictionnary( species, bed_filepath, extension_5p, extension_3p):
        
        sequence_dic = {}
        
        try:
            input_file = open( bed_filepath)
            for line in input_file:
                tokens = line.split()
                if len( tokens) > BEDParser._endindex_col:
                    chrom = tokens[ BEDParser._chrom_col].lower()
                    #if chrom[ 0:3] == "chr":
		    if chrom[0:1] != "#" and chrom[ 0:3] == "chr":
			if len(chrom) < 4:
			    chrom = "chr" + chrom
                        start = BEDParser.getTokenAsint( tokens[ BEDParser._startindex_col])
                        end = BEDParser.getTokenAsint( tokens[ BEDParser._endindex_col])
                        if start < end:
                            start = start - extension_5p;
                            if start < 0:
                                start = 0;
                            end = end + extension_3p;
                            bedsequence = BEDSequence( species, chrom, start, end )
                            if len( tokens) > BEDParser._id_col:
                                bedsequence.id = tokens[ BEDParser._id_col]
                            bedsequence_key = bedsequence.getKey()
                            if not sequence_dic.has_key( bedsequence_key):
                                sequence_dic[ bedsequence_key] = []
                            sequence_dic[ bedsequence_key].append( bedsequence)
                        else:
                            Log.log( "BEDParser.getBEDSequenceDictionnary : A sequence has inversed start and end coordinates : " + line)
		else:
		    Log.log( "No 'chr' in line :" + line)
        except ParsingException, par_exce:
            raise ParsingException( "BEDParser.getBEDSequenceDictionnary : Some attributes are mor numbers. From:\n\t-->  " +  str( par_exce))
        except IOError, io_exce:
            raise ParsingException( "BEDParser.getBEDSequenceDictionnary : Unable to open the file '" + bed_filepath +  "'. From:\n\t-->  " +  str( io_exce))
            
        return sequence_dic

    
    # --------------------------------------------------------------------------------------
   # Return the value of the token as int if possible and raise an exception if the value is required but not converted
    @staticmethod
    def getTokenAsint( token, required=True):
        
        try:
            att_value = int( token)
            return att_value
        except (TypeError, ValueError), val_exce:
            if required:
                raise ParsingException( "BEDParser.getTokenAsint : Unable to convert the token to int :'" + token + "'. From:\n\t---> " + str( val_exce))
            else:
                return None



# eflag: FileType = Python2
