
import math

from utils.exception.ExecutionException import ExecutionException

from utils.Constants import Constants

# This class corresponds to a PWM (position weight matrix)

class PWM:

    # --------------------------------------------------------------------------------------
    def __init__( self):

        self.nbSequences = 0
        self.totalLength = 0
        self.matrix = None
        self.ratioMatrix = None
        self.informationMatrix = None
        self.correctionValue = self.nbSequences
        self.informationLimits = None


    # --------------------------------------------------------------------------------------
    # Insert a new column of value in the PWM. Values must be provided in residu alphabet order
    def addColumn(self, values):
        
        if len(values) != len( Constants.DNA_ALPHABET):
            raise ExecutionException( "PWM.addColumn : Incorrect number of residu values : " + str( values))
        
        if self.matrix == None:
            self.matrix = {}
            for letter in Constants.DNA_ALPHABET:
                self.matrix[ letter] = []
            self.matrix[ Constants.MAX_INDEX] = []
        
        letter_index = 0
        max = 0
        letter_max = None
        for letter in sorted( Constants.DNA_ALPHABET):
            value =  values[letter_index]
            self.matrix[ letter] .append( value)
            if value >= max:
                max = value
                letter_max = letter
            letter_index += 1
        
        self.matrix[ Constants.MAX_INDEX].append( letter_max)
        
        self.totalLength += 1
        


    # --------------------------------------------------------------------------------------
    # Initialize the PWM using the given alignment, computing the score and ratio matrix
    def initFromAlignment(self, alignment, desired_species):

        self.totalLength = alignment.totalLength
        
        self.initMatrix( alignment, desired_species)
        self.computeRatioMatrix()
        self.computeInformationMatrix()
    
    
    # --------------------------------------------------------------------------------------
    # Compute the number of repetition of a letter at a position through the
    # sequences of the alignement
    def initMatrix( self, alignment, desired_species):
        
        #Initialise the matrix
        self.matrix ={}
        if self.totalLength > 0:
            for letter in Constants.DNA_ALPHABET:
                self.matrix[ letter] = self.totalLength*[0]
            self.matrix[ Constants.MAX_INDEX] = self.totalLength*[ Constants.SEQUENCE_INSERTION_CHAR]
            
        self.nbSequences = 0
        
        # Count the residu occurences for each position
        for species in alignment.sequences.keys():
            if desired_species == None or len( desired_species) == 0 or species in desired_species:
                self.nbSequences += 1
                sequence = alignment.sequences[ species]
                for index in range( self.totalLength):
                    letter = sequence[index]
                    if letter in Constants.DNA_ALPHABET:
                        self.matrix[ letter][ index] += 1
        
        # Determine the residu having the max value for each position
        for index in range( self.totalLength):
            max = 0
            for letter in Constants.DNA_ALPHABET:
                if self.matrix[ letter][ index] > max:
                    self.matrix[ Constants.MAX_INDEX][ index] = letter

        # TO REMOVE ONCE ALGORITHM OK
        self.correctionValue = self.nbSequences

    # --------------------------------------------------------------------------------------
    # Return the residu with the greater occurence at the given position
    def getMostConservedResidu(self, index):
        
        letter_max = self.matrix[ Constants.MAX_INDEX][ index]
        if letter_max in Constants.DNA_ALPHABET:
            return letter_max
        else:
            return None


    # --------------------------------------------------------------------------------------
    # Compute the conservation ratio of each letter at each position
    def computeRatioMatrix( self):
        
        self.ratioMatrix ={}
        for letter in Constants.DNA_ALPHABET:
            self.ratioMatrix[ letter] = self.totalLength*[0]
        self.ratioMatrix[ Constants.MAX_INDEX] = self.totalLength*[ Constants.SEQUENCE_INSERTION_CHAR]
        
        for index in range( self.totalLength):
            ratio_max = 0.0
            letter_max = ""
            for letter in Constants.DNA_ALPHABET:
                ratio = self.matrix[ letter][ index] / float( self.nbSequences)
                self.ratioMatrix[ letter][ index] = ratio
                if ratio > ratio_max:
                    ratio_max = ratio
                    letter_max = letter
            self.ratioMatrix[ Constants.MAX_INDEX][ index] = letter_max
   

    # --------------------------------------------------------------------------------------
    # Compute the information value of each letter at each position
    def computeInformationMatrix( self):
        
        self.informationMatrix ={}
        for letter in Constants.DNA_ALPHABET:
            self.informationMatrix[ letter] = self.totalLength*[0]
        self.informationMatrix[ Constants.MAX_INDEX] = self.totalLength*[ Constants.SEQUENCE_INSERTION_CHAR]
        
        for index in range( self.totalLength):
            info_max = -1000000
            letter_max = ""
            for letter in Constants.DNA_ALPHABET:
                probi = float( Constants.HG_BACKGROUND_MODEL[letter])
                fij = (self.matrix[ letter][ index] + self.correctionValue * probi) / float( self.nbSequences + self.correctionValue)
                Iij = fij * math.log( fij / probi)
                self.informationMatrix[ letter][index] = Iij
                if Iij > info_max :
                    info_max = Iij
                    letter_max = letter
            self.informationMatrix[ Constants.MAX_INDEX][ index] = letter_max
        
        self.computeInformationLimits()


    # --------------------------------------------------------------------------------------
    # Compute the lower and upper limits the information value can reach for each letter
    def computeInformationLimits( self):
        
        self.informationLimits ={}

        for letter in Constants.DNA_ALPHABET:
            self.informationLimits[ letter] =2*[0]
            probi = float( Constants.HG_BACKGROUND_MODEL[letter])
            fmin = (self.correctionValue * probi) / float( self.nbSequences + self.correctionValue)
            self.informationLimits[ letter][0] = fmin * math.log( fmin / probi)
            fmax = (self.nbSequences + self.correctionValue * probi) / float( self.nbSequences + self.correctionValue)
            self.informationLimits[ letter][1] = fmax * math.log( fmax / probi)


    # --------------------------------------------------------------------------------------
    # Return the PWM located between the two indexes as a dict key = letter/ value = tuple of counts
    def getPWM(self, index_start, index_end):
        
        result = PWM()
        result.matrix={}
        
        for letter in Constants.DNA_ALPHABET:
            result.matrix[ letter] = ( index_end - index_start)*[0]
        result.matrix[ Constants.MAX_INDEX] = ( index_end - index_start)*[ Constants.SEQUENCE_INSERTION_CHAR]
            
        for index in range( index_start, index_end):
            for letter in Constants.DNA_ALPHABET:
                result.matrix[letter][index - index_start] = ( self.matrix[ letter][ index])
            result.matrix[ Constants.MAX_INDEX][index - index_start] = self.matrix[ Constants.MAX_INDEX][ index]
                    
        result.nbSequences = self.nbSequences
        result.totalLength = index_end - index_start
        result.correctionValue = self.correctionValue
        
        result.computeRatioMatrix()
        result.computeInformationMatrix()
        
        return result


    # --------------------------------------------------------------------------------------
    # Add the given PWM at the end or at the beginning ot the current matrix
    def mergeMatrix( self, added_pwm, after):
        
        if after:
            for letter in Constants.DNA_ALPHABET:
                self.matrix[ letter][ self.totalLength :] = added_pwm.matrix[letter]
            self.matrix[ Constants.MAX_INDEX][ self.totalLength :] = self.matrix[Constants.MAX_INDEX]
        else:
            for letter in Constants.DNA_ALPHABET:
                self.matrix[ letter][0:0] = added_pwm.matrix[letter]
            self.matrix[ Constants.MAX_INDEX][0:0] = added_pwm.matrix[Constants.MAX_INDEX]
        
        self.totalLength += added_pwm.totalLength
        
        self.computeRatioMatrix()
        self.computeInformationMatrix()


    # --------------------------------------------------------------------------------------
    # Convert the PWM as a string in TRANSFAC Format
    def convertToTransfac(self):
        
        result = "PO\tA\tC\tG\tT\n"
        
        for position in range( self.totalLength):
            result += str( position)
            for letter in sorted( Constants.DNA_ALPHABET):
                result += "\t" + str( self.matrix[ letter][position])
            result += "\n"
        result += "XX\n"
        result += "BA\t" + str( self.nbSequences) + "\n"
            
        return result


    # --------------------------------------------------------------------------------------
    # Convert the PWM to a horizontal tabbed format
    def convertToHorizontaltab( self):
        
        result = ""
        for letter in sorted( Constants.DNA_ALPHABET):
            result += letter + "\t"
            for value in self.matrix[letter]:
                result += str( value) + "\t"
            result +="\n"
        
        return result
        

    # --------------------------------------------------------------------------------------
    # Convert the PWM to a vertical tabbed format
    def convertToVerticaltab( self):
        
        result = "A\tC\tG\tT\n"
        for position in range( self.totalLength):
            for letter in sorted( Constants.DNA_ALPHABET):
                result += str( self.matrix[ letter][position]) + "\t"
            result += "\n"
        
        return result    


    # --------------------------------------------------------------------------------------
    # Convert the PWM to MEME minimal format
    def convertToMEME(self, name):
        
        result = ""
        result += "MOTIF " + name + " " + name + "\n\n"
        result += "letter-probability matrix: alength= 4 w= " + str( self.totalLength) + " nsites= " + str( self.nbSequences) + " E= 0\n"
        for position in range( self.totalLength):
            result += "  "+ str( self.matrix['A'][position]/float( self.nbSequences))
            result += "\t  " + str( self.matrix['C'][position]/float( self.nbSequences))
            result += "\t  " + str( self.matrix['G'][position]/float( self.nbSequences))
            result += "\t  " + str( self.matrix['T'][position]/float( self.nbSequences))
            result += "\n"
            
        return result
    
    
    # --------------------------------------------------------------------------------------
    # Compare the given PWM to the instance. Return True if PWM are the same
    def equals(self, other_pwm):
        
        try:
            for letter in self.matrix.keys():
                for index in range( len( self.matrix[letter])):
                    if self.matrix[letter][index] !=  other_pwm.matrix[letter][index]:
                        return False
        except (KeyError, IndexError):
            return False
            
        return True
            

    # --------------------------------------------------------------------------------------
    # Return the PWM variables as a string
    def toString( self, all_info = True):
        
        result = "\n"
        
        if all_info:
            result += " Frequencies: \n"
        
        #result += self.convertToMEME( "test")
        
        result += "\n"
        
        result += self.convertToHorizontaltab()
        
        result += "\n"
        
        if all_info:
            result += "\n Ratios: \n"
            for letter in sorted( Constants.DNA_ALPHABET):
                result += letter + "\t"
                for value in self.ratioMatrix[ letter]:
                    result += str( value) + "\t"
                result +="\n"
                
            result += "\n Max ratios: \n"
            for value in self.ratioMatrix[ Constants.MAX_INDEX]:
                    result += str( value) + "\t"
                    
        return result
    
# eflag: FileType = Python2
