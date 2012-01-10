
import os, random, shutil, commands, math

from processor.Processor import Processor
from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct

from common.Motif import Motif

from manager.ProgressionManager import ProgressionManager

from utils.Constants import Constants
from utils.RSATUtils import RSATUtils
from utils.log.Log import Log
from utils.MotifUtils import MotifUtils
from utils.exception.ExecutionException import ExecutionException

class ImplantSitesProcessor( Processor):
    
    
    MOTIF_LIST_PARAM = "MotifList"
    OPTIMIZE_MOTIF_PARAM = "OptimizeMotif"
    DATABASE_FILE_PATH_PARAM = "DatabaseFilePath"
    SITE_NUMBER_PARAM = "SiteNumber"
    DISTRIBUTION_MODE_PARAM = "DistributionMode"
    UNIFORM_DISTRIBUTION_MODE_VALUE = "uniform"
    CENTERED_DISTRIBUTION_MODE_VALUE = "normal"
    
    
    # ---------------------------------------------------------------------------------------------
    def __init__( self):
        Processor.__init__( self)


    # ---------------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as input 
    # (None if no input CommStruct)
    @staticmethod
    def getInputCommStructClass():
        
        return ( BedSeqAlignmentStatsCommStruct, )


    # ---------------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as output
    # (None if no output CommStruct)
    @staticmethod
    def getOutputCommStructClass():
        
        return ( BedSeqAlignmentStatsCommStruct, )
        
        

    # --------------------------------------------------------------------------------------
    # Returns a name that will be used as display name in the user friendly outputs
    @staticmethod
    def getDisplayName():
        
        return "Implantation of motif sites in MSA"
        



    #------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return ( ImplantSitesProcessor.SITE_NUMBER_PARAM, ImplantSitesProcessor.MOTIF_LIST_PARAM, ImplantSitesProcessor.OPTIMIZE_MOTIF_PARAM, ImplantSitesProcessor.DATABASE_FILE_PATH_PARAM, ImplantSitesProcessor.DISTRIBUTION_MODE_PARAM)



    
    # ---------------------------------------------------------------------------------------------
    # Execute the processor
    def execute( self, input_commstructs):

        if input_commstructs == None or len( input_commstructs) == 0:
            raise ExecutionException( "ImplantSitesProcessor.execute : No inputs")
        
        input_commstruct = input_commstructs[0]

        # Implant TF Motif binding sites in mSA Sequences
        Log.trace( "ImplantSitesProcessor.execute : Implanting motif sites")
        ProgressionManager.setTaskProgression( "Implanting motif sites", self.component, 0.0)
        self.addSites( input_commstruct)
        ProgressionManager.setTaskProgression( "Implanting motif sites", self.component, 1.0)
        
        return input_commstruct



    # ---------------------------------------------------------------------------------------------
    # Add biding site randomly distributed in the sequences of MSA
    def addSites(self, output_commstruct):
        
        # Retrieve the algorithm parameters
        site_number = self.getParameterAsint( ImplantSitesProcessor.SITE_NUMBER_PARAM)

        if site_number <= 0:
            Log.trace( "ImplantSitesProcessor.addSites : Motif sites implantation not requested")
            return
            
        motif_list_line = self.getParameter( ImplantSitesProcessor.MOTIF_LIST_PARAM)
        motif_name_list = motif_list_line.split()
        
        optimize_motif = (self.getParameter( ImplantSitesProcessor.OPTIMIZE_MOTIF_PARAM).lower() == "true")
        
        database_file_path = self.getParameter( ImplantSitesProcessor.DATABASE_FILE_PATH_PARAM)
        
        distribution_mode = self.getParameter( ImplantSitesProcessor.DISTRIBUTION_MODE_PARAM)
        distribution_mode = distribution_mode.lower()

        # Retrieve the motifs PWM
        motif_def_list = self.getMotifDefinitions( motif_name_list, database_file_path)

        # Prepare output directory
        dir_path = os.path.join( self.component.outputDir, self.component.getComponentPrefix())
        shutil.rmtree( dir_path, True)
        os.mkdir( dir_path)


        # Generate the motif sites
        motif_sites = {}
        for motif in motif_def_list:
            if optimize_motif == False:
                motif_file_path = self.outputMotifDefinition( motif, dir_path)
                motif_sites[motif] = self.generateRandomSites( motif, motif_file_path, site_number)
            else:
                motif_sites[motif] = self.generateOptimalSites( motif, site_number)
                
        # Implant sites in the MSA
        self.implantSites( motif_sites, distribution_mode, output_commstruct, dir_path)



    
    # ---------------------------------------------------------------------------------------------
    # Retrieve the definition of given motifs in given database file
    def getMotifDefinitions(self, motif_name_list, database_file_path):
        
        motif_list = []
        for name in motif_name_list:
            motif_list.append( Motif( 0, 0, name, None))
            
        MotifUtils.getMotifsPWMFromJasparTF( motif_list, database_file_path)
                
        return motif_list



    # ---------------------------------------------------------------------------------------------
    # Output the transfac motif definition to file
    def outputMotifDefinition(self, motif, dir_path):
        
        file_path = os.path.join( dir_path, motif + ".tab")
        
        definition =  motif.pwm.convertToHorizontaltab()
        
        try:
            motif_file = open( file_path, "w")
            motif_file.write( definition)
            motif_file.flush()
            motif_file.close()
        except IOError, io_exce:
            raise ExecutionException( "ImplantSitesProcessor.outputMotifDefinition : Unable to write motif definition to file '" + file_path + "'. From:\n\t---> " + str( io_exce))

        return file_path



    # ---------------------------------------------------------------------------------------------
    # Output the transfac motif definition to file
    def generateRandomSites( self, motif, motif_file_path, site_number):
        
        # Retrieve method required parameters
        RSAT_PATH = self.component.getParameter( Constants.RSAT_DIR_PARAM)
        dir_path = os.path.join( self.component.outputDir, self.component.getComponentPrefix())
        output_path = os.path.join( dir_path, motif + "_sites.fasta")
        
        # Execute the RSAT random-seq command
        cmd = os.path.join( RSAT_PATH , "python-scripts/random-sites")
        cmd += " -m " + motif_file_path
        cmd += " -n " + str( site_number)
        cmd += " -o " + output_path
        
        # Execute the command
        cmd_result = commands.getstatusoutput( cmd)
        if cmd_result[0] != 0:
            Log.log( "ImplantSitesProcessor.generateSites : status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
            Log.log( "ImplantSitesProcessor.generateSites : command output is = \n" + str( cmd_result[1]))
            raise ExecutionException( "ImplantSitesProcessor.generateSites : Cannot execute random-sites commands. See logs for more details")
        
        # Parse the result of the command
        sites = []
        try:
            site_file = open( output_path, "r")
            for line in site_file:
                if not line.isspace() and line[0] != ">":
                    sites.append( line.split()[0].upper())
            site_file.close()
        except IOError, io_exce:
            raise ExecutionException( "ImplantSitesProcessor.generateSites : Unable to read motif sites from file '" + output_path + "'. From:\n\t---> " + str( io_exce))

        return sites


    # ---------------------------------------------------------------------------------------------
    # Generate the site corresponding to the maximum occurence residu at each position
    def generateOptimalSites( self, motif, site_number):
        
        site = ""
        for index in range( motif.pwm.totalLength):
            site += motif.pwm.matrix[ Constants.MAX_INDEX][index]
        
        site = site.upper()
        
        return [site] * site_number


    # ---------------------------------------------------------------------------------------------
    # Implant the sites in the MSA in position chosen according the chosen distribution mode
    def implantSites(self, site_list, distribution_mode, output_commstruct, dir_path):
        
        # collect the binding places for motif sites
        implantations = {}
        for motif in site_list.keys():
            self.chooseBindingPoints( motif, site_list[ motif], distribution_mode, implantations, output_commstruct, dir_path)
        
        # implants the motif site to the chosen biding points
        for bedseq in implantations.keys():
            for place in implantations[ bedseq]:
                for msa in output_commstruct.bedToMA[bedseq]:
                    for sequence in msa.sequences.values():
                        index_start = place[1] - bedseq.indexStart
                        index_end = place[2] - bedseq.indexStart
                        sequence[ index_start:index_end] = place[0]



    # ---------------------------------------------------------------------------------------------
    # Randomly choose a position to place sites through the BEDSequences according to the chosen distribution mode
    def chooseBindingPoints(self, motif, sites, distribution_mode, implantations, output_commstruct, dir_path):
        
        bedseq_list = output_commstruct.bedToMA.keys()
        bedseq_list_length = len( bedseq_list)
        
        chosen_distances_signed = []
        
        # Case of normal distribution
        # ...........................
        if distribution_mode == ImplantSitesProcessor.CENTERED_DISTRIBUTION_MODE_VALUE:
            for site in sites:
                tries = 0
                while True:
                    # Draw a bedseq with uniform probability
                    chosen_bedseq = bedseq_list[ int( random.uniform(0, bedseq_list_length))]
                    # Choose start index using normal distribution around peak reference index
                    chosen_middle_index = int( random.normalvariate( chosen_bedseq.referenceIndex, 30.0))
                    chosen_start_index = chosen_middle_index - int(len( site) / float( 2))
                    chosen_end_index = chosen_middle_index + int( math.ceil(len( site) / float( 2)))
                    chosen_distances_signed.append( chosen_middle_index - chosen_bedseq.referenceIndex )
                    # Test if any other site previously placed intersect the chosen one
                    intersect = False
                    if chosen_bedseq in implantations.keys():
                        for previous_indexes in implantations[ chosen_bedseq]:
                            if chosen_start_index < previous_indexes[2] and chosen_end_index > previous_indexes[1]:
                                intersect = True
                                break
                    # TO REMOVE
                    intersect = False
                    # END TO REMOVE
                    # If position is free, add the chosen position to the list of implantations
                    if not intersect:
                        if not chosen_bedseq in implantations.keys():
                            implantations[ chosen_bedseq] = []
                        implantations[ chosen_bedseq].append( (site, chosen_start_index, chosen_end_index))
                        break
                    else:
                        tries += 1
                        if tries > 50:
                            Log.trace("GenerateMSAProcessor.chooseBindingPoints : No place found for site : " + site + ". Bypassing site")
                            break
        
        # Case of uniform distribution
        # ............................
        elif distribution_mode == ImplantSitesProcessor.UNIFORM_DISTRIBUTION_MODE_VALUE:
            
            # Build a table that will permit to easily find a bedseq when drawing will be done
            total_length = 0
            bedseq_limits =[]
            for bedseq in bedseq_list:
                length = bedseq.getLength()
                bedseq_limits.append( total_length + length)
                total_length += length
                
            for site in sites:
                tries = 0
                while True:
                    # draw a number (uniform) over all BED Sequence indexes
                    drawen_index = random.randint( 0, total_length-1)
                    # find to which BED sequence and to which index correspond the drawed number
                    chosen_bedseq = None
                    for index in range( bedseq_list_length):
                        if drawen_index < bedseq_limits[ index]:
                            chosen_bedseq = bedseq_list[ index]
                            if index == 0:
                                chosen_middle_index = chosen_bedseq.indexStart + drawen_index + int( len( site) / float(2))
                            else:
                                chosen_middle_index = chosen_bedseq.indexStart + drawen_index - bedseq_limits[ index -1 ]
                            break
                    # if the BEDseq is correctly found, check if the place is good for the site
                    if chosen_bedseq != None:
                        # if the index is to near from the sequence endpoints, draw a new index
                        if chosen_middle_index < (chosen_bedseq.indexStart + int( len( site) / float(2))) or chosen_middle_index > (chosen_bedseq.indexEnd - int( math.ceil( len(site) / float(2)))):
                            #print "-------------------------------"
                            #print "chosen_middle_index = " + str( chosen_middle_index)
                            #print "chosen_bedseq.indexStart = " + str( chosen_bedseq.indexStart)
                            #print "chosen_bedseq.indexEnd = " + str( chosen_bedseq.indexEnd)
                            #print "int( len( site) / float(2)) = " + str( int( len( site) / float(2)))
                            #print "int( math.ceil( len(site) / float(2))) = " + str( int( math.ceil( len(site) / float(2))))
                            continue
                        chosen_start_index = chosen_middle_index - int(len( site) / float( 2))
                        chosen_end_index = chosen_middle_index + int( math.ceil(len( site) / float( 2)))
                        #chosen_distances.append( int( math.fabs( chosen_middle_index - chosen_bedseq.referenceIndex )))
                        chosen_distances_signed.append( chosen_middle_index - chosen_bedseq.referenceIndex)
                        # Test if any other site previously placed intersect the chosen one
                        intersect = False
                        if chosen_bedseq in implantations.keys():
                            for previous_indexes in implantations[ chosen_bedseq]:
                                if chosen_start_index < previous_indexes[2] and chosen_end_index > previous_indexes[1]:
                                    intersect = True
                                    break
                        # If position is free, add the chosen position to the list of implantations
                        if not intersect:
                            if not chosen_bedseq in implantations.keys():
                                implantations[ chosen_bedseq] = []
                            implantations[ chosen_bedseq].append( (site, chosen_start_index, chosen_end_index))
                            break
                        else:
                            tries += 1
                            if tries > 50:
                                Log.trace("ImplantSitesProcessor.chooseBindingPoints : No place found for site : " + site + ". Bypassing site")
                                break
                    else:
                        print "No bedseq found"
                    
        # Case of unknown distribution
        # ............................
        else:
            raise ExecutionException( "ImplantSitesProcessor.chooseBindingPoint : The chosen distribution mode is unknown '" + distribution_mode)    

        # Compute the histogram of sites distances and graph it
        RSATUtils.outputHistogram( chosen_distances_signed, 5, dir_path, motif.name + "Sites", self.component.pipelineName, "Global distribution of " + motif.name + " sites over peaks",  "Distance from peak maximum",  "Number of occurence" , None)
        

        

        
