
import os, shutil

from xml.etree.ElementTree import Element
from xml.etree import cElementTree as ET

from processor.Processor import Processor
from processor.io.BedSeqAlignmentStatsCommStruct import BedSeqAlignmentStatsCommStruct
from processor.io.BedSeqAlignmentStatsCommStruct import MotifStatistics

from utils.Constants import Constants
from utils.FileUtils import FileUtils
from utils.MotifUtils import MotifUtils
from utils.RSATUtils import RSATUtils
from utils.log.Log import Log
from utils.exception.ExecutionException import ExecutionException


# This processor aims to create an output file that contain all the information user needs at the end of the analysis
# The output is an XML file that would be displayed using an XSL stylesheet.
#
# Parameters:
#   MotifDatabasePath : the path to the motifs db files
#   MotifDatabaseFileList : list of the name of the motif db files to be used to produce motif logo
#   CustomMotifDatabaseFile (optional) : path to the motif database file provided by the user
#   DisplayLimitValue : the limit valye of the hypergeometric P-value above which an identified motif is not added in the
#                        final result

class FinalOutputProcessor( Processor):
    
    MOTIF_DATABASE_PATH_PARAM = "MotifDatabasePath"
    MOTIF_DATABASE_FILE_LIST_PARAM = "MotifDatabaseFileList"
    CUSTOM_MOTIF_DATABASE_FILE_PARAM = "CustomMotifDatabaseFile"
    DISPLAY_LIMIT_VALUE = "DisplayLimitValue"
    
    LOGOS_DIR_NAME = "Logos"
    
    # --------------------------------------------------------------------------------------
    def __init__( self):
        Processor.__init__( self)
        self.outPath = ""


    # --------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as input 
    # (None if no input CommStruct)
    @staticmethod
    def getInputCommStructClass():
        
        return ( BedSeqAlignmentStatsCommStruct, )


    # --------------------------------------------------------------------------------------
    # Returns the name of the CommStruct class used as output
    # (None if no output CommStruct)
    @staticmethod
    def getOutputCommStructClass():
        
        return None



    # --------------------------------------------------------------------------------------
    # Returns a name that will be used as display name in the user friendly outputs
    @staticmethod
    def getDisplayName():
        
        return "Generation of final result package"
        


    #------------------------------------------------------------------------------------
    # Returns a list of parameters names that are required parameters for the corresponding processor
    @staticmethod
    def getRequiredParameters():
        
        return None



    # --------------------------------------------------------------------------------------
    # Execute the processor
    def execute( self, input_commstructs):
    
        if input_commstructs == None or len( input_commstructs) == 0:
            raise ExecutionException( "FinalOutputProcessor.execute : No inputs")
        
        input_commstruct = input_commstructs[0]
        
        # Get all the pipeline parameters
        parameter_dic = self.collectAllParameters()
        
        # Retrieve the processor parameters
        self.dbPath = self.getParameter( FinalOutputProcessor.MOTIF_DATABASE_PATH_PARAM)
        
        # Retrieve the list of motif database files to use and count motifs in the list
        database_file_line = self.getParameter( FinalOutputProcessor.MOTIF_DATABASE_FILE_LIST_PARAM)
        if database_file_line != None and not database_file_line.isspace():
            file_list = database_file_line.split()
            self.dbFiles = []
            count_motif_in_db = 0
            for file_path in file_list:
                current_path = os.path.join( self.dbPath, file_path)
                current_size = MotifUtils.getMotifsNumberFromTransfac( current_path)
                if current_size != None:
                    count_motif_in_db = count_motif_in_db + current_size
                self.dbFiles.append( current_path)
            parameter_dic[ FinalOutputProcessor.PARAM_MotifDatabaseFileListSize] = str(count_motif_in_db)
        else:
            raise ExecutionException( "FinalOutputProcessor.getMethodParameters : No motif database file specified in parameter '" + FinalOutputProcessor.MOTIF_DATABASE_FILE_LIST_PARAM + "'")

        # Add the custom motif database files if any
        custom_database_file_line = self.getParameter( FinalOutputProcessor.CUSTOM_MOTIF_DATABASE_FILE_PARAM, False)
        if custom_database_file_line != None and not custom_database_file_line.isspace():
            self.dbFiles.append( custom_database_file_line)
                   
        # Retrieve the ID and the name of the reference motif
        motif_name = parameter_dic[ FinalOutputProcessor.PARAM_ReferenceMotif]
        for current_database in self.dbFiles:
            motif_ID = MotifUtils.getMotifIDFromTransfac( motif_name, current_database)
            if motif_ID != None:
                break
        parameter_dic[ FinalOutputProcessor.PARAM_ReferenceMotifID] = motif_ID           
        
        # Retrieve the Hypergeometric e/p-value threshold
        limit_value = self.getParameter( FinalOutputProcessor.DISPLAY_LIMIT_VALUE, False)
        if limit_value == None:
            limit_value = 1.0
        
        # Prepare the processor output dir
        self.outPath = os.path.join( self.component.outputDir, self.component.getComponentPrefix())
        shutil.rmtree( self.outPath, True)
        FileUtils.createDirectory( self.outPath, 0777)
        
        # Copy motif graph and stats files
        analysis = self.AnalyseMotifStats( input_commstruct)
        
        # Create motif logos
        self.createLogos( input_commstruct)
                
        # Copy input BED and custom motif files
        input_bed_file_destination_path = os.path.join(self.outPath, FinalOutputProcessor.PARAM_BEDFile)
        FileUtils.createDirectory( input_bed_file_destination_path)
        FileUtils.copyFile( parameter_dic[ FinalOutputProcessor.PARAM_BEDFile],input_bed_file_destination_path) 
        parameter_dic[ FinalOutputProcessor.PARAM_BEDFile] = os.path.join( input_bed_file_destination_path, os.path.basename(parameter_dic[ FinalOutputProcessor.PARAM_BEDFile]))

        # Copy the custom motif file (if any) to output location and count motifs in this db
        if FinalOutputProcessor.PARAM_CustomMotifDatabaseFile in parameter_dic.keys() and parameter_dic[ FinalOutputProcessor.PARAM_CustomMotifDatabaseFile] != None:
            custom_db_file_destination_path = os.path.join(self.outPath, FinalOutputProcessor.PARAM_CustomMotifDatabaseFile)
            FileUtils.createDirectory( custom_db_file_destination_path)
            FileUtils.copyFile( parameter_dic[ FinalOutputProcessor.PARAM_CustomMotifDatabaseFile],custom_db_file_destination_path) 
            parameter_dic[ FinalOutputProcessor.PARAM_CustomMotifDatabaseFile] = os.path.join( custom_db_file_destination_path, os.path.basename( parameter_dic[ FinalOutputProcessor.PARAM_CustomMotifDatabaseFile]))
            count_custom_motif = MotifUtils.getMotifsNumberFromTransfac( parameter_dic[ FinalOutputProcessor.PARAM_CustomMotifDatabaseFile])
            if( count_custom_motif != None):
                parameter_dic[ FinalOutputProcessor.PARAM_CustomMotifDatabaseFileSize] = str(count_custom_motif)
            else:
                parameter_dic[ FinalOutputProcessor.PARAM_CustomMotifDatabaseFileSize] = "0"
                
        # Output Results
        self.outputClassification( input_commstruct, analysis, limit_value, parameter_dic)
        
        # Copy other information
        FileUtils.copyFile( os.path.join( self.component.outputDir, Constants.PROGRESSION_XSL_FILE), self.outPath) 
        FileUtils.copyFile( os.path.join( self.component.outputDir, Constants.PROGRESSION_XML_FILE), self.outPath)


    # --------------------------------------------------------------------------------------
    # Execute the processor
    def AnalyseMotifStats(self, input_commstruct):
        
        
        # retrieve the global information
        file_pathes = {}
        stats_param = input_commstruct.paramStatistics

        # copy the BED sequences size histogram and graph
        self.copyWorkflowResultFileToFinalOutput( stats_param, BedSeqAlignmentStatsCommStruct.BED_SEQUENCES_SIZE_PATH, FinalOutputProcessor.BED_SEQUENCES_SIZE_PATH_ATT, file_pathes)
        self.copyWorkflowResultFileToFinalOutput( stats_param, BedSeqAlignmentStatsCommStruct.BED_SEQUENCES_SIZE_GRAPH_PATH, FinalOutputProcessor.BED_SEQUENCES_SIZE_GRAPH_PATH_ATT, file_pathes)

        # copy the conserved regions size histogram and graph
        self.copyWorkflowResultFileToFinalOutput( stats_param, BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_SIZE_PATH, FinalOutputProcessor.CONSERVED_BLOCKS_SIZE_PATH_ATT, file_pathes)
        self.copyWorkflowResultFileToFinalOutput( stats_param, BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_SIZE_GRAPH_PATH, FinalOutputProcessor.CONSERVED_BLOCKS_SIZE_GRAPH_PATH_ATT, file_pathes)

        # copy the MSA size histogram and graph
        self.copyWorkflowResultFileToFinalOutput( stats_param, BedSeqAlignmentStatsCommStruct.MSA_SIZE_PATH, FinalOutputProcessor.MSA_SIZE_PATH_ATT, file_pathes)
        self.copyWorkflowResultFileToFinalOutput( stats_param, BedSeqAlignmentStatsCommStruct.MSA_SIZE_GRAPH_PATH, FinalOutputProcessor.MSA_SIZE_GRAPH_PATH_ATT, file_pathes)

        # copy the BED output file
        self.copyWorkflowResultFileToFinalOutput( stats_param, BedSeqAlignmentStatsCommStruct.BED_OUTPUT_PATH, FinalOutputProcessor.BED_OUTPUT_ATT, file_pathes)
        self.copyWorkflowResultFileToFinalOutput( stats_param, BedSeqAlignmentStatsCommStruct.BIGBED_OUTPUT_PATH, FinalOutputProcessor.BIGBED_OUTPUT_ATT, file_pathes)
        
        # Retrieve the motifs statistics information
        identified_motif_names = input_commstruct.motifStatistics.keys()
        family_classification = {}
        motif_classification = {}
        
        for motif_name in identified_motif_names:
            motif_stats = input_commstruct.motifStatistics[ motif_name]
            
            # copy the distance histogram file
            self.copyMotifResultFileToFinalOutput( motif_stats, MotifStatistics.MOTIF_DISTANCE_HISTOGRAM, FinalOutputProcessor.MOTIF_DISTANCE_HISTOGRAM_ATT, file_pathes)
            # copy the distance histogram graph file
            self.copyMotifResultFileToFinalOutput( motif_stats, MotifStatistics.MOTIF_DISTANCE_HISTOGRAM_GRAPH, FinalOutputProcessor.MOTIF_DISTANCE_HISTOGRAM_GRAPH_ATT, file_pathes)
            
            # copy the peak score histogram file
            self.copyMotifResultFileToFinalOutput( motif_stats, MotifStatistics.MOTIF_PEAK_SCORE_HISTOGRAM, FinalOutputProcessor.MOTIF_PEAK_SCORE_HISTOGRAM_ATT, file_pathes)
            # copy the peak score histogram graph file
            self.copyMotifResultFileToFinalOutput( motif_stats, MotifStatistics.MOTIF_PEAK_SCORE_HISTOGRAM_GRAPH, FinalOutputProcessor.MOTIF_PEAK_SCORE_HISTOGRAM_GRAPH_ATT, file_pathes)

            # copy the co-location histogram file
            self.copyMotifResultFileToFinalOutput( motif_stats, MotifStatistics.MOTIF_COLOCATION_HISTOGRAM, FinalOutputProcessor.MOTIF_COLOCATION_HISTOGRAM_ATT, file_pathes)
            # copy the co-location histogram graph file
            self.copyMotifResultFileToFinalOutput( motif_stats, MotifStatistics.MOTIF_COLOCATION_HISTOGRAM_GRAPH, FinalOutputProcessor.MOTIF_COLOCATION_HISTOGRAM_GRAPH_ATT, file_pathes)

            # retrieve the motif classification 
            if motif_stats.hasAttribute( MotifStatistics.MOTIF_FAMILY_RANK):
                motif_family = motif_stats.motifFamily
                family_rank = motif_stats.getAttributeAsint( MotifStatistics.MOTIF_FAMILY_RANK)
                family_classification[ family_rank] = motif_family
                
                motif_rank_in_family = motif_stats.getAttributeAsint( MotifStatistics.MOTIF_RANK_IN_FAMILY)
                if not motif_family in motif_classification.keys():
                    motif_classification[ motif_family] = {}
                motif_classification[ motif_family][ motif_rank_in_family] = motif_name
        
        return ( family_classification, motif_classification, file_pathes)


    # --------------------------------------------------------------------------------------
    # Copy the given file to the out path in the provided sub-directory
    def copyWorkflowResultFileToFinalOutput(self, stats_param, stats_param_attribute, local_attribute, file_pathes):
        
        if stats_param_attribute in stats_param.keys():

            dir_path = os.path.join( self.outPath, stats_param_attribute)
            file_path =  stats_param[ stats_param_attribute]
            
            FileUtils.createDirectory( dir_path)
            FileUtils.copyFile( file_path, dir_path)
            
            final_path = os.path.join( stats_param_attribute, os.path.basename( file_path))
            
            file_pathes[ local_attribute] = final_path        



    # --------------------------------------------------------------------------------------
    # Copy the given file to the out path in the provided sub-directory
    def copyMotifResultFileToFinalOutput(self, motif_stats, motifStatistics_attribute, local_attribute, file_pathes):
        
        
        if motif_stats.hasAttribute( motifStatistics_attribute):
            if not local_attribute in file_pathes.keys():
                file_pathes[ local_attribute] = {}

            dir_path = os.path.join( self.outPath, motifStatistics_attribute)
            file_path =  motif_stats.getAttribute( motifStatistics_attribute)
            
            FileUtils.createDirectory( dir_path)
            FileUtils.copyFile( file_path, dir_path)
            
            final_path = os.path.join( motifStatistics_attribute, os.path.basename( file_path))
            
            file_pathes[ local_attribute][ motif_stats.motifName] = final_path     



    # --------------------------------------------------------------------------------------
    # Create the logo image file for each motif in the motif statistics
    def createLogos(self, input_commstruct):
        
        
        db_file_path = []
        for index in range( len( self.dbFiles)):
            db_file_path.append( os.path.join( self.dbPath, self.dbFiles[index]));
        
        motif_name_list = input_commstruct.motifStatistics.keys()
        motif_definition = MotifUtils.getMotifsDefinitionFromTF( motif_name_list, db_file_path)
        logos_path = os.path.join( self.outPath, FinalOutputProcessor.LOGOS_DIR_NAME)
        FileUtils.createDirectory( logos_path)
        
        for motif_name in motif_name_list:
            if motif_name in motif_definition.keys():
                file_name = motif_name + ".tf"
                def_file_path = os.path.join( logos_path, file_name)
                def_file = open( def_file_path, "w")
                for line in motif_definition[ motif_name]:
                    def_file.write( line)
                    def_file.flush
                def_file.close()
                RSATUtils.createLogoFromTF( logos_path, file_name, motif_name)
            else:
                Log.log( "FinalOutputProcessor.createLogos : No definition found to create logo for motif : " + motif_name)


    # --------------------------------------------------------------------------------------
    # Output the motif classification
    def outputClassification(self, input_commstruct, analysis, limit_value, parameter_dic):

        try:
            # Create and write to file the XML element
            root_element = self.toXML( input_commstruct, analysis, limit_value, parameter_dic)
            self.indent( root_element, 0)
            # Output the XML to file
            doc = ET.ElementTree( root_element)
            classification_file_path = os.path.join( self.outPath, self.component.pipelineName + "_MotifClassification.xml")
            outfile = open( classification_file_path, 'w')
            outfile.write('<?xml version="1.0" encoding="utf-8"?>\n')
            outfile.write('<?xml-stylesheet type="text/xsl" href="classification.xsl"?>\n')
            doc.write( outfile)
            outfile.close()
            # Copy the XSL file in the same directory than the XML
            shutil.copy( os.path.join( self.component.getParameter( Constants.INSTALL_DIR_PARAM), "resources/xsl/classification/classification.xsl"), self.outPath)
            shutil.copy( os.path.join( self.component.getParameter( Constants.INSTALL_DIR_PARAM), "resources/xsl/classification/jquery.dataTables.js"), self.outPath)
            shutil.copy( os.path.join( self.component.getParameter( Constants.INSTALL_DIR_PARAM), "resources/xsl/classification/peak-footprints.css"), self.outPath)
        except IOError, exce:
            Log.log( "ClassificationProcessor.outputClassification : Unable to write classification to XML file. From:\n\t---> " + str( exce))



    # --------------------------------------------------------------------------------------
    # Write the Classification to XML file
    def toXML( self, input_commstruct, analysis, limit_value, parameter_dic):
        
        # Retrieve the data from the analyzed statistics
        classified_families = analysis[ 0]
        motif_classification = analysis[ 1]
        file_pathes = analysis[2]
        
        
        # Create the root element with its attributes
        classification_element = Element( FinalOutputProcessor.CLASSIFICATION_TAG)
        classification_element.attrib[ FinalOutputProcessor.PIPELINE_NAME_ATT] = self.component.pipelineName
        classification_element.attrib[ FinalOutputProcessor.REFERENCE_SPECIES_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.REFERENCE_SPECIES]
        classification_element.attrib[ FinalOutputProcessor.ALIGNED_SPECIES_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.ALIGNED_SPECIES]
        classification_element.attrib[ FinalOutputProcessor.REFERENCE_MOTIF_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.REFERENCE_MOTIF]

        classification_element.attrib[ FinalOutputProcessor.BEDSEQUENCES_NUMBER_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.BEDSEQUENCES_NUMBER]
        classification_element.attrib[ FinalOutputProcessor.BEDSEQUENCES_MIN_SIZE_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.BEDSEQUENCES_MIN_SIZE]
        classification_element.attrib[ FinalOutputProcessor.BEDSEQUENCES_MAX_SIZE_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.BEDSEQUENCES_MAX_SIZE]
        classification_element.attrib[ FinalOutputProcessor.BEDSEQUENCES_MEAN_SIZE_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.BEDSEQUENCES_MEAN_SIZE]
        classification_element.attrib[ FinalOutputProcessor.BEDSEQUENCES_TOTAL_SIZE_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.BEDSEQUENCES_TOTAL_SIZE]

        classification_element.attrib[ FinalOutputProcessor.MSA_NUMBER_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.MSA_NUMBER]
        classification_element.attrib[ FinalOutputProcessor.MSA_MIN_SIZE_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.MSA_MIN_SIZE]
        classification_element.attrib[ FinalOutputProcessor.MSA_MAX_SIZE_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.MSA_MAX_SIZE]
        classification_element.attrib[ FinalOutputProcessor.MSA_MEAN_SIZE_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.MSA_MEAN_SIZE]
        classification_element.attrib[ FinalOutputProcessor.MSA_TOTAL_SIZE_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.MSA_TOTAL_SIZE]

        classification_element.attrib[ FinalOutputProcessor.CONSERVED_BLOCKS_NUMBER_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_NUMBER]
        classification_element.attrib[ FinalOutputProcessor.CONSERVED_BLOCKS_MIN_SIZE_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_MIN_SIZE]
        classification_element.attrib[ FinalOutputProcessor.CONSERVED_BLOCKS_MAX_SIZE_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_MAX_SIZE]
        classification_element.attrib[ FinalOutputProcessor.CONSERVED_BLOCKS_MEAN_SIZE_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_MEAN_SIZE]
        classification_element.attrib[ FinalOutputProcessor.CONSERVED_BLOCKS_TOTAL_SIZE_ATT] = input_commstruct.paramStatistics[ BedSeqAlignmentStatsCommStruct.CONSERVED_BLOCKS_TOTAL_SIZE]

        classification_element.attrib[ FinalOutputProcessor.BED_SEQUENCES_SIZE_PATH_ATT] = file_pathes[ FinalOutputProcessor.BED_SEQUENCES_SIZE_PATH_ATT]
        classification_element.attrib[ FinalOutputProcessor.BED_SEQUENCES_SIZE_GRAPH_PATH_ATT] = file_pathes[ FinalOutputProcessor.BED_SEQUENCES_SIZE_GRAPH_PATH_ATT]

        # Add all other parameters values
        for param_name in parameter_dic.keys():
            if not param_name in classification_element.attrib.keys():  
                classification_element.attrib[ param_name] = parameter_dic[ param_name]

        # Insert the path to the BED sequences sizes histogram and graph
        self.addFilePathAttribute( classification_element, FinalOutputProcessor.BED_SEQUENCES_SIZE_PATH_ATT, file_pathes)
        self.addFilePathAttribute( classification_element, FinalOutputProcessor.BED_SEQUENCES_SIZE_GRAPH_PATH_ATT, file_pathes)

        # Insert the path to the Conserved Regions sizes histogram and graph if any
        self.addFilePathAttribute( classification_element, FinalOutputProcessor.CONSERVED_BLOCKS_SIZE_PATH_ATT, file_pathes)
        self.addFilePathAttribute( classification_element, FinalOutputProcessor.CONSERVED_BLOCKS_SIZE_GRAPH_PATH_ATT, file_pathes)

        # Insert the path to the MSA sizes histogram and graph if any
        self.addFilePathAttribute( classification_element, FinalOutputProcessor.MSA_SIZE_PATH_ATT, file_pathes)
        self.addFilePathAttribute( classification_element, FinalOutputProcessor.MSA_SIZE_GRAPH_PATH_ATT, file_pathes)

        # Insert the path to the BED output file if any
        self.addFilePathAttribute( classification_element, FinalOutputProcessor.BED_OUTPUT_ATT, file_pathes)
        self.addFilePathAttribute( classification_element, FinalOutputProcessor.BIGBED_OUTPUT_ATT, file_pathes)

        # Create a root son element for each family
        for family_rank in sorted( classified_families):
            # fill the family element attributes
            family_element = Element( FinalOutputProcessor.FAMILY_TAG)
            classification_element.append( family_element)
            
            family = classified_families[ family_rank]
            family_element.attrib[ FinalOutputProcessor.FAMILY_NAME_ATT] = family
            
            # create a family son element for each motif in a family
            for motif_rank in sorted( motif_classification[ family].keys()):
                
                motif_name = motif_classification[ family][ motif_rank]
                motif_stats = input_commstruct.motifStatistics[ motif_name]
                
                # If the motif has its hypergeometric p-value above the limit, it is ignored
                #if not motif_stats.hasAttribute( MotifStatistics.MOTIF_HYP_PVALUE) or motif_stats.getAttributeAsfloat( MotifStatistics.MOTIF_HYP_PVALUE) > limit_value:
                #    continue
                # If the motif has its hypergeometric e-value above the limit, it is ignored
                if not motif_stats.hasAttribute( MotifStatistics.MOTIF_HYP_EVALUE) or motif_stats.getAttributeAsfloat( MotifStatistics.MOTIF_HYP_EVALUE) > limit_value:
                    continue
                
                # Create the motif XML element
                motif_element = Element( FinalOutputProcessor.MOTIF_TAG)
                family_element.append( motif_element)
                
                # fill the motif element attributes
                motif_element.attrib[ FinalOutputProcessor.MOTIF_NAME_ATT] = motif_name
                motif_element.attrib[ FinalOutputProcessor.MOTIF_FAMILY_ATT] = family
                motif_element.attrib[ FinalOutputProcessor.MOTIF_ID_ATT] = motif_stats.motifID
                motif_element.attrib[ FinalOutputProcessor.MOTIF_CLASS_ATT] = motif_stats.motifClass
                motif_element.attrib[ FinalOutputProcessor.MOTIF_TYPE_ATT] = motif_stats.motifType

                motif_element.attrib[ FinalOutputProcessor.MOTIF_RANK_ATT] = motif_stats.getAttribute( MotifStatistics.MOTIF_RANK)

                # fill the motif hit score attribute
                motif_element.attrib[ FinalOutputProcessor.MOTIF_HIT_SCORE_ATT] = motif_stats.getAttribute( MotifStatistics.MOTIF_HIT_SCORE)
                               
                # fill the motifhypergeometric score and p-value attribute
                if motif_stats.hasAttribute( MotifStatistics.MOTIF_HYP_HIT_SCORE):
                    motif_element.attrib[ FinalOutputProcessor.MOTIF_HYP_HIT_SCORE_ATT] = motif_stats.getAttribute( MotifStatistics.MOTIF_HYP_HIT_SCORE)
                if motif_stats.hasAttribute( MotifStatistics.MOTIF_HYP_PVALUE):
                    motif_element.attrib[ FinalOutputProcessor.MOTIF_HYP_PVALUE_ATT] = motif_stats.getAttribute( MotifStatistics.MOTIF_HYP_PVALUE)
                if motif_stats.hasAttribute( MotifStatistics.MOTIF_HYP_EVALUE):
                    motif_element.attrib[ FinalOutputProcessor.MOTIF_HYP_EVALUE_ATT] = motif_stats.getAttribute( MotifStatistics.MOTIF_HYP_EVALUE)

                # fill the motif chi2 and chi2 p-value attributes
                if motif_stats.hasAttribute( MotifStatistics.MOTIF_CHI2):
                    motif_element.attrib[ FinalOutputProcessor.MOTIF_CHI2_ATT] = motif_stats.getAttribute( MotifStatistics.MOTIF_CHI2)
                if motif_stats.hasAttribute( MotifStatistics.MOTIF_CHI2_PVALUE):
                    motif_element.attrib[ FinalOutputProcessor.MOTIF_CHI2_PVALUE_ATT] = motif_stats.getAttribute( MotifStatistics.MOTIF_CHI2_PVALUE)


                # fill the ratio overlap attribute
                motif_element.attrib[ FinalOutputProcessor.MOTIF_RATIO_HOMOLOCATION] = str( int( motif_stats.getAttributeAsfloat( MotifStatistics.MOTIF_RATIO_HOMOLOCATION)*1000.0)/float(10)) + "%"
                
                # Fill the motif logo and matrix attributes
                motif_element.attrib[ FinalOutputProcessor.MOTIF_LOGO_ATT] = os.path.join( FinalOutputProcessor.LOGOS_DIR_NAME, motif_name + "_m1.png")
                motif_element.attrib[ FinalOutputProcessor.MOTIF_MATRIX_ATT] = os.path.join( FinalOutputProcessor.LOGOS_DIR_NAME, motif_name + ".tf")

                # fill the motif element graphs
                self.addMotifFilePathAttribute( motif_element, FinalOutputProcessor.MOTIF_DISTANCE_HISTOGRAM_ATT, motif_name, file_pathes)
                self.addMotifFilePathAttribute( motif_element, FinalOutputProcessor.MOTIF_DISTANCE_HISTOGRAM_GRAPH_ATT, motif_name, file_pathes)

                self.addMotifFilePathAttribute( motif_element, FinalOutputProcessor.MOTIF_PEAK_SCORE_HISTOGRAM_ATT, motif_name, file_pathes)
                self.addMotifFilePathAttribute( motif_element, FinalOutputProcessor.MOTIF_PEAK_SCORE_HISTOGRAM_GRAPH_ATT, motif_name, file_pathes)
                
                self.addMotifFilePathAttribute( motif_element, FinalOutputProcessor.MOTIF_COLOCATION_HISTOGRAM_ATT, motif_name, file_pathes)
                self.addMotifFilePathAttribute( motif_element, FinalOutputProcessor.MOTIF_COLOCATION_HISTOGRAM_GRAPH_ATT, motif_name, file_pathes)

        return classification_element


    # --------------------------------------------------------------------------------------
    # Add to the element an attribute containing the file path associated to given attribute
    def addFilePathAttribute(self, element, attribute, file_pathes):
        
        if attribute in file_pathes.keys():
            element.attrib[ attribute] = file_pathes[ attribute]



    # --------------------------------------------------------------------------------------
    # Add to the element an attribute containing the file path associated to given motif and attribute
    def addMotifFilePathAttribute(self, element, attribute, name, file_pathes):
        
        if attribute in file_pathes.keys():
            if name in file_pathes[ attribute].keys():
                element.attrib[ attribute] = file_pathes[ attribute][ name]


    # --------------------------------------------------------------------------------------
    # Collect all the parameters used by the various processors of the executed pipeline
    def collectAllParameters(self):
        
        option_dic ={}
        
        for component in self.pipeline.getComponentList():
            parameter_dic = component.getParameterDic()
            for parameter_name in parameter_dic.keys():
                option_dic[ parameter_name] = parameter_dic[ parameter_name]
                
        return option_dic


    # --------------------------------------------------------------------------------------
    # Constants used for the XML file creation (tag and attribute names)

    CLASSIFICATION_TAG = "classification"
    PIPELINE_NAME_ATT = "name"
    REFERENCE_SPECIES_ATT = "referenceSpecies"
    ALIGNED_SPECIES_ATT = "alignedSpecies"
    REFERENCE_MOTIF_ATT = "referenceMotif"
    
    PARAM_ReferenceMotif = "ReferenceMotif"
    PARAM_ReferenceMotifID = "ReferenceMotifID"
    PARAM_BEDFile = "BEDFile"
    PARAM_ResiduConservationLimit= "ResiduConservationLimit"
    PARAM_WindowSize= "WindowSize"
    PARAM_WindowConservationLimit= "WindowConservationLimit"
    PARAM_CustomMotifDatabaseFile= "CustomMotifDatabaseFile"
    PARAM_CustomMotifDatabaseFileSize= "CustomMotifDatabaseFileSize"
    PARAM_MotifDatabasePath= "MotifDatabasePath"
    PARAM_MotifDatabaseFileList= "MotifDatabaseFileList"
    PARAM_MotifDatabaseFileListSize= "MotifDatabaseFileListSize"
    PARAM_MaxHypergeometricEValue= "MaxHypergeometricEValue"
    PARAM_MaxChi2EValue= "MaxChi2EValue"
    PARAM_MaxMotifByFamily= "MaxMotifByFamily"
    PARAM_MaxMotifNumber= "MaxMotifNumber"
    
    BEDSEQUENCES_NUMBER_ATT = "bedSequencesNumber"
    BEDSEQUENCES_MIN_SIZE_ATT = "bedSequencesMinSize"
    BEDSEQUENCES_MAX_SIZE_ATT = "bedSequencesMaxSize"
    BEDSEQUENCES_MEAN_SIZE_ATT = "bedSequencesMeanSize"
    BEDSEQUENCES_TOTAL_SIZE_ATT = "bedSequencesTotalSize"

    MSA_NUMBER_ATT = "MSANumber"
    MSA_MIN_SIZE_ATT = "MSAMinSize"
    MSA_MAX_SIZE_ATT = "MSAMaxSize"
    MSA_MEAN_SIZE_ATT = "MSAMeanSize"
    MSA_TOTAL_SIZE_ATT = "MSATotalSize"

    CONSERVED_BLOCKS_NUMBER_ATT = "conservedBlocksNumber"
    CONSERVED_BLOCKS_MIN_SIZE_ATT = "conservedBlocksMinSize"
    CONSERVED_BLOCKS_MAX_SIZE_ATT = "conservedBlocksMaxSize"
    CONSERVED_BLOCKS_MEAN_SIZE_ATT = "conservedBlocksMeanSize"
    CONSERVED_BLOCKS_TOTAL_SIZE_ATT = "conservedBlocksTotalSize"
    
    BED_SEQUENCES_SIZE_PATH_ATT = "BEDSequencesSizeHistogram"
    BED_SEQUENCES_SIZE_GRAPH_PATH_ATT = "BEDSequencesSizeHistogramGraph"
    
    CONSERVED_BLOCKS_SIZE_PATH_ATT = "ConservedBlocksSizeHistogram"
    CONSERVED_BLOCKS_SIZE_GRAPH_PATH_ATT = "ConservedBlocksSizeGraph"
    MSA_SIZE_PATH_ATT = "MSASizeHistogram"
    MSA_SIZE_GRAPH_PATH_ATT = "MSASizeHistogramGraph"
    BED_OUTPUT_ATT = "bedOutput"
    BIGBED_OUTPUT_ATT = "bigbedOutput"
    
    FAMILY_TAG = "family"
    FAMILY_NAME_ATT = "name"
    FAMILY_SCORE_ATT = "score"
    
    MOTIF_TAG = "motif"
    MOTIF_NAME_ATT = "name"
    MOTIF_FAMILY_ATT = "family"
    MOTIF_ID_ATT = "id"
    MOTIF_CLASS_ATT = "class"
    MOTIF_TYPE_ATT = "type"
    MOTIF_RANK_ATT = "rank"
    MOTIF_HIT_SCORE_ATT = "hitscore"
    MOTIF_CHI2_ATT ="chi2"
    MOTIF_CHI2_PVALUE_ATT ="chi2pvalue"
    MOTIF_HYP_HIT_SCORE_ATT = "hyphitsscore"
    MOTIF_HYP_PVALUE_ATT ="hyppvalue"
    MOTIF_HYP_EVALUE_ATT ="hypevalue"
    MOTIF_RATIO_HOMOLOCATION = "ratioHomoLocation"
    MOTIF_LOGO_ATT = "logo"
    MOTIF_MATRIX_ATT = "matrix"

    MOTIF_DISTANCE_HISTOGRAM_ATT = "distanceHistogram"
    MOTIF_DISTANCE_HISTOGRAM_GRAPH_ATT = "distanceHistogramGraph"
    MOTIF_PEAK_SCORE_HISTOGRAM_ATT ="peakScoreHistogram"
    MOTIF_PEAK_SCORE_HISTOGRAM_GRAPH_ATT = "peakScoreHistogramGraph"
    MOTIF_COLOCATION_HISTOGRAM_ATT = "coLocationHistogram"
    MOTIF_COLOCATION_HISTOGRAM_GRAPH_ATT = "coLocationHistogramGraph"

