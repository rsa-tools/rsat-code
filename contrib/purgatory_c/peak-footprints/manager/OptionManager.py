#!/usr/bin/python

## BEDProcessor
#   BEDFile : the path to the BED file
#   Species : the species of the described sequences
#   PeakFile (optional): the path to the file giving information on the maximum of the peak 
#                        for each sequence in the BED file. Used in case of BED file represent ChIP-seq peaks result.
#   PeakNumber (optional): indicates the number of peak to consider in the given BED file. Peaks are uniformly 
#                          randomly chosen through all the peaks.

## MAFProcessor
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

## BlockProcessor
#   WindowSize : the size of the window used to parse the sequences when searching conserved blocks
#   ResiduConservationLimit : the ratio limit to consider a residu is conserved
#   WindowConservationLimit : the ratio limit to consider a block is conserved
#   DesiredSpeciesList (optional): the list of species used in the multiple alignement to research conserved blocks
#   Algorithm : the chosen algorithm. Could be "OccurenceRatio", "InformationRatio" or "null". null mean that the whole
#               peak region is considered as a unique conserved region

## MotifProcessor
#   Method : the method used to obtain the result. Value should be "RSAT compare-matrices" or "TOMTOM"
#   For method "RSAT compare-matrices":
#       MotifDatabasePath : the path to the motifs files
#       MotifDatabaseFileList : list of the name of the motif files to be used as comparison
#       MotifDatabaseFormatList : list of the format of the motif files listed in the MotifDatabaseFileList parameter (one to one)
#       CustomMotifDatabaseFile (optional) : path to the motif database file provided by the user
#       CustomMotifDatabaseFormat (optional) : format of the motif database file provided by the user
#       DesiredSpeciesList (optional): the species to take into account in the multiple alignments
#       CorrelationLimit (optional) : the selection minimum threshold of the normalized correlation between reference and query motif (float number beetween 0 and 1)
#       CommandOptions (optional) : specify option to pass to the compare-matrices command
#   For method "TOMTOM":
#       MotifDatabasePath : the path to the motifs files
#       MotifDatabaseFileList : list of the name of the motif files to be used as comparison
#       CustomMotifDatabaseFile (optional) : path to the motif database file provided by the user
#       DesiredSpeciesList (optional): the species to take into account in the multiple alignments
#       CommandOptions (optional) : specify option to pass to the TOMTOM command

## HistogramProcessor
#   HistogramInterval : the interval length of value used to compute the histogram
#   ReferenceMotif    : the motif used as reference in the peaks

## ClassificationProcessor
#   MaxMotifNumber (optional): the maximum number of motifs to display in the result
#   MaxMotifByFamily (optional): the maximum number of motifs to display in the same family
#   MaxHypergeometricEValue (optional): the maximum value of the hypergeometric e-value the displayed motif must have
#   MaxChi2EValue (optional): the maximum value of the chi2 e-value the displayed motif must have

## ColocationAnalysisProcessor
#   ReferenceMotif = Motif that can be used as reference (should be the motif of the TF studied by the chIP-seq)
#   HistogramInterval = Class size of the computed histograms
#   MaximalDistance = The distance between the reference motif and the considered motif above which the hit is not reported

## BEDOutputProcessor
#   ReferenceMotif: the motif used as reference
#   Method : the method used to assign a color to the motifs. Value should be "family" or "score"
#            - the "family" method assign a color for each motif family. Note that following the UCSC recommendation only
#              5 colors are assigned , so some families can have the same colors.
#            - the "score" method assign a color to the motif according to its score. The greater the score, the more the color
#              got from the blue to red in the rainbow order
#            Note: in both cases, the reference motif is assigned to black color
#   ScoreMin : the minimum value the score can reach
#   ScoreMax : the maximum value the score can reach

## FinalOutputProcessor
#   MotifDatabasePath : the path to the motifs db files
#   MotifDatabaseFileList : list of the name of the motif db files to be used to produce motif logo
#   CustomMotifDatabaseFile (optional) : path to the motif database file provided by the user
#   DisplayLimitValue : the limit value of the hypergeometric P-value above which an identified motif is not added in the
#                        final result

from processor.BEDProcessor import BEDProcessor
from processor.MAFProcessor import MAFProcessor
from processor.BlockProcessor import BlockProcessor
from processor.MotifProcessor import MotifProcessor
from processor.HistogramProcessor import HistogramProcessor
from processor.CoLocationAnalysisProcessor import CoLocationAnalysisProcessor
from processor.BEDOutputProcessor import BEDOutputProcessor
from processor.FinalOutputProcessor import FinalOutputProcessor
from processor.ClassificationProcessor import ClassificationProcessor

class OptionManager:

    PIPELINE_NAME = "pipeline_name"
    INPUT_PEAKS = "input_peaks"
    PEAK_NUMBER = "peak_number"
    REF_SPECIES = "ref_species"
    MAF_FILE_PATH = "maf_path"
    ALIGN_SPECIES = "align_species"
    DB_ROOT_PATH = "db_root_path"
    DB_FILE_LIST = "db_file_list"
    DB_FORMAT_LIST = "db_format_list"
    CUSTO_DB_FILE_PATH = "custo_db_file_path"
    REF_MOTIF = "ref_motif"
    WINDOW_SIZE = "window_size"
    RESIDU_CONSERVATION_LIMIT = "conservation_threshold"
    WINDOW_CONSERVATION_LIMIT = "window_conservation_threshold"
    MAX_MOTIF_NUMBER = "max_motif_number"
    MAX_MOTIF_BY_FAMILY = "max_motif_by_family"
    MAX_HYP_EVALUE = "max_hyp_pvalue"
    MAX_CHI2_EVALUE = "max_chi2_pvalue"
    EXTENSION_5P = "5pext";
    EXTENSION_3P = "3pext";

    PIPELINE = "pipeline"   
    OUTPUT = "output"
    RSAT = "rsat"       
    VERBOSITY = "v"         
    FORCE = "force"        
    SERVER = "server"       
    HELP = "h"              
    
    # --------------------------------------------------------------------------------------
    # Apply the provided options to the pipelines components
    @staticmethod
    def applyOptions( pipelines, opts):
    
        for pipeline in pipelines:
    
            for opt, arg in opts:
                
                # Remove the dashed at the beginning of the option tag
                while opt[0] == "-":
                    opt = opt[1:]
                
                print "option : " + str( opt) + " = " + str( arg)
                
                # Set the name of the pipeline (to the pipeline itself and to its components)
                if opt == OptionManager.PIPELINE_NAME:
                    pipeline.name = arg
                    for component in pipeline.componentList:
                        component.pipelineName = arg
                
                # Add the BED file path
                elif opt == OptionManager.INPUT_PEAKS:
                    OptionManager.addParam( pipeline, "processor.BEDProcessor.BEDProcessor", BEDProcessor.INPUT_BED_FILE_PARAM, arg)
                
                # Add the value of the number of peaks taken into account
                elif opt == OptionManager.PEAK_NUMBER:
                    OptionManager.addParam( pipeline, "processor.BEDProcessor.BEDProcessor", BEDProcessor.PEAK_NUMBER, arg)
                
                # Add the reference species
                elif opt == OptionManager.REF_SPECIES:
                    OptionManager.addParam( pipeline, "processor.BEDProcessor.BEDProcessor", BEDProcessor.REFERENCE_SPECIES_PARAM, arg)
                    OptionManager.addParam( pipeline, "processor.MAFProcessor.MAFProcessor", MAFProcessor.REFERENCE_SPECIES_PARAM, arg)
                    OptionManager.addParam( pipeline, "processor.BlockProcessor.BlockProcessor", BlockProcessor.REFERENCE_SPECIES_PARAM, arg)
                    OptionManager.addParam( pipeline, "processor.MotifProcessor.MotifProcessor", MotifProcessor.REFERENCE_SPECIES_PARAM, arg)
       
                # Add the path to the MAF files
                elif opt == OptionManager.MAF_FILE_PATH:
                    OptionManager.addParam( pipeline, "processor.MAFProcessor.MAFProcessor", MAFProcessor.INPUT_MAF_FILE_PARAM, arg)
       
                # Add the list of species to be aligned with the reference one
                elif opt == OptionManager.ALIGN_SPECIES:
                    OptionManager.addParam( pipeline, "processor.MAFProcessor.MAFProcessor", MAFProcessor.DESIRED_SPECIES_LIST_PARAM, arg)
                    OptionManager.addParam( pipeline, "processor.BlockProcessor.BlockProcessor", BlockProcessor.DESIRED_SPECIES_LIST_PARAM, arg)
                    OptionManager.addParam( pipeline, "processor.MotifProcessor.MotifProcessor", MotifProcessor.DESIRED_SPECIES_LIST_PARAM, arg)
                
                # Add root path to the TF databases
                elif opt == OptionManager.DB_ROOT_PATH:
                    OptionManager.addParam( pipeline, "processor.MotifProcessor.MotifProcessor", MotifProcessor.MOTIF_DATABASE_PATH_PARAM, arg)
                    OptionManager.addParam( pipeline, "processor.FinalOutputProcessor.FinalOutputProcessor", FinalOutputProcessor.MOTIF_DATABASE_PATH_PARAM, arg)
                
                # Add TF databases file name list (names could contains missing path pieces from root path)
                elif opt == OptionManager.DB_FILE_LIST:
                    OptionManager.addParam( pipeline, "processor.MotifProcessor.MotifProcessor", MotifProcessor.MOTIF_DATABASE_FILE_LIST_PARAM, arg)
                    OptionManager.addParam( pipeline, "processor.FinalOutputProcessor.FinalOutputProcessor", FinalOutputProcessor.MOTIF_DATABASE_FILE_LIST_PARAM, arg)
                
                # Add TF databases format list
                elif opt == OptionManager.DB_FORMAT_LIST:
                    OptionManager.addParam( pipeline, "processor.MotifProcessor.MotifProcessor", MotifProcessor.MOTIF_DATABASE_FORMAT_LIST_PARAM, arg)
                    #OptionManager.addParam( pipeline, "processor.FinalOutputProcessor.FinalOutputProcessor", FinalOutputProcessor.MOTIF_DATABASE_FORMAT_LIST_PARAM, arg)
                
                # Add Custom TF database file path
                elif opt == OptionManager.CUSTO_DB_FILE_PATH:
                    OptionManager.addParam( pipeline, "processor.MotifProcessor.MotifProcessor", MotifProcessor.CUSTOM_MOTIF_DATABASE_FILE_PARAM, arg)
                    OptionManager.addParam( pipeline, "processor.FinalOutputProcessor.FinalOutputProcessor", FinalOutputProcessor.CUSTOM_MOTIF_DATABASE_FILE_PARAM, arg)
                
                # Add reference motif
                elif opt == OptionManager.REF_MOTIF:
                    OptionManager.addParam( pipeline, "processor.HistogramProcessor.HistogramProcessor", HistogramProcessor.REFERENCE_MOTIF, arg)
                    OptionManager.addParam( pipeline, "processor.CoLocationAnalysisProcessor.CoLocationAnalysisProcessor", CoLocationAnalysisProcessor.REFERENCE_MOTIF_PARAM, arg)
                    OptionManager.addParam( pipeline, "processor.BEDOutputProcessor.BEDOutputProcessor", BEDOutputProcessor.REFERENCE_MOTIF, arg)

                # Add the initial size of the sliding window used for conservation detection
                elif opt == OptionManager.WINDOW_SIZE:
                    OptionManager.addParam( pipeline, "processor.BlockProcessor.BlockProcessor", BlockProcessor.WINDOW_SIZE_PARAM, arg)
                    
                # Add the value of the threshold indicating if a residue is conserved
                elif opt == OptionManager.RESIDU_CONSERVATION_LIMIT:
                    OptionManager.addParam( pipeline, "processor.BlockProcessor.BlockProcessor", BlockProcessor.RESIDU_CONSERVATION_LIMIT_PARAM, arg)
                
                # Add the value of the threshold indicating if a block is conserved
                elif opt == OptionManager.WINDOW_CONSERVATION_LIMIT:
                    OptionManager.addParam( pipeline, "processor.BlockProcessor.BlockProcessor", BlockProcessor.WINDOW_CONSERVATION_LIMIT_PARAM, arg)
                    
                # Add the maximal number of motif finally reported
                elif opt == OptionManager.MAX_MOTIF_NUMBER:
                    OptionManager.addParam( pipeline, "processor.ClassificationProcessor.ClassificationProcessor", ClassificationProcessor.MAX_MOTIF_NUMBER, arg)
                    
                # Add the maximal number of motif reported in each motif family
                elif opt == OptionManager.MAX_MOTIF_BY_FAMILY:
                    OptionManager.addParam( pipeline, "processor.ClassificationProcessor.ClassificationProcessor", ClassificationProcessor.MAX_MOTIF_BY_FAMILY, arg)

                # Add the threshold on the hypergeometric e-value
                elif opt == OptionManager.MAX_HYP_EVALUE:
                    OptionManager.addParam( pipeline, "processor.ClassificationProcessor.ClassificationProcessor", ClassificationProcessor.MAX_HYP_EVALUE, arg)

                # Add the threshold on the chi2 e-value
                elif opt == OptionManager.MAX_CHI2_EVALUE:
                    OptionManager.addParam( pipeline, "processor.ClassificationProcessor.ClassificationProcessor", ClassificationProcessor.MAX_CHI2_EVALUE, arg)

                # Add the peak 5' extensions
                elif opt == OptionManager.EXTENSION_5P:
                    OptionManager.addParam( pipeline, "processor.BEDProcessor.BEDProcessor", BEDProcessor.EXTENSION_5P, arg)

                # Add the peak 3' extensions
                elif opt == OptionManager.EXTENSION_3P:
                    OptionManager.addParam( pipeline, "processor.BEDProcessor.BEDProcessor", BEDProcessor.EXTENSION_3P, arg)

    # --------------------------------------------------------------------------------------
    # Add the given parameter with the given value to the component with the given name
    @staticmethod
    def addParam( pipeline, component_name, param_name, param_value):
        
        component = pipeline.getComponent( component_name)
        if component != None:
            component.addParameters( param_name, param_value)
            
            
    # --------------------------------------------------------------------------------------
    # Add the given parameter value to the current parameter value of the component with the given name
    @staticmethod
    def extendParam( pipeline, component_name, param_name, param_value):
        
        component = pipeline.getComponent( component_name)
        if component != None:
            old_param_value = component.getParameter( param_name, False)
            if old_param_value != None:
                component.addParameters( param_name, old_param_value + param_value)
