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
#   ConservationLimit : the ratio limit to consider a block is conserved
#   DesiredSpeciesList (optionnal): the list of species used in the multiple alignement to research conserved blocks
#   Algorithm : the chosen algorithm. Could be "OccurenceRatio", "InformationRatio" or "null". null mean that the whole
#               peak region is considered as a unique conserved region

## MotifProcessor
#   Method : the method used to obtain the result. Value should be "RSAT compare-matrices" or "TOMTOM"
#   For method "RSAT compare-matrices":
#       MotifDatabasePath : the path to the motifs files
#       MotifDatabaseFileList : list of the name of the motif files to be used as comparison
#       MotifDatabaseFormatList : list of the format of the motif files listed in the MotifDatabaseFileList parameter (one to one)
#       DesiredSpeciesList (optional): the species to take into account in the multiple alignments
#       CorrelationLimit (optional) : the selection minimum threshold of the normalized correlation between reference and query motif (float number beetween 0 and 1)
#       CommandOptions (optional) : specify option to pass to the compare-matrices command
#   For method "TOMTOM":
#       MotifDatabasePath : the path to the motifs files
#       MotifDatabaseFileList : list of the name of the motif files to be used as comparison
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
# No parameters

from processor.BEDProcessor import BEDProcessor
from processor.MAFProcessor import MAFProcessor
from processor.BlockProcessor import BlockProcessor
from processor.MotifProcessor import MotifProcessor
from processor.HistogramProcessor import HistogramProcessor
from processor.CoLocationAnalysisProcessor import CoLocationAnalysisProcessor
from processor.BEDOutputProcessor import BEDOutputProcessor

class OptionManager:

    INPUT_PEAKS = "input_peaks"         ## RESERVE THE LETTER i
    REF_SPECIES = "ref_species"
    ALIGN_SPECIES = "align_species"
    DB_ROOT_PATH = "db_root_path"
    DB_FILE_LIST =  "db_file_list"
    DB_FORMAT_LIST = "db_format_list"
    REF_MOTIF = "ref_motif"

    PIPELINE = "pipeline"   ## RESERVE THE LETTER p
    OUTPUT = "output"       ## RESERVE THE LETTER o
    VERBOSITY = "v"         ## RESERVE THE LETTER v
    FORCE = "force"         ## RESERVE THE LETTER f
    SERVER = "server"       ## RESERVE THE LETTER s
    HELP = "h"              ## RESERVE THE LETTER h
    
    # --------------------------------------------------------------------------------------
    # Apply the provided options to the pipelines components
    @staticmethod
    def applyOptions( pipelines, opts):
    
        for pipeline in pipelines:
    
            for opt, arg in opts:
                
                # Add the BED file path
                if opt == OptionManager.INPUT_PEAKS:
                    OptionManager.addParam( pipeline, "processor.BEDProcessor.BEDProcessor", BEDProcessor.INPUT_BED_FILE_PARAM, arg)

                if opt == OptionManager.REF_SPECIES:
                    OptionManager.addParam( pipeline, "processor.BEDProcessor.BEDProcessor", BEDProcessor.SPECIES_PARAM, arg)
                    OptionManager.addParam( pipeline, "processor.MAFProcessor.MAFProcessor", MAFProcessor.REFERENCE_SPECIES_PARAM, arg)
        
                if opt == OptionManager.ALIGN_SPECIES:
                    OptionManager.addParam( pipeline, "processor.MAFProcessor.MAFProcessor", MAFProcessor.DESIRED_SPECIES_LIST_PARAM, arg)
                    OptionManager.addParam( pipeline, "processor.BlockProcessor.BlockProcessor", BlockProcessor.DESIRED_SPECIES_LIST_PARAM, arg)
                    OptionManager.addParam( pipeline, "processor.MotifProcessor.MotifProcessor", MotifProcessor.DESIRED_SPECIES_LIST_PARAM, arg)
                    
                if opt == OptionManager.DB_ROOT_PATH:
                    OptionManager.addParam( pipeline, "processor.MotifProcessor.MotifProcessor", MotifProcessor.MOTIF_DATABASE_PATH_PARAM, arg)
                    
                if opt == OptionManager.DB_FILE_LIST:
                    OptionManager.addParam( pipeline, "processor.MotifProcessor.MotifProcessor", MotifProcessor.MOTIF_DATABASE_FILE_LIST_PARAM, arg)
                    
                if opt == OptionManager.DB_FORMAT_LIST:
                    OptionManager.addParam( pipeline, "processor.MotifProcessor.MotifProcessor", MotifProcessor.MOTIF_DATABASE_FORMAT_LIST_PARAM, arg)
                
                if opt == OptionManager.REF_MOTIF:
                    OptionManager.addParam( pipeline, "processor.HistogramProcessor.HistogramProcessor", HistogramProcessor.REFERENCE_MOTIF, arg)
                    OptionManager.addParam( pipeline, "processor.CoLocationAnalysisProcessor.CoLocationAnalysisProcessor", CoLocationAnalysisProcessor.REFERENCE_MOTIF_PARAM, arg)
                    OptionManager.addParam( pipeline, "processor.BEDOutputProcessor.BEDOutputProcessor", BEDOutputProcessor.REFERENCE_MOTIF, arg)
            



    # --------------------------------------------------------------------------------------
    # Add the given parameter with the given value to the component with the given name
    @staticmethod
    def addParameter( pipeline, component_name, param_name, param_value):
        
        component = pipeline.getComponent( component_name)
        if component != None:
            component.addParameters( param_name, param_value)