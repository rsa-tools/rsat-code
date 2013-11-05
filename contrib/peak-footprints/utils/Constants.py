
class Constants:

    # --------------------------------------------------------------------------------------
    #Init parameters
    PROJECT_NAME = "peak-footprints"
    RSAT_PATH_ENV_VAR = "RSAT"
    PROJECT_PATH_IN_RSAT = "contrib/peak-footprints"
    PROJECT_INSTALL_PATH_ENV_VAR = "FOOTPRINTING_PATH"
    MANAGER_CONFIG_FILE_NAME = "manager.props"
    OUTPUT_DIR_NAME = "output"
    QUEUE_DIR_NAME = "queue"
    LOG_FILE_NAME = "logs.txt"
    TRACE_FILE_NAME = "traces.txt"
    SERVER_QUEUE_FILE_NAME = "serverQueue.txt"
    

    # --------------------------------------------------------------------------------------
    #Manager Config parameters
    INSTALL_DIR_PARAM = "install_dir"
    BASE_OUTPUT_DIR_PARAM = "base_output_dir"
    SERVER_QUEUE_DIR_PARAM = "server_queue_dir"
    LISTENING_DIR_PARAM = "listening_dir"
    OUTPUT_DIR_PARAM = "output_dir"
    QUEUE_DIR_PARAM = "queue_dir"
    PIPELINE_DTD_PARAM = "pipeline_dtd"
    RSAT_DIR_PARAM = "rsat_dir"
    RSAT_JASPAR_MOTIF_DATABASE_PARAM = "rsat_default_motif_database"
    MEME_DIR_PARAM = "meme_dir"
    
    CLUSTALW_COMMAND_PARAM = "clustalw_command"
    MAFFT_COMMAND_PARAM = "mafft_command"


    # File constants
    PROGRESSION_XSL_PATH = "resources/xsl/progression"
    PROGRESSION_XSL_FILE = "progression.xsl"
    PROGRESSION_XML_FILE = "progression.xml"

    # --------------------------------------------------------------------------------------
    # 'biology' constants
    DNA_ALPHABET =['A', 'C', 'G', 'T']
    MAX_INDEX = "Max"
    
    HG_BACKGROUND_MODEL = { 'A' : 0.325 , 'C' : 0.175, 'G' : 0.175, 'T' : 0.325}
    
    SEQUENCE_INIT_CHAR = "."
    SEQUENCE_INSERTION_CHAR = "-"
    
    NEGATIVE_STRAND = "-"
    POSITIVE_STRAND = "+"
    
    DIRECT_STRAND = "D"
    REVERSE_STRAND = "R"

    # --------------------------------------------------------------------------------------
    # various constants
    
    COMMENT_CHAR = "#"
    ORDERED = "Ordered"
    MIXED = "Mixed"
