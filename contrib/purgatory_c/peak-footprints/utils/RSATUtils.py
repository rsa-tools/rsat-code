
import os, commands

from log.Log import Log

from utils.FileUtils import FileUtils
from utils.exception.ExecutionException import ExecutionException

class RSATUtils:
    
    # Variables initialized by the PipelineMnaager at the tool start
    RSAT_PATH = ""
    RSAT_JASPAR_MOTIF_DATABASE = ""
    
    # ---------------------------------------------------------------------------------------------
    # Execute a chi2 test over the data listed in the two provided columns of the given histogram file
    @staticmethod
    def executeChi2Test( histogram_filepath, data_col, homogen_col):
        
        cmd = "grep -v '^#' " + histogram_filepath;
        cmd += " | cut -f " + str( data_col) + "," + str( homogen_col);
        cmd += " | " + os.path.join( RSATUtils.RSAT_PATH, "perl-scripts/transpose-table")
        cmd += " | cut -f 2-10000"
        cmd += " | " + os.path.join( RSATUtils.RSAT_PATH, "perl-scripts/chi-square -test goodness")
        
        cmd_result = commands.getstatusoutput( cmd)
        if cmd_result[0] != 0:
            Log.log( "RSATUtils.executeChi2Test : status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
            Log.log( " chi-square output is = \n" + str( cmd_result[1]))
            return None
            
        output = str( cmd_result[1])
        
        index_chi2 = output.find( '#chi2')
        if index_chi2 >= 0 :
            final_header = 'Sgroup'
            index_pval = output.find( final_header, index_chi2)
            if index_pval >= 0:
                results_string = output[ index_pval + len( final_header):]
                results = results_string.split()
                if len( results) >= 3:
                    chi2 = results[0]
                    pvalue = results[2]
                    return ( chi2, pvalue)
        
        Log.log( "RSATUtils.executeChi2Test : chi-square output does not contains desired results :" + str( output))
        
        return None
        
        
    
    # ---------------------------------------------------------------------------------------------
    # Permute the matrix in the given file using permute-matrix
    @staticmethod
    def permuteMatrix( matrix_file, file_format, destination_path, permutation_number):
    
        destination_files = []
    
        for index in range( permutation_number):
            
            # compose the destination file path
            destination_file = os.path.join( destination_path, os.path.basename( matrix_file) + "_" + str(index))
        
            # Compose the permute-matrix command
            cmd = os.path.join( RSATUtils.RSAT_PATH, "perl-scripts/permute-matrix")
            cmd += " -i '" + matrix_file + "'"
            cmd += " -in_format " + file_format
            cmd += " -o " + destination_file
            cmd += " -out_format " + file_format
            
            # Execute the command
            cmd_result = commands.getstatusoutput( cmd)
            if cmd_result[0] != 0:
                Log.log( "RSATUtils.permuteMatrix : status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
                Log.log( " permute-matrix output is = \n" + str( cmd_result[1]))
                return None
            
            destination_files.append( destination_file)
        
        return destination_files
            

    # ---------------------------------------------------------------------------------------------
    # Compute the histogram of chosen indexes and graph it in PNG
    @staticmethod
    def outputHistogram( table, histogram_interval, dir_path, prefix, title1, title2,  legendx,  legendy, other_columns, fhisto):
        return RSATUtils.outputHistogramFormat(table, histogram_interval, dir_path, prefix, title1, title2, legendx, legendy, other_columns, fhisto, "png")
    
    
    # ---------------------------------------------------------------------------------------------
    # Compute the histogram of chosen indexes and graph it
    @staticmethod
    def outputHistogramFormat( table, histogram_interval, dir_path, prefix, title1, title2,  legendx,  legendy, other_columns, fhisto, out_format):
        
        # save the stats to a tabbed file for classfreq command
        input_path = os.path.join( dir_path, prefix + ".tab")
        RSATUtils.outputTable( table, input_path)

        # execute the classfreq command
        histo_path = os.path.join( dir_path, prefix + "_histogram.tab")  
        
        cmd = os.path.join( RSATUtils.RSAT_PATH, "perl-scripts/classfreq")
        cmd += " -v 1"
        cmd += " -i '" + input_path + "'"
        cmd += " -col 1"
        cmd += " -ci " + str( histogram_interval)
        cmd += " -o '" + histo_path + "'"
        
        cmd_result = commands.getstatusoutput( cmd)
        if cmd_result[0] != 0:
            Log.log( "RSATUtils.outputHistogram : status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
            Log.log( " classfreq output is = \n" + str( cmd_result[1]))
            return


        # Build the graph corresponding to the histogram using RSAT XYGraph command
        graph_path = os.path.join( dir_path, prefix + "." + out_format)
        cmd = os.path.join( RSATUtils.RSAT_PATH, "perl-scripts/XYgraph")
        cmd += " -i '" + histo_path + "'"
        cmd += " -title1 '" + title1 + "'" 
        cmd += " -title2 '" + title2 + "'" 
        cmd += " -xcol 3 -ycol 4"
        if other_columns != None and len( other_columns) > 0:
            cmd += "," + ",".join( other_columns)
        cmd += " -xleg1 '" + legendx + "'"
        cmd += " -yleg1 '" + legendy + "'"
        cmd += " -legend -header -format "+ out_format
        if fhisto:
            cmd += " -fhisto"
        else:
            cmd += " -lines"
        cmd += " -o '" + graph_path + "'"
        
        cmd_result = commands.getstatusoutput( cmd)
        if cmd_result[0] != 0:
            Log.log( "RSATUtils.outputHistogram : status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
            Log.log( "  XYgraph output is = \n" + str( cmd_result[1]))

        return (histo_path, graph_path)



    

    # --------------------------------------------------------------------------------------
    # create the logo image of the motif which definition is provided in the given Transfac file
    @staticmethod
    def createLogoFromTF( input_path, input_file, motif_name):
        
        # /home/rsat/rsa-tools/perl-scripts/convert-matrix -logo_dir /home/rsat/rsa-tools/public_html/tmp  -v 1  -return logo -o /home/rsat/rsa-tools/public_html/tmp/convert-matrix_2011_09_29.124903_OtxJEc.res
        matrix_path = os.path.join( input_path, input_file)
        logo_path = os.path.join( input_path, motif_name)
        cmd = os.path.join( RSATUtils.RSAT_PATH, "perl-scripts/convert-matrix")
        cmd += " -i '" + matrix_path + "'"
        cmd += " -from transfac"
        cmd += " -logo_format png"
        cmd += " -logo_file '" + logo_path + "'"
        cmd += " -logo_opt '-e' -logo_opt '-M'"
        cmd += " -logo_dir " + input_path
        cmd += " -return logo"

        cmd_result = commands.getstatusoutput( cmd)
        if cmd_result[0] != 0:
            Log.log( "RSATUtils.createLogoFromTF : status returned is :" + str( cmd_result[0]) + " for command '" + cmd + "'" )
            Log.log( "  convert-matrix output is = \n" + str( cmd_result[1]))
        
        # Try to remove the reverse complement logo 
        rc_logo_path = os.path.join( input_path, motif_name + "_m1_rc.png")
        if os.path.exists( rc_logo_path):
            FileUtils.removeFile( rc_logo_path)
        
        

    # --------------------------------------------------------------------------------------
    # output the given table to file
    @staticmethod
    def outputTable( table, path):

        try:
            out_file = FileUtils.openFile( path, "w")
            for number in table:
                out_file.write( str( number) + "\n")
                out_file.flush()
            out_file.close()
        except IOError, io_exce:
            raise ExecutionException( "HistogramProcessor.outputMotifStatistics : Unable to build statistics out_file. From:\n\t---> " +str( io_exce))
