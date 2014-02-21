import os, commands
import rpy2.robjects as robjects

from rpy2.robjects.packages import importr

from log.Log import Log

from utils.FileUtils import FileUtils
from utils.exception.ExecutionException import ExecutionException

class RUtils:

    # ---------------------------------------------------------------------------------------------
    # Compute the histogram of chosen indexes and graph it with R
    @staticmethod
    def plotHistogram( table, histogram_interval, dir_path, prefix, title1, title2,  legendx, legendy, out_format):
        
        graph_path = os.path.join( dir_path, prefix + "." + out_format)
        
        graphics = importr('graphics')
        grdevices = importr('grDevices')
            
        r = robjects.r
        
        # Open the output file in the right format
        if out_format == "png" :
            grdevices.png(file=graph_path, width=512, height=512)
        elif out_format == "ps" :
            grdevices.postscript(file=graph_path, width=512, height=512)
        elif out_format == "ps" :
            grdevices.pdf(file=graph_path, width=512, height=512)
            
        # Get the values to plot
        table_x = robjects.FloatVector( table)

        # Prepare the class limits of the histogram        
        p_breaks1 = range(0, max( table)+histogram_interval, histogram_interval)
        p_breaks2 = range(-histogram_interval, min( table)-histogram_interval, -histogram_interval)
        p_breaks2.reverse()
        p_breaks = p_breaks2 + p_breaks1
        table_breaks = robjects.IntVector( p_breaks)
        
        # Plot the histogram
        graphics.hist( table_x, breaks=table_breaks, col = "blue", main = title1 + "\n" + title2, xlab = legendx, ylab = legendy)
        
        # Close the file
        grdevices.dev_off()

        histo_path = ""
        
        return (histo_path, graph_path)
    
        # ---------------------------------------------------------------------------------------------
    # Compute the histogram of chosen indexes and graph it with R
    @staticmethod
    def plotTwoHistograms( table1, table2, histogram_interval, dir_path, prefix, title1, title2,  legendx, legendy, out_format):
        
        graph_path = os.path.join( dir_path, prefix + "." + out_format)
        
        graphics = importr('graphics')
        grdevices = importr('grDevices')
        stats = importr('stats')
            
        r = robjects.r
        
        # Open the output file in the right format
        if out_format == "png" :
            grdevices.png(file=graph_path, width=512, height=512)
        elif out_format == "pdf" :
            grdevices.pdf(file=graph_path, width=512, height=512)
       
        # Get the values to plot
        table_x1 = robjects.FloatVector( table1)
        table_x2 = robjects.FloatVector( table2)

        # Prepare the class limits of the histogram        
        p_breaks1 = range(0, max( table)+histogram_interval, histogram_interval)
        p_breaks2 = range(-histogram_interval, min( table)-histogram_interval, -histogram_interval)
        p_breaks2.reverse()
        p_breaks = p_breaks2 + p_breaks1
        table_breaks = robjects.IntVector( p_breaks)
        
        # Plot the histograms
        graphics.hist( table_x, breaks=table_breaks, col = "blue", main = title1 + "\n" + title2, xlab = legendx, ylab = legendy)
        graphics.lines(stats.density( table_x2), col = "green")
        
        # Close the file
        grdevices.dev_off()

        histo_path = ""
        
        return (histo_path, graph_path)