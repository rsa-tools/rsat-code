import datetime

def time_warn (message, level=1):
    """This function prints a message to the standard error, and prepends
    the date and time in a format ensuring alphanumerical sorting.
    
    The message can be a string or a list of strings, which will
    be concatenated with tabulations as separators.
    
    """	
    if not "verbosity" in globals():
        verbosity = 1

    # Return if main verbosity smaller than message level
    if (verbosity < level):
        return()

    # Get current date and time
    alphadate = str(datetime.datetime.now())
    
    # Print date + time + message
    if isinstance(message, str):
        sys.stderr.write("%s\t%s\n"% (alphadate,message))
    else:
        sys.stderr.write("%s\t%s\n"% (alphadate,"\t".join(message)))
