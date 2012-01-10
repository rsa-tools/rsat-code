
class Pipeline:

    # --------------------------------------------------------------------------------------
    def __init__( self):
        
        self.name = ""
        self.componentList = []
        self.firstComponents = []
        self.lastComponents = []
    

    # --------------------------------------------------------------------------------------
    # Add the given component to the pipeline
    def addComponent(self, component):
        
        if not component in self.componentList:
            self.componentList.append( component)


    # --------------------------------------------------------------------------------------
    # Add the given component to the list of first executed components of the pipeline
    def addFirstComponent(self, component):

        if not component in self.firstComponents:
            self.firstComponents.append( component)


    # --------------------------------------------------------------------------------------
    # Add the given component to the list of last executed components of the pipeline
    def addLastComponent(self, component):

        if not component in self.lastComponents:
            self.lastComponents.append( component)


    # --------------------------------------------------------------------------------------
    # Return a string representation of the pipeline
    def toString( self):
        
        result = "Pipeline '" + self.name + "'\n"
        for component in self.componentList:
            result += "|--" + component.toString()
            
        return result
