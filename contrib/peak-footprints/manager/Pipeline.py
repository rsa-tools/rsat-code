
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
    # Return the component which has the given name as processor name. Return None if no
    # suitable component is found
    def getComponent(self, component_name):
        
        for component in self.componentList:
            if component.processorName == component_name:
                return component
        
        return None

    # --------------------------------------------------------------------------------------
    # Return the component list of the pipeline
    def getComponentList(self):
        
        return self.componentList

    # --------------------------------------------------------------------------------------
    # Return a string representation of the pipeline
    def toString( self):
        
        result = "Pipeline '" + self.name + "'\n"
        for component in self.componentList:
            result += "|--" + component.toString()
            
        return result

