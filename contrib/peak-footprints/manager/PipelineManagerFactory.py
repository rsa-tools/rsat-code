
from manager.PipelineManager import PipelineManager

class PipelineManagerFactory:
    
    managerInstance = None
    
    @staticmethod
    def getManager():
        
        if PipelineManagerFactory.managerInstance == None:
            PipelineManagerFactory.managerInstance = PipelineManager()
        
        return PipelineManagerFactory.managerInstance
        
    
            
    
