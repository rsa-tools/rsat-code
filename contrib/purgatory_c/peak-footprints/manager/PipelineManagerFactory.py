
from manager.PipelineManager import PipelineManager

class PipelineManagerFactory:
    
    managerInstance = None
    
    @staticmethod
    def getManager( output_dir, rsat_path):
        
        if PipelineManagerFactory.managerInstance == None:
            PipelineManagerFactory.managerInstance = PipelineManager( output_dir, rsat_path)
        
        return PipelineManagerFactory.managerInstance
        
    
            
    
