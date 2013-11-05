
import RSATWS.RSATWSPortType;
import RSATWS.RSATWebServicesLocator;
import RSATWS.RandomSequenceRequest;
import RSATWS.RandomSequenceResponse;


public class RSATRandomSeqClient {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		 try
	        {
			 	/* Get the location of the service */
			 	RSATWebServicesLocator service = new RSATWebServicesLocator();
			 	RSATWSPortType proxy = service.getRSATWSPortType();

	         	            
	            /* prepare the parameters */
	            RandomSequenceRequest parameters = new RandomSequenceRequest();
	            
	            //Name of the query organism
	            parameters.setOrganism("Mus_musculus_EnsEMBL");
	            //Number of sequences to generate
	            parameters.setRepetition(5);
	            //Length of sequences to generate
	            parameters.setSequence_length(500);
	            //Type of sequence
	            parameters.setType("dna");
	            
	            // Precalculated background model : 
	            // type of bg model
	            parameters.setBg_model("upstream");
	            // length of oligomer (=markov order + 1)
	            parameters.setOligo_length(3);
	            
	            
          
	            
	        	/* Call the service */
	            System.out.println("Calling RSAT server...");
	            RandomSequenceResponse res = proxy.random_seq(parameters);
	           


	            /* Process results  */  
	            //Report the remote command
	            System.out.println("Command used on the server:"+ res.getCommand());
	            //Report the result
	            System.out.println("Result:\n"+ res.getClient());
	        }
	        catch(Exception e) { System.out.println(e.toString()); 
	        }

	}

}

