import RSATWS.RSATWSPortTypeProxy;
import RSATWS.SupportedOrganismsRequest;
import RSATWS.SupportedOrganismsResponse;


public class RSATSupportedOrganismClient {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		 try
	        {
			 	/* Get the location of the service */
			 	RSATWSPortTypeProxy proxy = new RSATWSPortTypeProxy();

	         	            
	            /* prepare the parameters */
			 	SupportedOrganismsRequest parameters = new SupportedOrganismsRequest();

          
	            
	        	/* Call the service */
	            System.out.println("Calling RSAT server...");
	            SupportedOrganismsResponse res = proxy.supported_organisms(parameters);
	           


	            /* Process results  */  
	            //Report the remote command
	            System.out.println("Command used on the server:"+ res.getCommand());
	            //Report the result
	            System.out.println("Gene(s) info(s):\n"+ res.getClient());
	        }
	        catch(Exception e) { System.out.println(e.toString()); 
	        }

	}

}
