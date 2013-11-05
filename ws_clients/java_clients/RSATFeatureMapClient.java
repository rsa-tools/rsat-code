
import RSATWS.RSATWSPortType;
import RSATWS.RSATWebServicesLocator;
import RSATWS.FeatureMapRequest;
import RSATWS.FeatureMapResponse;


public class RSATFeatureMapClient {

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
	            FeatureMapRequest parameters = new FeatureMapRequest();

		    // Output
		    parameters.setOutput("server");
	            
	            // Features
	            parameters.setFeatures("PHO5    dnapat  acgtgc|gcacgt   R       438     443     ACGTGC  4.37\nPHO8    dnapat  acgtgc|gcacgt   D       268     273     ACGTGC  4.37");

	            // Legend
	            parameters.setLegend(1);
	            // Scalebar
	            parameters.setScalebar(1);
	            // Scalestep
	            parameters.setScalestep(50);
	            // Scorethick
	            parameters.setScorethick(1);
	            // Picture format
	            parameters.setFormat("jpg");            
	            
	        	/* Call the service */
	            System.out.println("Calling RSAT server...");
	            FeatureMapResponse res = proxy.feature_map(parameters);

	            /* Process results  */  
	            //Report the remote command
	            System.out.println("Command used on the server:"+ res.getCommand());
	            //Report the result
	            System.out.println("Result:\n"+ res.getServer());
	        }
	        catch(Exception e) { System.out.println(e.toString()); 
	        }

	}

}

