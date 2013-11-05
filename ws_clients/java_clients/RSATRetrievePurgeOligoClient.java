

import RSATWS.OligoAnalysisRequest;
import RSATWS.OligoAnalysisResponse;
import RSATWS.PurgeSequenceRequest;
import RSATWS.PurgeSequenceResponse;
import RSATWS.RSATWSPortType;
import RSATWS.RSATWebServicesLocator;
import RSATWS.RetrieveSequenceRequest;
import RSATWS.RetrieveSequenceResponse;


public class RSATRetrievePurgeOligoClient {

	/**
	 * This script runs a simple demo of the web service interface to the
	 * RSAT tools retrieve-seq, purge-sequence and oligo-analysis linked in a workflow.
	 * It sends a request to the server for discovering 6 letter words
	 * in upstream sequences of 5 yeast genes. The sequences are first
	 * retrieved and purged for repeated segments
	 */
	public static void main(String[] args) {
		 try
	        {
			 
			 System.out.println("This demo script illustrates a work flow combining three requests to the RSAT web services:\n\tretrieve-seq | purge-sequence | oligo-analysis");
			 
			 	String organism = "Saccharomyces_cerevisiae";
			 	
			 	/* Get the location of the service */
			 	RSATWebServicesLocator service = new RSATWebServicesLocator();
			 	RSATWSPortType proxy = service.getRSATWSPortType();

	
			 	/** Retrieve-seq part **/
			 	
	            /* prepare the parameters */
	            RetrieveSequenceRequest retrieveSeqParams = new RetrieveSequenceRequest();
	            	            
	            //Name of the query organism
	            retrieveSeqParams.setOrganism(organism);
	            //List of query genes
	            String[] q=  { "PHO5", "PHO8", "PHO11", "PHO81", "PHO84" };	            
	            retrieveSeqParams.setQuery(q);
	            // Clip sequences to avoid upstream ORFs
	            retrieveSeqParams.setNoorf(1);
	            retrieveSeqParams.setNocom(0);
	            // The result will stay in a file on the server
	            retrieveSeqParams.setOutput("server");
	            
	        	/* Call the service */
	            System.out.println("Retrieve-seq: sending request to the server...");
	            RetrieveSequenceResponse res = proxy.retrieve_seq(retrieveSeqParams);
	           
	            /* Process results  */  
	            //Report the remote command
	            System.out.println("Command used on the server:"+ res.getCommand());
	            //Report the server file location
	            String retrieveSeqFileServer = res.getServer();
	            System.out.println("Result file on the server::\n"+ res.getServer());
	            
	            /** Purge-sequence part **/
	            
	            /* prepare the parameters */
	            PurgeSequenceRequest purgeSeqParams = new PurgeSequenceRequest();
	            // The result will stay in a file on the server
	            purgeSeqParams.setOutput("server");
	            // Output from retrieve-seq part is used as input here
	            purgeSeqParams.setTmp_infile(retrieveSeqFileServer);
	            
	            /* Call the service */
	            System.out.println("Purge-sequence: sending request to the server...");
	            PurgeSequenceResponse res2 = proxy.purge_seq(purgeSeqParams);
	            
	            /* Process results  */  
	            //Report the remote command
	            System.out.println("Command used on the server:"+ res2.getCommand());
	            //Report the server file location
	            String purgeSeqFileServer = res2.getServer();
	            System.out.println("Result file on the server::\n"+ res2.getServer());
	            
	            /** Oligo-analysis part **/
	            
	            /* prepare the parameters */
	            OligoAnalysisRequest oligoParams = new OligoAnalysisRequest();
	            // Output from purge-seq part is used as input here
	            oligoParams.setTmp_infile(purgeSeqFileServer);
	            oligoParams.setOrganism(organism);
	            // Length of patterns to be discovered
	            oligoParams.setLength(6);
	            // Type of background used
	            oligoParams.setBackground("upstream-noorf");
	            // Returned statistics
	            oligoParams.setStats("occ,proba,rank");
	            // Do not allow overlapping patterns
	            oligoParams.setNoov(1);
	            // Search on both strands
	            oligoParams.setStr(2);
	            // Sort the result according to score
	            oligoParams.setSort(1);
	            // Lower limit to score is 0, less significant patterns are not displayed
	            String[] lth_values = {"occ_sig 0"};
	            oligoParams.setLth(lth_values);
	            
	            
	            /* Call the service */
	            System.out.println("Oligo-analysis: sending request to the server...");
	            OligoAnalysisResponse res3 = proxy.oligo_analysis(oligoParams);
	            
	            /* Process results  */  
	            //Report the remote command
	            System.out.println("Command used on the server:"+ res3.getCommand());
	            //Report the result
	            System.out.println("Discovered oligo(s):\n"+ res3.getClient());
	            //Report the server file location
	            System.out.println("Result file on the server::\n"+ res3.getServer());
	            
	        }
	        catch(Exception e) { System.out.println(e.toString()); 
	        }

	}

}
