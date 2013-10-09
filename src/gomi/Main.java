/**
 * 
 */
package gomi;

import jp.go.enri.prml.BayesianUpdate;

import org.apache.log4j.PropertyConfigurator;

/**
 * @author m-fujita
 *
 */
public class Main {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try{
			PropertyConfigurator.configure( "log4j.properties" );
			gNN();
			gNDE();
			gDDE();
			emNN();
			emNDE();
			emDDE();
			vbNN();
			vbNDE();
			vbDDE();
		}
		catch(Exception ex){
			ex.printStackTrace();
		}

	}
	
	private static void gNN() throws Exception{
		String args[] = {"onde","rg","sample_ONDE/ONN.properties","10000","result_ONDE/ONNdata.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void gNDE() throws Exception{
		String args[] = {"onde","rg","sample_ONDE/ONDE.properties","10000","result_ONDE/ONDEdata.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void gDDE() throws Exception{
		String args[] = {"onde","rg","sample_ONDE/ODDE.properties","10000","result_ONDE/ODDEdata.csv"};
		BayesianUpdate.run(args);
	}

	private static void emNN() throws Exception{
		String args[] = {"onde","em","sample_ONDE/ONNEM.properties","result_ONDE/ONNdata.csv","result_ONDE/ONNEM_result.properties"};
		BayesianUpdate.run(args);
	}
	
	private static void emNDE() throws Exception{
		String args[] = {"onde","em","sample_ONDE/ONDEEM.properties","result_ONDE/ONDEdata.csv","result_ONDE/ONDEEM_result.properties"};
		BayesianUpdate.run(args);
	}
	
	private static void emDDE() throws Exception{
		String args[] = {"onde","em","sample_ONDE/ODDEEM.properties","result_ONDE/ODDEdata.csv","result_ONDE/ODDEEM_result.properties"};
		BayesianUpdate.run(args);
	}
	
	private static void vbNN() throws Exception{
		String args[] = {"onde","vb","sample_ONDE/ONN_prior.properties","result_ONDE/ONNdata.csv","result_ONDE/ONNVB_result.properties"};
		BayesianUpdate.run(args);
	}
	
	private static void vbNDE() throws Exception{
		String args[] = {"onde","vb","sample_ONDE/ONDE_prior.properties","result_ONDE/ONDEdata.csv","result_ONDE/ONDEVB_result.properties"};
		BayesianUpdate.run(args);
	}
	
	private static void vbDDE() throws Exception{
		String args[] = {"onde","vb","sample_ONDE/ODDE_prior.properties","result_ONDE/ODDEdata.csv","result_ONDE/ODDEVB_result.properties"};
		BayesianUpdate.run(args);
	}
	
}
