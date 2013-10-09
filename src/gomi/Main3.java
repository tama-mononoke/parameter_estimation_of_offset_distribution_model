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
public class Main3 {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try{
			PropertyConfigurator.configure( "log4j.properties" );
			vbNN();
			vbNDE();
			vbDDE();
			vbNN2();
			vbNDE2();
			vbDDE2();
			vbNN3();
			vbNDE3();
			vbDDE3();
		}
		catch(Exception ex){
			ex.printStackTrace();
		}

	}
	
	private static void vbNN() throws Exception{
		String args[] = {"onde","vb","test_result/bmlm2_1/sample_ONDE/ONN_prior.properties","test_result/bmlm2_1/result_ONDE/ONNdata.csv","test_result/bmlm2/result_ONDE/ONNVB2_result.properties"};
		BayesianUpdate.run(args);
	}
	
	private static void vbNDE() throws Exception{
		String args[] = {"onde","vb","test_result/bmlm2_1/sample_ONDE/ONDE_prior.properties","test_result/bmlm2_1/result_ONDE/ONDEdata.csv","test_result/bmlm2/result_ONDE/ONDEVB2_result.properties"};
		BayesianUpdate.run(args);
	}
	
	private static void vbDDE() throws Exception{
		String args[] = {"onde","vb","test_result/bmlm2_1/sample_ONDE/ODDE_prior.properties","test_result/bmlm2_1/result_ONDE/ODDEdata.csv","test_result/bmlm2/result_ONDE/ODDEVB2_result.properties"};
		BayesianUpdate.run(args);
	}
	
	private static void vbNN2() throws Exception{
		String args[] = {"onde","vb","test_result/bmlm2/sample_ONDE/ONN_prior.properties","test_result/bmlm2/result_ONDE/ONNdata.csv","test_result/bmlm2/result_ONDE/ONNVB3_result.properties"};
		BayesianUpdate.run(args);
	}
	
	private static void vbNDE2() throws Exception{
		String args[] = {"onde","vb","test_result/bmlm2/sample_ONDE/ONDE_prior.properties","test_result/bmlm2/result_ONDE/ONDEdata.csv","test_result/bmlm2/result_ONDE/ONDEVB3_result.properties"};
		BayesianUpdate.run(args);
	}
	
	private static void vbDDE2() throws Exception{
		String args[] = {"onde","vb","test_result/bmlm2/sample_ONDE/ODDE_prior.properties","test_result/bmlm2/result_ONDE/ODDEdata.csv","test_result/bmlm2/result_ONDE/ODDEVB3_result.properties"};
		BayesianUpdate.run(args);
	}
	
	private static void vbNN3() throws Exception{
		String args[] = {"onde","vb","test_result/bmlm2/sample_ONDE2/ONN_prior.properties","test_result/bmlm2/result_ONDE/ONNdata.csv","test_result/bmlm2/result_ONDE/ONNVB4_result.properties"};
		BayesianUpdate.run(args);
	}
	
	private static void vbNDE3() throws Exception{
		String args[] = {"onde","vb","test_result/bmlm2/sample_ONDE2/ONDE_prior.properties","test_result/bmlm2/result_ONDE/ONDEdata.csv","test_result/bmlm2/result_ONDE/ONDEVB4_result.properties"};
		BayesianUpdate.run(args);
	}
	
	private static void vbDDE3() throws Exception{
		String args[] = {"onde","vb","test_result/bmlm2/sample_ONDE2/ODDE_prior.properties","test_result/bmlm2/result_ONDE/ODDEdata.csv","test_result/bmlm2/result_ONDE/ODDEVB4_result.properties"};
		BayesianUpdate.run(args);
	}
	
}
