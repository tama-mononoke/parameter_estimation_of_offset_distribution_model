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
public class Main2 {
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try{
			PropertyConfigurator.configure( "log4j.properties" );
			gNN();
			gNDE();
			gDDE();
			gNN2();
			gNDE2();
			gDDE2();
			oedNN();
			oedNDE();
			oedDDE();
			oebNN();
			oebNDE();
			oebDDE();
			oedNN2();
			oedNDE2();
			oedDDE2();
			oebNN2();
			oebNDE2();
			oebDDE2();
			oebNN3();
			oebNDE3();
			oebDDE3();
		}
		catch(Exception ex){
			ex.printStackTrace();
		}
	}
	
	private static void gNN() throws Exception{
		String args[] = {"onde","rg","result_ONDE/ONN.properties","1000","result_ONDE_offset/ONNdata.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void gNDE() throws Exception{
		String args[] = {"onde","rg","result_ONDE/ONDE.properties","1000","result_ONDE_offset/ONDEdata.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void gDDE() throws Exception{
		String args[] = {"onde","rg","result_ONDE/ODDE.properties","1000","result_ONDE_offset/ODDEdata.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void gNN2() throws Exception{
		String args[] = {"onde","rg","result_ONDE/ONN2.properties","1000","result_ONDE_offset/ONN2data.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void gNDE2() throws Exception{
		String args[] = {"onde","rg","result_ONDE/ONDE2.properties","1000","result_ONDE_offset/ONDE2data.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void gDDE2() throws Exception{
		String args[] = {"onde","rg","result_ONDE/ODDE2.properties","1000","result_ONDE_offset/ODDE2data.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void oedNN() throws Exception{
		String args[] = {"onde","oe","d","result_ONDE/ONNEM_result.properties","result_ONDE_offset/ONNdata.csv","result_ONDE_offset/ONNresult1.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void oedNDE() throws Exception{
		String args[] = {"onde","oe","d","result_ONDE/ONDEEM_result.properties","result_ONDE_offset/ONDEdata.csv","result_ONDE_offset/ONDEresult1.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void oedDDE() throws Exception{
		String args[] = {"onde","oe","d","result_ONDE/ODDEEM_result.properties","result_ONDE_offset/ODDEdata.csv","result_ONDE_offset/ODDEresult1.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void oebNN() throws Exception{
		String args[] = {"onde","oe","b","result_ONDE/ONNVB2_result.properties","result_ONDE_offset/ONNresult1.csv","result_ONDE_offset/ONNresult2.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void oebNDE() throws Exception{
		String args[] = {"onde","oe","b","result_ONDE/ONDEVB2_result.properties","result_ONDE_offset/ONDEresult1.csv","result_ONDE_offset/ONDEresult2.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void oebDDE() throws Exception{
		String args[] = {"onde","oe","b","result_ONDE/ODDEVB2_result.properties","result_ONDE_offset/ODDEresult1.csv","result_ONDE_offset/ODDEresult2.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void oedNN2() throws Exception{
		String args[] = {"onde","oe","d","result_ONDE/ONNEM2_result.properties","result_ONDE_offset/ONN2data.csv","result_ONDE_offset/ONN2result1.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void oedNDE2() throws Exception{
		String args[] = {"onde","oe","d","result_ONDE/ONDEEM2_result.properties","result_ONDE_offset/ONDE2data.csv","result_ONDE_offset/ONDE2result1.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void oedDDE2() throws Exception{
		String args[] = {"onde","oe","d","result_ONDE/ODDEEM2_result.properties","result_ONDE_offset/ODDE2data.csv","result_ONDE_offset/ODDE2result1.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void oebNN2() throws Exception{
		String args[] = {"onde","oe","b","result_ONDE/ONNVB3_result.properties","result_ONDE_offset/ONN2result1.csv","result_ONDE_offset/ONN2result2.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void oebNDE2() throws Exception{
		String args[] = {"onde","oe","b","result_ONDE/ONDEVB3_result.properties","result_ONDE_offset/ONDE2result1.csv","result_ONDE_offset/ONDE2result2.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void oebDDE2() throws Exception{
		String args[] = {"onde","oe","b","result_ONDE/ODDEVB3_result.properties","result_ONDE_offset/ODDE2result1.csv","result_ONDE_offset/ODDE2result2.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void oebNN3() throws Exception{
		String args[] = {"onde","oe","b","result_ONDE/ONNVB4_result.properties","result_ONDE_offset/ONN2result2.csv","result_ONDE_offset/ONN2result3.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void oebNDE3() throws Exception{
		String args[] = {"onde","oe","b","result_ONDE/ONDEVB4_result.properties","result_ONDE_offset/ONDE2result2.csv","result_ONDE_offset/ONDE2result3.csv"};
		BayesianUpdate.run(args);
	}
	
	private static void oebDDE3() throws Exception{
		String args[] = {"onde","oe","b","result_ONDE/ODDEVB4_result.properties","result_ONDE_offset/ODDE2result2.csv","result_ONDE_offset/ODDE2result3.csv"};
		BayesianUpdate.run(args);
	}
}
