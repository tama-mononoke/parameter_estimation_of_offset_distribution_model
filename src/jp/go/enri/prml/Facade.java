/**
 * 
 */
package jp.go.enri.prml;

import java.util.Arrays;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.log4j.PropertyConfigurator;

/**
 * Run the software.
 * @author Masatoo Fujita
 *
 */
public class Facade {
	/**
	 * Log
	 */
	public static Log log = LogFactory.getLog(Facade.class);
	/**
	 * run
	 * @param args
	 * arg[0]: "u" Bayesian update programs
	 */
	public static void main(String[] args) {
		PropertyConfigurator.configure( "log4j.properties" );
		try{
			String tmp[] = Arrays.copyOfRange(args, 1, args.length);
			if(args[0].equals("u")) BayesianUpdate.run(tmp);
			else throw new IllegalArgumentException();
		}
		catch(Exception ex){
			log.error("Error!", ex);
			ex.printStackTrace();
		}
	}

}
