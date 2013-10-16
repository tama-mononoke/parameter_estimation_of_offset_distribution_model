/**
 * 
 */
package jp.go.enri.prml.bf;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * Parameters used for brute-force search parameter estimation (mixture distribution)
 * @author Masato Fujita
 *
 */
public class NDEParameter {
	/**
	 * Log
	 */
	public static Log log = LogFactory.getLog(NDEParameter.class);
	/**
	 * Division number
	 */
	int D;
	/**
	 * Possible range of sigma
	 */
	double initial_sigma_bound[][];
	/**
	 * Possible range of lambda
	 */
	double initial_lambda_bound[][];
	/**
	 * Constructor
	 * @param d division number 
	 * @param initial_sigma_bound Possible range of sigma
	 * @param initial_lambda_bound Possible range of lambda
	 */
	public NDEParameter(int d, double[][] initial_sigma_bound,
			double[][] initial_lambda_bound) {
		super();
		D = d;
		if(D<=0 || initial_lambda_bound==null || initial_sigma_bound==null){
			log.error("invalid argument");
			throw new IllegalArgumentException();
		}
		this.initial_sigma_bound = initial_sigma_bound;
		this.initial_lambda_bound = initial_lambda_bound;
	}
	
	
	
}
