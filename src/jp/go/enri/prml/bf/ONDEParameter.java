/**
 * 
 */
package jp.go.enri.prml.bf;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * Parameters used for brute-force search parameter estimation (offset mixture distribution)
 * @author Masato Fujita
 *
 */
public class ONDEParameter {
	/**
	 * Log
	 */
	public static Log log = LogFactory.getLog(ONDEParameter.class);
	/**
	 * Division number
	 */
	int D;
	/**
	 * Feasible offset values
	 */
	double offset[];
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
	 * @param d Division number
	 * @param offset Feasible offset values
	 * @param initial_sigma_bound Possible range of sigma
	 * @param initial_lambda_bound Possible range of lambda
	 */
	public ONDEParameter(int d, double offset[], double[][] initial_sigma_bound,
			double[][] initial_lambda_bound) {
		super();
		D = d;
		if(D<=0 || initial_lambda_bound==null || initial_sigma_bound==null){
			log.error("invalid argument");
			throw new IllegalArgumentException();
		}
		this.offset = offset;
		this.initial_sigma_bound = initial_sigma_bound;
		this.initial_lambda_bound = initial_lambda_bound;
	}
	
	
	
}
