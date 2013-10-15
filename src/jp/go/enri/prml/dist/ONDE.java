/**
 * 
 */
package jp.go.enri.prml.dist;

import java.util.Arrays;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * Parameter of offset mixture distribution with Gaussian and Laplace components
 * @author Masato Fujita（Electronic Navigation Research Institute）
 * @version 1.0.1　(Last update: 06/12/2011)
 *
 */
public class ONDE {
	/**
	 * Get the feasible offset values.
	 * @return the feasible offset values
	 */
	public double[] getOffset() {
		return offset;
	}
	/**
	 * Get the offset mixing coefficients.
	 * @return the offset mixing coefficients.
	 */
	public double[] getOmega() {
		return omega;
	}

	/**
	 * Log
	 */
	public static Log log = LogFactory.getLog(ONDE.class);
	/**
	 * Get the number of Gaussian components in the mixture distribution..
	 * @return the number of Gaussian components in the mixture distribution.
	 */
	public int getM() {
		return m;
	}
	/**
	 * Get the number of Laplace components in the mixture distribution.
	 * @return the number of Laplace components in the mixture distribution.
	 */
	public int getN() {
		return n;
	}
	/**
	 * Get the mixing coefficients.
	 * @return the mixing coefficients.
	 */
	public double[] getPi() {
		return Arrays.copyOf(pi, pi.length);
	}
	/**
	 * Get the scale parameter of Laplace distributions.
	 * @return the scale parameter of Laplace distributions
	 */
	public double[] getLambda() {
		return Arrays.copyOf(lambda, lambda.length);
	}
	/**
	 * Get the standard deviation of Gaussian distributions.
	 * @return sigma the standard deviation of Gaussian distributions
	 */
	public double[] getSigma() {
		return Arrays.copyOf(sigma, sigma.length);
	}
	/**
	 * Get the number of the feasible offset values.
	 * @return the number of the feasible offset values
	 */
	public int getL(){
		return offset.length;
	}
	/**
	 * the number of Gaussian components in the mixture distribution.
	 */
	int m;
	/**
	 * the number of Laplace components in the mixture distribution.
	 */
	int n;
	/**
	 * the mixing coefficients
	 */
	double pi[];
	/**
	 * the scale parameter of Laplace distributions
	 */
	double lambda[];
	/**
	 * the standard deviation of Gaussian distributions
	 */
	double sigma[];
	/**
	 * the feasible offset values
	 */
	double offset[];
	/**
	 * the offset mixing coefficients
	 */
	double omega[];
	/**
	 * Constructor
	 */
	private ONDE(){}
	/**
	 * Constructor
	 * @param pi the mixing coefficients
	 * @param sigma the standard deviation of Gaussian distributions
	 * @param lambda the scale parameter of Laplace distributions
	 * @param offset the feasible offset values
	 * @param omega the offset mixing coefficients
	 */
	public ONDE(double[] pi, double[] sigma, double[] lambda, double[] offset, double[] omega) {
		super();
		if(sigma==null){
			log.error("The length of array sigma is not correct.");
			throw new IllegalArgumentException();
		}
		this.m = sigma.length;
		if(lambda==null){
			log.error("The length of array lambda is not correct.");
			throw new IllegalArgumentException();
		}
		this.n = lambda.length;
		int K = m+n;
		if(pi==null || pi.length!=K){
			log.error("The length of array pi is not correct.");
			throw new IllegalArgumentException();
		}
		for(int i=0;i<pi.length;i++){
			if(pi[i]<0){
				log.error("The value pi[" + i + "] should not be negative. pi[" + i + "]=" + pi[i]);
				throw new IllegalArgumentException();
			}
		}
		this.pi = Arrays.copyOf(pi, pi.length);
		for(int i=0;i<sigma.length;i++){
			if(sigma[i]<=0){
				log.error("The value sigma[" + i + "] should be positive. sigma[" + i + "]=" + sigma[i]);
				throw new IllegalArgumentException();
			}
		}
		this.sigma = Arrays.copyOf(sigma, sigma.length);
		for(int i=0;i<lambda.length;i++){
			if(lambda[i]<=0){
				log.error("The value lambda[" + i + "] should be positive. lambda[" + i + "]=" + lambda[i]);
				throw new IllegalArgumentException();
			}
		}
		this.lambda = Arrays.copyOf(lambda, lambda.length);
		if(offset==null || offset.length<1){
			log.error("The length of array offfset is not correct.");
			throw new IllegalArgumentException();
		}
		this.offset = Arrays.copyOf(offset,offset.length);
		if(omega==null || omega.length!=offset.length){
			log.error("The length of array omega is not correct.");
			throw new IllegalArgumentException();
		}
		this.omega = Arrays.copyOf(omega, omega.length);
	}
	
	/**
	 * @see java.lang.Object#clone()
	 */
	@Override
	public ONDE clone(){
		ONDE tmp = new ONDE();
		tmp.m = m;
		tmp.n = n;
		tmp.pi = Arrays.copyOf(pi, pi.length);
		tmp.sigma = Arrays.copyOf(sigma, sigma.length);
		tmp.lambda = Arrays.copyOf(lambda, lambda.length);
		tmp.offset = Arrays.copyOf(offset,offset.length);
		tmp.omega = Arrays.copyOf(omega, omega.length);
		return tmp;
	}
	/**
	 * Estimate the probability that the given offset value is applied.
	 * @param x Observation. 
	 * @return the probability that each offset value is applied.
	 */
	public double[] estimateOffsetProbability(double x){
		double ans[] = new double[offset.length];
		double sum = 0;
		for(int l=0;l<offset.length;l++){
			double tmp2 = 0;
			for(int k=0;k<m;k++){
				tmp2 += pi[k]*Math.exp(-0.5*Math.pow(x-offset[l], 2)/(sigma[k]*sigma[k]))/(Math.sqrt(2*Math.PI)*sigma[k]);
			}
			for(int k=0;k<n;k++){
				tmp2 += pi[k+m]*Math.exp(-Math.abs(x-offset[l])/lambda[k])/(2*lambda[k]);
			}
			tmp2 *= omega[l];
			ans[l] = tmp2;
			sum += tmp2;
		}
		for(int l=0;l<offset.length;l++) ans[l] /= sum;
		return ans;
	}
	
}
