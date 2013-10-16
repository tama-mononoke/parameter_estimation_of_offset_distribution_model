/**
 * 
 */
package jp.go.enri.prml.vb;

import java.util.Arrays;

import jp.go.enri.prml.dist.Dirichlet;
import jp.go.enri.prml.dist.Gamma;
import jp.go.enri.prml.dist.NDE;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * Prior/posterior distribution parameters of the mixture of Gaussian and Laplace distributions
 * @author Masato Fujita (Electronic Navigation Research Institute)
 * @version 1.0.1 (Last update: 30/11/2011)
 */
public class NDEParameterDistribution {
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
	 * Get the parameters of Dirichlet prior/posterior distribution of the mixing coefficients \pi. 
	 * @return the parameters of Dirichlet prior/posterior distribution of the mixing coefficients \pi 
	 */
	public double[] getAlpha() {
		return Arrays.copyOf(alpha, alpha.length);
	}

	/**
	 * Get the parameters `a' of Gamma prior/posterior distribution of the accuracy parameters \eta.
	 * @return the parameters `a' of Gamma prior/posterior distribution of the accuracy parameters \eta 
	 */
	public double[] getA() {
		return Arrays.copyOf(a, a.length);
	}

	/**
	 * Get the parameters `b' of Gamma prior/posterior distribution of the accuracy parameters \eta.
	 * @return the parameters `b' of Gamma prior/posterior distribution of the accuracy parameters \eta 
	 */
	public double[] getB() {
		return Arrays.copyOf(b, b.length);
	}

	/**
	 * Log
	 */
	public static Log log = LogFactory.getLog(NDEParameterDistribution.class);
	/**
	 * the number of Gaussian components in the mixture distribution.
	 */
	int m;
	/**
	 * the number of Laplace components in the mixture distribution.
	 */
	int n;
	/**
	 * the parameters of Dirichlet prior/posterior distribution of the mixing coefficients \pi 
	 */
	double alpha[];
	/**
	 * the parameters `a' of Gamma prior/posterior distribution of the accuracy parameters \eta 
	 */
	double a[];
	/**
	 * the parameters `b' of Gamma prior/posterior distribution of the accuracy parameters \eta 
	 */
	double b[];
	/**
	 * Gamma distribution for random sample generation
	 */
	private Gamma gamma[];
	/**
	 * Dirichlet distribution for random sample generation
	 */
	private Dirichlet dirichlet;
	/**
	 * @see java.lang.Object#clone()
	 */
	@Override
	public NDEParameterDistribution clone(){
		double[] alpha_tmp = Arrays.copyOf(alpha, alpha.length);
		double[] a_tmp = Arrays.copyOf(a, a.length);
		double[] b_tmp = Arrays.copyOf(b, b.length);
		return new NDEParameterDistribution(m, n, alpha_tmp, a_tmp, b_tmp);
	}
	
	/**
	 * Constructor
	 * @param m the number of Gaussian components in the mixture
	 * @param n the number of Laplace components in the mixture
	 * @param alpha the parameters of Dirichlet prior/posterior distribution of the mixing coefficients \pi ((m+n)-dimensional array)
	 * @param a the parameters `a' of Gamma prior/posterior distribution of the accuracy parameters \eta ((m+n)-dimensional array)
	 * @param b the parameters `b' of Gamma prior/posterior distribution of the accuracy parameters \eta  ((m+n)-dimensional array)
	 */
	public NDEParameterDistribution(int m, int n, double[] alpha, double[] a, double[] b){
		initialize(m, n, alpha, a, b);
	}
	
	/**
	 * Constructor
	 * @param m the number of Gaussian components in the mixture
	 * @param n the number of Laplace components in the mixture
	 * @param alpha0 the parameters of Dirichlet prior/posterior distribution of the mixing coefficients \pi (Common)
	 * @param a the parameters `a' of Gamma prior/posterior distribution of the accuracy parameters \eta ((m+n)-dimensional array)
	 * @param b the parameters `b' of Gamma prior/posterior distribution of the accuracy parameters \eta  ((m+n)-dimensional array)
	 */
	public NDEParameterDistribution(int m, int n, double alpha0, double[] a, double[] b){
		double tmpAlpha[] = new double[m+n];
		Arrays.fill(tmpAlpha, alpha0);
		initialize(m, n, tmpAlpha, a, b);
	}
	

	/**
	 * Initialization
	 * @param m the number of Gaussian components in the mixture
	 * @param n the number of Laplace components in the mixture
	 * @param alpha the parameters of Dirichlet prior/posterior distribution of the mixing coefficients \pi ((m+n)-dimensional array)
	 * @param a the parameters `a' of Gamma prior/posterior distribution of the accuracy parameters \eta ((m+n)-dimensional array)
	 * @param b the parameters `b' of Gamma prior/posterior distribution of the accuracy parameters \eta  ((m+n)-dimensional array)
	 */
	private void initialize(int m, int n, double[] alpha, double[] a, double[] b){
		if(m<0 || n<0){
			log.error("The number m and n should not be negative.");
			throw new IllegalArgumentException();
		}
		this.m = m;
		this.n = n;
		int K = m+n;
		if(alpha==null || alpha.length != K || a==null || a.length != K || b==null || b.length != K){
			log.error("The length of arrays for is illegal.");
			throw new IllegalArgumentException();
		}
		for(int i=0;i<K;i++){
			if(alpha[i]<=0){
				log.error("alpha[" + i + "] should be positive.");
				throw new IllegalArgumentException();
			}
			if(a[i]<=0){
				log.error("a[" + i + "] should be positive.");
				throw new IllegalArgumentException();
			}
			if(b[i]<=0){
				log.error("b[" + i + "] should be positive.");
				throw new IllegalArgumentException();
			}
		}
		this.alpha = Arrays.copyOf(alpha, alpha.length);
		this.a = Arrays.copyOf(a, a.length);
		this.b = Arrays.copyOf(b, b.length);
		// Update the other parameters.
		updateOtherParameters();
	}
	/**
	 * Update the other parameters.
	 */
	void updateOtherParameters(){
		gamma = new Gamma[m+n];
		for(int i=0;i<m+n;i++) gamma[i] = new Gamma(a[i], b[i]);
		dirichlet = new Dirichlet(alpha);
	}
	
	/**
	 * MAP estimation of parameter values.
	 * @return MAP estimation of parameter values.
	 */
	public NDE MAPEstimation(){
		log.info("Variational Bayesian MAP estimation process started.");
		int K = m+n;
		double pi[] = new double[K];
		double sum = 0;
		for(int i=0;i<K;i++) sum += alpha[i];
		for(int i=0;i<K;i++) pi[i] = (alpha[i]-1)/(sum-K);
		double eta[] = new double[K];
		for(int i=0;i<K;i++) eta[i] = (a[i]-1)/b[i];
		double sigma[] = new double[m];
		for(int i=0;i<m;i++){
			sigma[i] = Math.sqrt(1/eta[i]);
		}
		double lambda[] = new double[n];
		for(int i=0;i<n;i++){
			lambda[i] = 1/eta[i+m];
		}
		NDE para = new NDE(pi,sigma,lambda);
		log.info("Variational Bayesian MAP estimation process finished.");
		if(log.isInfoEnabled()){
			StringBuilder sb = new StringBuilder();
			for(int i=0;i<K;i++){
				sb.append((i+1));
				sb.append(":");
				sb.append("pi=");
				sb.append(pi[i]);
				if(i<m){
					sb.append(", sigma=");
					sb.append(sigma[i]);
				}
				else{
					sb.append(", lambda=");
					sb.append(lambda[i-m]);
				}
				sb.append(";");
			}
			log.info(sb.toString());
		}
		return para;
	}
	
	/**
	 * Average estimation of parameter values.
	 * @return Average estimation of parameter values
	 */
	public NDE ptEstimation(){
		log.info("Variational Bayesian point estimation process started.");
		int K = m+n;
		double pi[] = new double[K];
		double sum = 0;
		for(int i=0;i<K;i++) sum += alpha[i];
		for(int i=0;i<K;i++) pi[i] = alpha[i]/sum;
		double eta[] = new double[K];
		for(int i=0;i<K;i++) eta[i] = a[i]/b[i];
		double sigma[] = new double[m];
		for(int i=0;i<m;i++){
			sigma[i] = Math.sqrt(1/eta[i]);
		}
		double lambda[] = new double[n];
		for(int i=0;i<n;i++){
			lambda[i] = 1/eta[i+m];
		}
		NDE para = new NDE(pi,sigma,lambda);
		log.info("Variational Bayesian point estimation process finished.");
		if(log.isInfoEnabled()){
			StringBuilder sb = new StringBuilder();
			for(int i=0;i<K;i++){
				sb.append((i+1));
				sb.append(":");
				sb.append("pi=");
				sb.append(pi[i]);
				if(i<m){
					sb.append(", sigma=");
					sb.append(sigma[i]);
				}
				else{
					sb.append(", lambda=");
					sb.append(lambda[i-m]);
				}
				sb.append(";");
			}
			log.info(sb.toString());
		}
		return para;
	}
	
	/**
	 * Sampling from the parameter distribution.
	 * @return sample
	 */
	public NDE sampling(){
		double sigma[] = new double[m];
		double lambda[] = new double[n];
		for(int i=0;i<m;i++) sigma[i] = 1/Math.sqrt(gamma[i].nextDouble());
		for(int i=0;i<n;i++) lambda[i] = 1/gamma[i+m].nextDouble();
		return new NDE(dirichlet.nextDouble(), sigma, lambda);
	}
	
}
