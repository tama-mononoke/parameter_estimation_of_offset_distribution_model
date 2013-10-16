/**
 * 
 */
package jp.go.enri.prml.vb;

import java.util.Arrays;

import jp.go.enri.prml.dist.Dirichlet;
import jp.go.enri.prml.dist.Gamma;
import jp.go.enri.prml.dist.ONDE;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * Prior/posterior distribution parameters of the offset mixture of Gaussian and Laplace distributions
 * @author Masato Fujita (Electronic Navigation Research Institute)
 * @version 1.0.1 (Last update: 06/12/2011)
 */
public class ONDEParameterDistribution {
	/**
	 * Get the parameters of Dirichlet prior/posterior distribution of the offset mixing coefficients \omega. 
	 * @return the parameters of Dirichlet prior/posterior distribution of the offset mixing coefficients \omega. 
	 */
	public double[] getP() {
		return p;
	}

	/**
	 * Get the feasible offset values. 
	 * @return the feasible offset values
	 */
	public double[] getOffset() {
		return offset;
	}

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
	public static Log log = LogFactory.getLog(ONDEParameterDistribution.class);
	/**
	 * the number of Gaussian components in the mixture distribution.
	 */
	int m;
	/**
	 * the number of Laplace components in the mixture distribution.
	 */
	int n;
	/**
	 * the parameters of Dirichlet prior/posterior distribution of the offset mixing coefficients \omega
	 */
	double p[];
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
	 * the feasible offset values
	 */
	double[] offset;
	/**
	 * Gamma distribution for random sample generation
	 */
	private Gamma gamma[];
	/**
	 * Dirichlet distribution for random sample generation
	 */
	private Dirichlet dirichlet_alpha;
	/**
	 * Dirichlet distribution for random sample generation
	 */
	private Dirichlet dirichlet_p;
	
	/**
	 * @see java.lang.Object#clone()
	 */
	@Override
	public ONDEParameterDistribution clone(){
		double[] alpha_tmp = Arrays.copyOf(alpha, alpha.length);
		double[] off_tmp = Arrays.copyOf(offset, offset.length);
		double[] p_tmp = Arrays.copyOf(p, p.length);
		double[] a_tmp = Arrays.copyOf(a, a.length);
		double[] b_tmp = Arrays.copyOf(b, b.length);
		return new ONDEParameterDistribution(m, n, off_tmp, p_tmp, alpha_tmp, a_tmp, b_tmp);
	}
	
	/**
	 * Constructor
	 * @param m the number of Gaussian components in the mixture
	 * @param n the number of Laplace components in the mixture
	 * @param offset Offset values (L-dimensional array)
	 * @param p the parameters of Dirichlet prior/posterior distribution of the offset mixing coefficients \omega (L-dimensional array)
	 * @param alpha the parameters of Dirichlet prior/posterior distribution of the mixing coefficients \pi ((m+n)-dimensional array)
	 * @param a the parameters `a' of Gamma prior/posterior distribution of the accuracy parameters \eta ((m+n)-dimensional array)
	 * @param b the parameters `b' of Gamma prior/posterior distribution of the accuracy parameters \eta  ((m+n)-dimensional array)
	 */
	public ONDEParameterDistribution(int m, int n, double offset[], double[] p, double[] alpha, double[] a, double[] b){
		initialize(m, n, offset, p, alpha, a, b);
	}
	
	/**
	 * Constructor
	 * @param m the number of Gaussian components in the mixture
	 * @param n the number of Laplace components in the mixture
	 * @param offset Offset values (L-dimensional array)
	 * @param p0 the parameters of Dirichlet prior/posterior distribution of the offset mixing coefficients \omega (common)
	 * @param alpha0 the parameters of Dirichlet prior/posterior distribution of the mixing coefficients \pi (Common)
	 * @param a the parameters `a' of Gamma prior/posterior distribution of the accuracy parameters \eta ((m+n)-dimensional array)
	 * @param b the parameters `b' of Gamma prior/posterior distribution of the accuracy parameters \eta  ((m+n)-dimensional array)
	 */
	public ONDEParameterDistribution(int m, int n,double offset[], double p0, double alpha0, double[] a, double[] b){
		double tmpAlpha[] = new double[m+n];
		Arrays.fill(tmpAlpha, alpha0);
		double tmpP[] = new double[offset.length];
		Arrays.fill(tmpP, p0);
		initialize(m, n, offset,tmpP,tmpAlpha, a, b);
	}
	

	/**
	 * Initialization
	 * @param m the number of Gaussian components in the mixture
	 * @param n the number of Laplace components in the mixture
	 * @param offset Offset values (L-dimensional array)
	 * @param p the parameters of Dirichlet prior/posterior distribution of the offset mixing coefficients \omega (L-dimensional array)
	 * @param alpha the parameters of Dirichlet prior/posterior distribution of the mixing coefficients \pi ((m+n)-dimensional array)
	 * @param a the parameters `a' of Gamma prior/posterior distribution of the accuracy parameters \eta ((m+n)-dimensional array)
	 * @param b the parameters `b' of Gamma prior/posterior distribution of the accuracy parameters \eta  ((m+n)-dimensional array)
	 */
	private void initialize(int m, int n, double offset[], double p[], double[] alpha, double[] a, double[] b){
		if(m<0 || n<0){
			log.error("The number m and n should not be negative.");
			throw new IllegalArgumentException();
		}
		this.m = m;
		this.n = n;
		int K = m+n;
		if(alpha==null || alpha.length != K || a==null || a.length != K || b==null || b.length != K || p==null || offset==null || p.length!=offset.length){
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
		for(int i=0;i<p.length;i++){
			if(p[i]<=0){
				log.error("p[" + i + "] should be positive.");
				throw new IllegalArgumentException();
			}
		}
		this.offset = Arrays.copyOf(offset, offset.length);
		this.alpha = Arrays.copyOf(alpha, alpha.length);
		this.p = Arrays.copyOf(p, p.length);
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
		dirichlet_alpha = new Dirichlet(alpha);
		dirichlet_p = new Dirichlet(p);
	}
	
	/**
	 * MAP estimation of parameter values.
	 * @return MAP estimation of parameter values.
	 */
	public ONDE MAPEstimation(){
		log.info("Variational Bayesian MAP estimation process started.");
		int K = m+n;
		double pi[] = new double[K];
		double sum = 0;
		for(int i=0;i<K;i++) sum += alpha[i];
		for(int i=0;i<K;i++) pi[i] = (alpha[i]-1)/(sum-K);
		double omega[] = new double[p.length];
		sum = 0;
		for(int i=0;i<p.length;i++) sum += p[i];
		for(int i=0;i<p.length;i++) omega[i] = (p[i]-1)/(sum-p.length);
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
		ONDE para = new ONDE(pi,sigma,lambda, offset,omega);
		log.info("Variational Bayesian MAP estimation process finished.");
		if(log.isInfoEnabled()){
			StringBuilder sb = new StringBuilder();
			for(int i=0;i<offset.length;i++){
				sb.append(i+1);
				sb.append(":offset=");
				sb.append(offset[i]);
				sb.append(",");
				sb.append(omega[i]);
			}
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
	public ONDE ptEstimation(){
		log.info("Variational Bayesian point estimation process started.");
		int K = m+n;
		double pi[] = new double[K];
		double sum = 0;
		for(int i=0;i<K;i++) sum += alpha[i];
		for(int i=0;i<K;i++) pi[i] = alpha[i]/sum;
		double omega[] = new double[p.length];
		sum = 0;
		for(int i=0;i<p.length;i++) sum += p[i];
		for(int i=0;i<p.length;i++) omega[i] = p[i]/sum;
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
		ONDE para = new ONDE(pi,sigma,lambda, offset,omega);
		log.info("Variational Bayesian point estimation process finished.");
		if(log.isInfoEnabled()){
			StringBuilder sb = new StringBuilder();
			for(int i=0;i<offset.length;i++){
				sb.append(i+1);
				sb.append(":offset=");
				sb.append(offset[i]);
				sb.append(",omega=");
				sb.append(omega[i]);
			}
			
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
	 * Estimate the probability that each feasible offset value is applied.
	 * @param x Observation
	 * @return Probability that each feasible offset value is applied.
	 */
	public double[] estimateOffsetProbability(double x){
		double data[] = new double[1];
		data[0] = x;
		int K = m + n;
		int N = data.length;
		int L = offset.length;
		double r[][][] = new double[N][][];
		for(int i=0;i<N;i++){
			r[i] = new double[K][];
			for(int k=0;k<K;k++){
				r[i][k] = new double[L];
				for(int l=0;l<L;l++) r[i][k][l] = 0;
			}
		}
		
		double Rk[] = new double[K];
		double lowerbound = Double.NEGATIVE_INFINITY;
		double lowerbound_pre = Double.NEGATIVE_INFINITY;
		ONDEParameterDistribution params = this.clone();
		ONDEVB engine = new ONDEVB();
		boolean flag = true;
		do{
			// E step
			engine.Estep(params,data,r);
			// M step
			engine.Mstep(this,data,r,Rk,params);
			// lower bound
			lowerbound_pre = lowerbound;
			lowerbound = engine.lowerbound(this,params,r,Rk);
			// Stop condition
			if(Math.abs(lowerbound - lowerbound_pre) < ONDEVB.threshold) flag = false;
			else if(lowerbound_pre > lowerbound) throw new ArithmeticException();
		}while(flag);
		double w[] = new double[L];
		for(int l=0;l<L;l++){
			w[l] = 0;
			for(int k=0;k<K;k++) w[l] += r[0][k][l];
		}
		return w;
	}
	
	/**
	 * Sampling from the parameter distribution.
	 * @return sample
	 */
	public ONDE sampling(){
		double sigma[] = new double[m];
		double lambda[] = new double[n];
		for(int i=0;i<m;i++) sigma[i] = 1/Math.sqrt(gamma[i].nextDouble());
		for(int i=0;i<n;i++) lambda[i] = 1/gamma[i+m].nextDouble();
		return new ONDE(dirichlet_alpha.nextDouble(), sigma, lambda,offset,dirichlet_p.nextDouble());
	}
	
}
