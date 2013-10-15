/**
 * 
 */
package jp.go.enri.prml.vb;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.math.special.Gamma;

/**
 * Parameter estimation of the mixture distribution by means of the variational Bayesian method.
 * @author Masato Fujita (Electronic Navigation Research Institute) 
 * @version 1.0.1　(Last update: 30/11/2011)
 *
 */
public class NDEVB {
	/**
	 * Log
	 */
	public static Log log = LogFactory.getLog(NDEVB.class);
	/**
	 * Very small number. Computation is terminated when the lower bound is smaller than this value.
	 */
	public static double threshold = 10E-3; 
	
	/**
	 * Result of variational Bayes
	 * @author Masato Fujita
	 *
	 */
	public static class Result{
		/**
		 * Get the parameter for variational Bayes.
		 * @return the parameter for variational Bayes.
		 */
		public NDEParameterDistribution getParam() {
			return param;
		}
		/**
		 * Get the lower bound.
		 * @return the lower bound
		 */
		public double getLowerbound() {
			return lowerbound;
		}
		/**
		 * the parameter for variational Bayes.
		 */
		NDEParameterDistribution param;
		/**
		 * the lower bound
		 */
		double lowerbound;
	}
	
	/**
	 * Parameter estimation by means of variational Bayes
	 * @param initial_params prior distribution
	 * @param dataMatrix data matrix
	 * @return estimation result
	 */
	public Result estimate(NDEParameterDistribution initial_params, double data[]){
		log.info("Variational Bayesian estimation process started.");
		// Innitialization
		log.info("Variational Bayesian initialization process started.");
		double r[][] = new double[data.length][];
		int K = initial_params.m + initial_params.n;
		int N = data.length;
		for(int i=0;i<r.length;i++){
			r[i] = new double[K];
		}
		double Nk[] = new double[K];
		double lowerbound = Double.NEGATIVE_INFINITY;
		double lowerbound_pre = Double.NEGATIVE_INFINITY;
		NDEParameterDistribution params = initial_params.clone();
		log.info("Variational Bayesian initialization process finished.");
		boolean flag = true;
		do{
			// E step
			Estep(params,data,r);
			// M step
			Mstep(initial_params,data,r,Nk,params);
			// Evaluate lower bound.
			lowerbound_pre = lowerbound;
			lowerbound = lowerbound(initial_params,params,r,Nk);
			// Stop condition
			if(Math.abs(lowerbound - lowerbound_pre) < threshold/N) flag = false;
			else if(lowerbound_pre > lowerbound) throw new ArithmeticException();
		}while(flag);
		params.updateOtherParameters();
		Result res = new Result();
		res.param = params;
		res.lowerbound = lowerbound;
		log.info("Variational Bayesian estimation process finished.");
		return res;
	}
	
	/**
	 * E-step
	 * @param param Parameter distribution(IN)
	 * @param data data matrix(IN)
	 * @param r Array for stocking the computation result of the variable `r.'(OUT)
	 */
	void Estep(NDEParameterDistribution param, double data[], double r[][]){
		log.info("Variational Bayesian E-step process started.");
		int K = param.m + param.n;
		// Compute the terms which are independent of the value of i.
		double tmp1[] = new double[K];
		for(int k=0;k<param.m;k++){
			tmp1[k] = Gamma.digamma(param.alpha[k])+0.5*(Gamma.digamma(param.a[k])-Math.log(param.b[k])-Math.log(2*Math.PI));
		}
		for(int k=param.m;k<K;k++){
			tmp1[k] = Gamma.digamma(param.alpha[k])+(Gamma.digamma(param.a[k])-Math.log(param.b[k])-Math.log(2));
		}
		int N = data.length;
		for(int i=0;i<N;i++){
			// Since rho_{nk} is too small, the maximum of log(rho_{nk}) is saved.
			double rhoMax = Double.NEGATIVE_INFINITY;
			for(int k=0;k<param.m;k++){
				double tmp = tmp1[k]-0.5*param.a[k]*Math.pow(data[i],2)/param.b[k];
				r[i][k] = tmp;
				rhoMax = Math.max(rhoMax, tmp);
			}
			for(int k=param.m;k<K;k++){
				double tmp = tmp1[k]-param.a[k]*Math.abs(data[i])/param.b[k]; 
				r[i][k] = tmp;
				rhoMax = Math.max(rhoMax, tmp);
			}
			double sum = 0;
			for(int j=0;j<K;j++){
				double tmp = Math.exp(r[i][j]-rhoMax);
				r[i][j] = tmp;
				sum += tmp;
			}
			for(int j=0;j<K;j++){
				r[i][j] /= sum; 
			}
		}
		if(log.isDebugEnabled()){
			log.debug("r values:" + r.length);
			for(int i=0;i<r.length;i++){
				StringBuilder sb = new StringBuilder();
				sb.append(data[i]);
				for(int j=0;j<r[i].length;j++){
					sb.append(", ");
					sb.append(r[i][j]);
				}
				log.debug(sb.toString());
			}
		}
		log.info("Variational Bayesian E-step process finished.");
	}
	
	/**
	 * M-step
	 * @param initial_param Prior distribution(IN)
	 * @param data data matrix(IN)
	 * @param r The variable `r'(IN)
	 * @param Nk \sum_{i=1}^N r_{i,k}(OUT)
	 * @param param Parameter distribution(OUT)
	 */
	void Mstep(NDEParameterDistribution initial_param, double data[], double r[][], double Nk[], NDEParameterDistribution param){
		log.info("Variational Bayesian M-step process started.");
		int N = data.length;
		int K = initial_param.m + initial_param.n;
		// \sum_{i=1}^N r_{i,k}
		for(int i=0;i<K;i++){
			Nk[i] = 0;
			for(int j=0;j<N;j++) Nk[i] += r[j][i];
		}
		// \sum_{i=1}^N r_{i,k}x_i^2と\sum_{i=1}^N r_{i,k}|x_i|
		double xk[] = new double[K];
		for(int k=0;k<param.m;k++){
			xk[k] = 0;
			for(int i=0;i<N;i++) xk[k] += r[i][k]*Math.pow(data[i], 2);
		}
		for(int k=param.m;k<K;k++){
			xk[k] = 0;
			for(int i=0;i<N;i++) xk[k] += r[i][k]*Math.abs(data[i]);
		}
		// Update Parameters.
		for(int k=0;k<param.m;k++){
			param.alpha[k] = initial_param.alpha[k] + Nk[k]; //(6)
			param.a[k] = initial_param.a[k] + Nk[k]/2; //(7)
			param.b[k] = initial_param.b[k] + xk[k]/2; //(8)
		}
		for(int k=param.m;k<K;k++){
			param.alpha[k] = initial_param.alpha[k] + Nk[k]; //(6)
			param.a[k] = initial_param.a[k] + Nk[k]; //(7)
			param.b[k] = initial_param.b[k] + xk[k]; //(8)
		}
		log.info("Variational Bayesian M-step process finished.");
	}
	
	/**
	 * Compute the lower bound.
	 * @param initial_param Prior distribution(IN)
	 * @param r The values of the parameter `r'(IN)
	 * @param Nk \sum_{i=1}^N r_{i,k}(IN)
	 * @param param Parameter distribution(IN)
	 * @return Lower bound
	 */
	double lowerbound(NDEParameterDistribution initial_param, NDEParameterDistribution param, double r[][], double Nk[]){
		log.info("Variational Bayesian lower bound calculation process started.");
		int N = r.length;
		int K = Nk.length;
		double alpha_hat = 0;
		double alpha0_hat = 0;
		for(int k=0;k<K;k++){
			alpha0_hat += initial_param.alpha[k];
			alpha_hat += param.alpha[k];
		}
		double tmp = Gamma.logGamma(alpha0_hat) - Gamma.logGamma(alpha_hat); //Terms of (9) independent of k
		for(int k=0;k<K;k++){
			double tmp2 = 0;
			for(int i=0;i<N;i++){
				if(r[i][k]!=0){
					tmp2 += r[i][k]*Math.log(r[i][k]);
				}
			}
			tmp += - Gamma.logGamma(initial_param.alpha[k]) + Gamma.logGamma(param.alpha[k])
				+initial_param.a[k]*Math.log(initial_param.b[k]) - param.a[k]*Math.log(param.b[k])
				-Gamma.logGamma(initial_param.a[k]) + Gamma.logGamma(param.a[k])-tmp2
				- ((k<initial_param.m) ? Math.log(2*Math.PI)*Nk[k]/2 : Math.log(2)*Nk[k]); //(9)
		}
		log.info("Variational Bayesian lower bound calculation process finished.");
		log.info("Variational lower bound: " + tmp);		
		return tmp;
	}
	
}
