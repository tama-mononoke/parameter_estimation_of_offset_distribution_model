/**
 * 
 */
package jp.go.enri.prml.em;


import jp.go.enri.prml.dist.NDE;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * Parameter estimation of the mixture distribution by means of the EM algorithm.
 * @author Masato Fujita (Electronic Navigation Research Institute) 
 * @version 1.0.1　(Last update: 30/11/2011)
 *
 */
public class NDEEM {
	/**
	 * Log
	 */
	public static Log log = LogFactory.getLog(NDEEM.class);
	/**
	 * Threshold for log-likelihood
	 */
	public static double logLiklihood_threshold = 10E-3;
	/**
	 * Result of EM algorithm
	 * @author Masato Fujita
	 *
	 */
	public static class Result{
		/**
		 * Get the parameters of distribution model.
		 * @return the parameters of distribution model
		 */
		public NDE getParam() {
			return param;
		}
		/**
		 * Get the log-likelihood.
		 * @return log-likelihood
		 */
		public double getLogLiklihood() {
			return logLiklihood;
		}
		/**
		 * the parameters of distribution model
		 */
		NDE param;
		/**
		 * log-likelihood
		 */
		double logLiklihood;
	}
	
	/**
	 * Parameter estimation by means of EM algorithm
	 * @param initial_params initial parameter
	 * @param data observations
	 * @return estimations
	 */
	public Result estimate(NDE initial_params, double[] data){
		log.info("EM method estimation process started.");
		// Initialization
		log.info("EM method initialization process started.");
		int N = data.length;
		int K = initial_params.getM() + initial_params.getN();
		double gamma[][] = new double[N][];
		for(int i=0;i<gamma.length;i++){
			gamma[i] = new double[K];
		}
		double logLiklihood = Double.NEGATIVE_INFINITY;
		double logLiklihood_pre = Double.NEGATIVE_INFINITY;
		NDE params = initial_params.clone();
		log.info("EM method initialization process finished.");
		boolean flag = true;
		do{
			// E step
			Estep(params,data,gamma);
			// M step
			params = Mstep(initial_params.getM(),initial_params.getN(),data,gamma);
			// log-likelihood
			logLiklihood_pre = logLiklihood;
			logLiklihood = logLiklihood(params,data);
			// stopping condition
			if(Math.abs(logLiklihood - logLiklihood_pre) < logLiklihood_threshold/N) flag = false;
			else if(logLiklihood_pre > logLiklihood) throw new ArithmeticException();
		}while(flag);
		Result res = new Result();
		res.param = params;
		res.logLiklihood = logLiklihood;
		log.info("EM method estimation process finished.");
		return res;
	}
	
	/**
	 * E-step. Compute the burden rate. 
	 * @param param Parameter of mixture distributions
	 * @param data observations
	 * @param gamma burden rates
	 */
	void Estep(NDE param, double[] data, double gamma[][]){
		log.info("EM method E-step process started.");
		int N = data.length;
		int K = param.getM() + param.getN();
		for(int i=0;i<N;i++){
			double sum = 0;
			double x = data[i];
			for(int k=0;k<param.getM();k++){
				gamma[i][k] = param.getPi()[k]*Math.exp(-0.5*x*x/(param.getSigma()[k]*param.getSigma()[k]))/(Math.sqrt(2*Math.PI)*param.getSigma()[k]);
				sum += gamma[i][k];
			}
			for(int k=0;k<param.getN();k++){
				gamma[i][k+param.getM()] = param.getPi()[k+param.getM()]*Math.exp(-Math.abs(x)/param.getLambda()[k])/(2*param.getLambda()[k]);
				sum += gamma[i][k+param.getM()];
			}
			for(int k=0;k<K;k++) gamma[i][k] /= sum;
		}
		if(log.isDebugEnabled()){
			log.debug("gamma values:" + gamma.length);
			for(int i=0;i<gamma.length;i++){
				StringBuilder sb = new StringBuilder();
				sb.append(data[i]);
				for(int j=0;j<gamma[i].length;j++){
					sb.append(", ");
					sb.append(gamma[i][j]);
				}
				log.debug(sb.toString());
			}
		}
		log.info("EM method E-step process finished.");
	}
	
	/**
	 * M-step
	 * @param m the number of Gaussian components in the mixture distribution.
	 * @param n the number of Laplace components in the mixture distribution.
	 * @param data observation
	 * @param gamma burden rate
	 * @return parameters of mixture distribution
	 */
	NDE Mstep(int m, int n, double data[],double gamma[][]){
		log.info("EM method M-step process started.");
		int N = data.length;
		int K = m+n;
		double Nk[] = new double[K];
		double pi[] = new double[K];
		for(int i=0;i<K;i++){
			Nk[i] = 0;
			for(int j=0;j<N;j++) Nk[i] += gamma[j][i];
			pi[i] = Nk[i]/N;
		}
		double sigma[] = new double[m];
		for(int k=0;k<m;k++){
			sigma[k] = 0;
			for(int i=0;i<N;i++){
				sigma[k] += gamma[i][k]*data[i]*data[i];
			}
			sigma[k] /= Nk[k];
			sigma[k] = Math.sqrt(sigma[k]);
		}
		double lambda[] = new double[n];
		for(int k=0;k<n;k++){
			lambda[k] = 0;
			for(int i=0;i<N;i++){
				lambda[k] += gamma[i][m+k]*Math.abs(data[i]);
			}
			lambda[k] /= Nk[m+k];
		}
		log.info("EM method M-step process finished.");
		return new NDE(pi,sigma,lambda);
	}
	/**
	 * Compute the log-likelihood.
	 * @param param parameters of mixture distribution
	 * @param data observations
	 * @return log-likelihood
	 */
	double logLiklihood(NDE param,double data[]){
		log.info("EM method log-liklihood evaluation process started.");
		int N = data.length;
		double tmp = 0;
		for(int i=0;i<N;i++){
			double tmp2 = 0;
			double x = data[i];
			for(int k=0;k<param.getM();k++){
				tmp2 += param.getPi()[k]*Math.exp(-0.5*x*x/(param.getSigma()[k]*param.getSigma()[k]))/(Math.sqrt(2*Math.PI)*param.getSigma()[k]);
			}
			for(int k=0;k<param.getN();k++){
				tmp2 += param.getPi()[k+param.getM()]*Math.exp(-Math.abs(x)/param.getLambda()[k])/(2*param.getLambda()[k]);
			}
			tmp += Math.log(tmp2);
		}
		log.info("EM method log-liklihood evaluation process finished.");
		return tmp;
	}

}
