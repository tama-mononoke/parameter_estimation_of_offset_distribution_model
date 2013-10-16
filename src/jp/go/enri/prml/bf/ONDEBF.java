/**
 * 
 */
package jp.go.enri.prml.bf;


import jp.go.enri.prml.dist.ONDE;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * Maximum likelihood estimation of the offset mixture distribution by means of brute-force search.
 * @author Masato Fujita (Electronic Navigation Research Institute)
 * @version 1.0.1 (Last update: 30/11/2011)
 *
 */
public class ONDEBF {
	/**
	 * Log
	 */
	public static Log log = LogFactory.getLog(ONDEBF.class);
	/**
	 * Threshold for log-likelihood
	 */
	public static double logLiklihood_threshold = 10E-3;
	/**
	 * Result of the brute-force search algorithm
	 * @author Masato Fujita
	 *
	 */
	public static class Result{
		/**
		 * Get the parameters of distribution model.
		 * @return the parameters of distribution model
		 */
		public ONDE getParam() {
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
		ONDE param;
		/**
		 * log-likelihood
		 */
		double logLiklihood;
	}
	
	/**
	 * Used for the parameter value increments
	 * @author Masato Fujita
	 */
	static class Increment{
		/**
		 * Log
		 */
		public static Log log = LogFactory.getLog(Increment.class);
		/**
		 * Division number
		 */
		private int D;
		/**
		 * the number of Gaussian components in the mixture distribution.
		 */
		private int m;
		/**
		 * the number of Laplace components in the mixture distribution.
		 */
		private int n;
		/**
		 * Range of omega (offset mixing coefficients)
		 */
		private double omega_bound[][];
		/**
		 * Range of alpha (mixing coefficients)
		 */
		private double alpha_bound[][];
		/**
		 * Range of sigma (standard deviation of Gaussian components)
		 */
		private double sigma_bound[][];
		/**
		 * Range of lambda (scale parameters of Laplace distributions)
		 */
		private double lambda_bound[][];
		/**
		 * Pointer
		 */
		private int pointer[][];
		/**
		 * Feasible offset values
		 */
		private double offset[];
		/**
		 * Constructor
		 * @param d Division number
		 * @param offset Feasible offset values
		 * @param omega_bound Range of omega (offset mixing coefficients)
		 * @param alpha_bound Range of alpha (mixing coefficients)
		 * @param sigma_bound Range of sigma (standard deviation of Gaussian components)
		 * @param lambda_bound Range of lambda (scale parameters of Laplace distributions)
		 */
		Increment(int d, double offset[], double[][] omega_bound, double[][] alpha_bound, double[][] sigma_bound,
				double[][] lambda_bound) {
			super();
			D = d;
			this.offset = offset;
			this.omega_bound = omega_bound;
			this.alpha_bound = alpha_bound;
			this.sigma_bound = sigma_bound;
			this.lambda_bound = lambda_bound;
			m = sigma_bound.length;
			n = lambda_bound.length;
			pointer = new int[4][];
			pointer[0] = new int[omega_bound.length];
			pointer[1] = new int[m+n-1];
			pointer[2] = new int[m];
			pointer[3] = new int[n];
			for(int i=0;i<4;i++){
				for(int j=0;j<pointer[i].length;j++) pointer[i][j] = 0;
			}
		}
		/**
		 * Move the pointer.
		 * @return false if the pointer arrives at the end point, and true otherwise.
		 */
		boolean increment(){
			while(increment2()){
				double sum = 0;
				for(int i=0;i<alpha_bound.length;i++){
					sum += alpha_bound[i][0]*(1-(2*pointer[1][i]+1)/(2*D)) + alpha_bound[i][1]*(2*pointer[1][i]+1)/(2*D);
				}
				if(sum<=1){
					sum = 0;
					for(int i=0;i<omega_bound.length;i++){
						sum += omega_bound[i][0]*(1-(2*pointer[0][i]+1)/(2*D)) + omega_bound[i][1]*(2*pointer[0][i]+1)/(2*D);
					}
					if(sum<=1) return true;
				}
			}
			return false;	
		}
		
		/**
		 * Move the pointer.
		 * @return false if the pointer arrives at the end point, and true otherwise.
		 */
		private boolean increment2(){
			for(int i=pointer.length-1;i>=0;i--){
				int pt[] = pointer[i];
				for(int j=pt.length-1;j>=0;j--){
					pt[j]++;
					if(pt[j]==D){
						pt[j] = 0;
					}
					else return true;
				}
			}
			return false;	
		}
		
		/**
		 * Parameter values at the position of the pointer.
		 * @return Parameter values
		 */
		ONDE getCurrent(){
			double omega[] = new double[offset.length];
			double sum = 0;
			for(int i=0;i<offset.length-1;i++){
				omega[i] = omega_bound[i][0]*(1-(2.0*pointer[0][i]+1)/(2*D)) + omega_bound[i][1]*(2.0*pointer[0][i]+1)/(2*D);
				sum += omega[i];
			}
			omega[offset.length-1] = 1-sum;
			double alpha[] = new double[m+n];
			sum = 0;
			for(int i=0;i<m+n-1;i++){
				alpha[i] = alpha_bound[i][0]*(1-(2.0*pointer[1][i]+1)/(2*D)) + alpha_bound[i][1]*(2.0*pointer[1][i]+1)/(2*D);
				sum += alpha[i];
			}
			alpha[m+n-1] = 1-sum;
			double sigma[] = new double[m];
			for(int i=0;i<m;i++){
				sigma[i] = sigma_bound[i][0]*(1-(2.0*pointer[2][i]+1)/(2*D)) + sigma_bound[i][1]*(2.0*pointer[2][i]+1)/(2*D);
			}
			double lambda[] = new double[n];
			for(int i=0;i<n;i++){
				lambda[i] = lambda_bound[i][0]*(1-(2.0*pointer[3][i]+1)/(2*D)) + lambda_bound[i][1]*(2.0*pointer[3][i]+1)/(2*D);
			}
			return new ONDE(alpha,sigma,lambda,offset,omega);
		}
		
		Increment getNewIncrement(){
			double tmp_omega_bound[][] = new double[omega_bound.length][];
			double tmp_alpha_bound[][] = new double[alpha_bound.length][];
			double tmp_sigma_bound[][] = new double[sigma_bound.length][];
			double tmp_lambda_bound[][] = new double[lambda_bound.length][];
			for(int i=0;i<alpha_bound.length;i++){
				tmp_omega_bound[i] = new double[2];
				tmp_omega_bound[i][0] = omega_bound[i][0]*(1- (double) pointer[0][i]/D) + omega_bound[i][1]*((double) pointer[0][i]/D);
				tmp_omega_bound[i][1] = omega_bound[i][0]*(1- (double) (pointer[0][i]+1)/D) + omega_bound[i][1]*((double) (pointer[0][i]+1)/D);
			}
			for(int i=0;i<alpha_bound.length;i++){
				tmp_alpha_bound[i] = new double[2];
				tmp_alpha_bound[i][0] = alpha_bound[i][0]*(1- (double) pointer[1][i]/D) + alpha_bound[i][1]*((double) pointer[1][i]/D);
				tmp_alpha_bound[i][1] = alpha_bound[i][0]*(1- (double) (pointer[1][i]+1)/D) + alpha_bound[i][1]*((double) (pointer[1][i]+1)/D);
			}
			for(int i=0;i<sigma_bound.length;i++){
				tmp_sigma_bound[i] = new double[2];
				tmp_sigma_bound[i][0] = sigma_bound[i][0]*(1- (double) pointer[2][i]/D) + sigma_bound[i][1]*((double) pointer[2][i]/D);
				tmp_sigma_bound[i][1] = sigma_bound[i][0]*(1- (double) (pointer[2][i]+1)/D) + sigma_bound[i][1]*((double) (pointer[2][i]+1)/D);
			}
			for(int i=0;i<lambda_bound.length;i++){
				tmp_lambda_bound[i] = new double[2];
				tmp_lambda_bound[i][0] = lambda_bound[i][0]*(1- (double) pointer[3][i]/D) + lambda_bound[i][1]*((double) pointer[3][i]/D);
				tmp_lambda_bound[i][1] = lambda_bound[i][0]*(1- (double) (pointer[3][i]+1)/D) + lambda_bound[i][1]*((double) (pointer[3][i]+1)/D);
			}
			return new Increment(D, offset,tmp_omega_bound,tmp_alpha_bound, tmp_sigma_bound, tmp_lambda_bound);
		}
		
		void print(){
			if(log.isInfoEnabled()){
				StringBuilder sb = new StringBuilder();
				sb.append("omega:");
				for(int i=0;i<omega_bound.length;i++){
					sb.append("[" + i + ":" + omega_bound[i][0] + "," + omega_bound[i][1] + "],");
				}
				sb.append("alpha:");
				for(int i=0;i<alpha_bound.length;i++){
					sb.append("[" + i + ":" + alpha_bound[i][0] + "," + alpha_bound[i][1] + "],");
				}
				sb.append("sigma:");
				for(int i=0;i<sigma_bound.length;i++){
					sb.append("[" + i + ":" + sigma_bound[i][0] + "," + sigma_bound[i][1] + "],");
				}
				sb.append("lambda:");
				for(int i=0;i<lambda_bound.length;i++){
					sb.append("[" + i + ":" + lambda_bound[i][0] + "," + lambda_bound[i][1] + "],");
				}
				log.info(sb.toString());
			}
		}
	}
	
	/**
	 * parameter estimation by means of the brute-force search. 
	 * @param initial_params inital parameter
	 * @param data observations
	 * @return estimations
	 */
	public Result estimate(ONDEParameter initial_params, double[] data){
		log.info("BF method estimation process started.");
		// initialization
		log.info("BF method initialization process started.");
		int N = data.length;
		double initial_alpha_bound[][] = new double[initial_params.initial_sigma_bound.length+initial_params.initial_lambda_bound.length-1][];
		for(int i=0;i<initial_alpha_bound.length;i++){
			initial_alpha_bound[i] = new double[2];
			for(int j=0;j<2;j++) initial_alpha_bound[i][j] = j;
		}
		double initial_omega_bound[][] = new double[initial_params.offset.length-1][];
		for(int i=0;i<initial_omega_bound.length;i++){
			initial_omega_bound[i] = new double[2];
			for(int j=0;j<2;j++) initial_omega_bound[i][j] = j;
		}
		double logLiklihood_pre = Double.NEGATIVE_INFINITY;
		Increment incre = new Increment(initial_params.D, initial_params.offset,initial_omega_bound,initial_alpha_bound, initial_params.initial_sigma_bound, initial_params.initial_lambda_bound);
		log.info("BF method initialization process finished.");
		boolean flag = true;
		do{
			incre.print();
			log.info("BF increment process started.");
			double maxlogLiklihood = Double.NEGATIVE_INFINITY;
			Increment next_incre = null;
			do{
				double logLiklihood = logLiklihood(incre.getCurrent(),data);
				if(logLiklihood > maxlogLiklihood){
					maxlogLiklihood = logLiklihood;
					next_incre = incre.getNewIncrement();
				}
			}
			while(incre.increment());
			// Stopping condition 
			if(Math.abs(maxlogLiklihood - logLiklihood_pre) < logLiklihood_threshold/N) flag = false;
			else{
				incre = next_incre;
			}
			logLiklihood_pre = maxlogLiklihood;
		}while(flag);
		Result res = new Result();
		res.param = incre.getCurrent();
		res.logLiklihood = logLiklihood_pre;
		log.info("BF method estimation process finished.");
		return res;
	}
	
	
	/**
	 * Compute the log-likelihood.
	 * @param param parameters of mixture distribution
	 * @param data observations
	 * @return log-likelihood
	 */
	double logLiklihood(ONDE param,double data[]){
		//log.info("BF method log-liklihood evaluation process started.");
		int N = data.length;
		double tmp = 0;
		double off[] = param.getOffset();
		double pi[] = param.getPi();
		double omega[] = param.getOmega();
		double lambda[] = param.getLambda();
		double sigma[] = param.getSigma();
		for(int i=0;i<N;i++){
			for(int l=0;l<param.getL();l++){
				double tmp2 = 0;
				double x = data[i];
				for(int k=0;k<param.getM();k++){
					tmp2 += omega[l]*pi[k]*Math.exp(-0.5*Math.pow(x-off[l], 2)/Math.pow(sigma[k],2))/(Math.sqrt(2*Math.PI)*sigma[k]);
				}
				for(int k=0;k<param.getN();k++){
					tmp2 += omega[l]*pi[k+param.getM()]*Math.exp(-Math.abs(x-off[l])/lambda[k])/(2*lambda[k]);
				}
				tmp += Math.log(tmp2);
			}
		}
		//log.info("BF method log-liklihood evaluation process finished.");
		return tmp;
	}

}
