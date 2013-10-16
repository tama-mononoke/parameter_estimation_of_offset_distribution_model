/**
 * 
 */
package jp.go.enri.prml.vb;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.math.special.Gamma;

/**
 * Parameter estimation of the offset mixture distribution by means of the variational Bayesian method.
 * @author Masato Fujita (Electronic Navigation Research Institute)
 * @version 1.0.1 (Last update: 06/12/2011)
 *
 */
public class ONDEVB {
	/**
	 * Log
	 */
	public static Log log = LogFactory.getLog(ONDEVB.class);
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
		public ONDEParameterDistribution getParam() {
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
		ONDEParameterDistribution param;
		/**
		 * the lower bound
		 */
		double lowerbound;
	}
	
	/**
	 * Parameter estimation by means of variational Bayes
	 * @param initial_params prior distribution
	 * @param data data matrix
	 * @return estimation result
	 */
	public Result estimate(ONDEParameterDistribution initial_params, double data[]){
		log.info("Variational Bayesian estimation process started.");
		// Initialization
		log.info("Variational Bayesian initialization process started.");
		int K = initial_params.m + initial_params.n;
		int N = data.length;
		int L = initial_params.offset.length;
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
		ONDEParameterDistribution params = initial_params.clone();
		log.info("Variational Bayesian initialization process finished.");
		boolean flag = true;
		do{
			// E step
			Estep(params,data,r);
			// M step
			Mstep(initial_params,data,r,Rk,params);
			// lower bound
			lowerbound_pre = lowerbound;
			lowerbound = lowerbound(initial_params,params,r,Rk);
			// stopping condition
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
	
	void Estep(ONDEParameterDistribution param, double data[], double r[][][]){
		log.info("Variational Bayesian E-step process started.");
		int K = param.m + param.n;
		int L = param.offset.length;
		double tmp_t[] = new double[L];
		double tmp_r[] = new double[K];
		// (4)
		// Evaluate the terms only dependent on l.
		for(int l=0;l<L;l++){
			tmp_t[l] = Gamma.digamma(param.p[l]); //(4)
		}
		// Compute the terms only dependent on k.
		for(int k=0;k<param.m;k++){
			tmp_r[k] = Gamma.digamma(param.alpha[k]) + 0.5*(Gamma.digamma(param.a[k])-Math.log(param.b[k])-Math.log(2*Math.PI)); //(4)
		}
		for(int k=param.m;k<K;k++){
			tmp_r[k] = Gamma.digamma(param.alpha[k]) + (Gamma.digamma(param.a[k])-Math.log(param.b[k])-Math.log(2)); //(4)
		}
		
		int N = data.length;
		for(int i=0;i<N;i++){
			// Since rho_{nk} is too small, the maximum of log(rho_{nk}) is saved.
			double rhoMax = Double.NEGATIVE_INFINITY;
			for(int l=0;l<L;l++){
				for(int k=0;k<param.m;k++){
					r[i][k][l] = tmp_t[l]+tmp_r[k]-0.5*param.a[k]*Math.pow(data[i]-param.offset[l],2)/param.b[k]; //(4)
					rhoMax = Math.max(rhoMax, r[i][k][l]);
				}
				for(int k=param.m;k<K;k++){
					r[i][k][l] = tmp_t[l]+tmp_r[k]-param.a[k]*Math.abs(data[i]-param.offset[l])/param.b[k]; //(4)
					rhoMax = Math.max(rhoMax, r[i][k][l]);
				}
			}
			double sum = 0;
			for(int l=0;l<L;l++){
				for(int k=0;k<K;k++){
					double tmp = Math.exp(r[i][k][l]-rhoMax);
					r[i][k][l] = tmp;
					sum += tmp;
				}
			}
			for(int l=0;l<L;l++){
				for(int k=0;k<K;k++){
					r[i][k][l] /= sum;
				}
			}
		}
		log.info("Variational Bayesian E-step process finished.");
	}
	
	
	/**
	 * M-step
	 * @param initial_param prior distribution(IN)
	 * @param data data matrix(IN)
	 * @param r Values of `r7(IN)
	 * @param Rk \sum_{i=1}^N \sum_{l=1}^L r_{i,k,l}(OUT)
	 * @param param Parameter distribution(OUT)
	 */
	void Mstep(ONDEParameterDistribution initial_param, double data[], double r[][][], double Rk[], ONDEParameterDistribution param){
		log.info("Variational Bayesian M-step process started.");
		int N = data.length;
		int K = initial_param.m + initial_param.n;
		int L = initial_param.offset.length;
		// R_k
		for(int k=0;k<K;k++){
			Rk[k] = 0;
			for(int l=0;l<L;l++){
				double tmp = 0;
				for(int i=0;i<N;i++) tmp += r[i][k][l];
				Rk[k] += tmp;
			}
		}
		double Tl[] = new double[L];
		// T_l
		for(int l=0;l<L;l++){
			Tl[l] = 0;
			for(int k=0;k<K;k++){
				double tmp = 0;
				for(int i=0;i<N;i++) tmp += r[i][k][l];
				Tl[l] += tmp;
			}
		}
		//Others
		double xk[] = new double[K];
		for(int k=0;k<param.m;k++){
			xk[k] = 0;
			for(int l=0;l<L;l++){
				double tmp = 0;
				for(int i=0;i<N;i++) tmp += r[i][k][l]*Math.pow(data[i]-initial_param.offset[l], 2);
				xk[k] += tmp;
			}
		}
		for(int k=param.m;k<K;k++){
			xk[k] = 0;
			for(int l=0;l<L;l++){
				double tmp = 0;
				for(int i=0;i<N;i++) tmp += r[i][k][l]*Math.abs(data[i]-initial_param.offset[l]);
				xk[k] += tmp;
			}
		}
		// Update the parameters.
		for(int l=0;l<param.offset.length;l++){
			param.p[l] = initial_param.p[l] + Tl[l];//(9)
		}
		for(int k=0;k<param.m;k++){
			param.alpha[k] = initial_param.alpha[k] + Rk[k]; //(11)
			param.a[k] = initial_param.a[k] + Rk[k]/2; //(12)
			param.b[k] = initial_param.b[k] + xk[k]/2; //(13)
		}
		for(int k=param.m;k<K;k++){
			param.alpha[k] = initial_param.alpha[k] + Rk[k]; //(11)
			param.a[k] = initial_param.a[k] + Rk[k]; //(12)
			param.b[k] = initial_param.b[k] + xk[k]; //(13)
		}
		log.info("Variational Bayesian M-step process finished.");
	}
	
	/**
	 * Lower bound
	 * @param initial_param Prior distribution(IN)
	 * @param r Values of `r'(IN)
	 * @param Rk \sum_{i=1}^N \sum_{l=1}^L r_{i,k,l}(IN)
	 * @param param Parameter distribution(IN)
	 * @return lower bound
	 */
	double lowerbound(ONDEParameterDistribution initial_param, ONDEParameterDistribution param, double r[][][], double Rk[]){
		log.info("Variational Bayesian lower bound calculation process started.");
		int N = r.length;
		int K = Rk.length;
		int L = initial_param.offset.length;
		double alpha_hat = 0;
		double alpha0_hat = 0;
		for(int k=0;k<K;k++){
			alpha0_hat += initial_param.alpha[k];
			alpha_hat += param.alpha[k];
		}
		double p_hat = 0;
		double p0_hat = 0;
		for(int k=0;k<L;k++){
			p0_hat += initial_param.p[k];
			p_hat += param.p[k];
		}
		double tmp = Gamma.logGamma(alpha0_hat) - Gamma.logGamma(alpha_hat) + Gamma.logGamma(p0_hat) - Gamma.logGamma(p_hat); //(12) independent of k and l
		// terms only  dependent on k
		for(int k=0;k<K;k++){
			tmp += - Gamma.logGamma(initial_param.alpha[k]) + Gamma.logGamma(param.alpha[k])
				+initial_param.a[k]*Math.log(initial_param.b[k]) - param.a[k]*Math.log(param.b[k])
				-Gamma.logGamma(initial_param.a[k]) + Gamma.logGamma(param.a[k])
				- ((k<initial_param.m) ? Math.log(2*Math.PI)*Rk[k]/2 : Math.log(2)*Rk[k]); //(12)
		}
		// terms only  dependent on l
		for(int l=0;l<L;l++){
			tmp += - Gamma.logGamma(initial_param.p[l]) + Gamma.logGamma(param.p[l]); //(9)
		}
		for(int k=0;k<K;k++){
			for(int l=0;l<L;l++){
				double tmp2 = 0;
				for(int i=0;i<N;i++){
					if(r[i][k][l]>0) tmp2 -= r[i][k][l]*Math.log(r[i][k][l]); 
				}
				tmp += tmp2;
			}
		}
		log.info("Variational Bayesian lower bound calculation process finished.");
		log.info("Variational lower bound: " + tmp);		
		return tmp;
	}
	
}
