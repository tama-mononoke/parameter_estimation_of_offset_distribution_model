/**
 * 
 */
package jp.go.enri.prml.vb;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.math.special.Gamma;

/**
 * 変分ベイズ法によるガウス分布とラプラス分布の混合分布のパラメータ推定。
 * アルゴリズムは「変分ベイズ法による正規分布と両側指数分布の混合分布のパラメータ推定」に記載。
 * @author 藤田雅人（電子航法研究所）
 * @version 1.0.1　(Last update: 30/11/2011)
 *
 */
public class NDEVB {
	/**
	 * ログの取得
	 */
	public static Log log = LogFactory.getLog(NDEVB.class);
	/**
	 * とても小さな正の数。変分下限の増加がこれ以下になったら、計算をやめる。
	 */
	public static double threshold = 10E-3; 
	
	/**
	 * 変分ベイズ法による計算値を格納
	 * @author 藤田雅人
	 *
	 */
	public static class Result{
		/**
		 * ガウス混合分布の変分ベイズ法パラメータを取得。
		 * @return ガウス混合分布の変分ベイズ法パラメータ
		 */
		public NDEParameterDistribution getParam() {
			return param;
		}
		/**
		 * 変分下限を取得
		 * @return 変分下限
		 */
		public double getLowerbound() {
			return lowerbound;
		}
		/**
		 * ガウス混合分布の変分ベイズ法パラメータ
		 */
		NDEParameterDistribution param;
		/**
		 * 変分下限
		 */
		double lowerbound;
	}
	
	/**
	 * 変分ベイズ法によるパラメータ推定
	 * @param initial_params 事前分布
	 * @param dataMatrix データ行列
	 * @return 推定結果
	 */
	public Result estimate(NDEParameterDistribution initial_params, double data[]){
		log.info("Variational Bayesian estimation process started.");
		// 変数の初期化
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
			// 変分下限の計算
			lowerbound_pre = lowerbound;
			lowerbound = lowerbound(initial_params,params,r,Nk);
			// 停止条件
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
	 * E-stepを実装
	 * @param param パラメータ分布(IN)
	 * @param data データ行列(IN)
	 * @param r rの計算値を格納する配列(OUT)
	 */
	void Estep(NDEParameterDistribution param, double data[], double r[][]){
		log.info("Variational Bayesian E-step process started.");
		int K = param.m + param.n;
		// 式(3)のiによらない部分をあらかじめ計算
		double tmp1[] = new double[K];
		for(int k=0;k<param.m;k++){
			tmp1[k] = Gamma.digamma(param.alpha[k])+0.5*(Gamma.digamma(param.a[k])-Math.log(param.b[k])-Math.log(2*Math.PI));
		}
		for(int k=param.m;k<K;k++){
			tmp1[k] = Gamma.digamma(param.alpha[k])+(Gamma.digamma(param.a[k])-Math.log(param.b[k])-Math.log(2));
		}
		int N = data.length;
		for(int i=0;i<N;i++){
			// rho_{nk}の計算値が小さくなりすぎるので、log(rho_{nk})の最大値の値を保持
			double rhoMax = Double.NEGATIVE_INFINITY;
			for(int k=0;k<param.m;k++){
				double tmp = tmp1[k]-0.5*param.a[k]*Math.pow(data[i],2)/param.b[k]; //式(3)
				r[i][k] = tmp;
				rhoMax = Math.max(rhoMax, tmp);
			}
			for(int k=param.m;k<K;k++){
				double tmp = tmp1[k]-param.a[k]*Math.abs(data[i])/param.b[k]; //式(3)
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
				r[i][j] /= sum; //式(4)
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
	 * M-stepを実装
	 * @param initial_param 事前分布(IN)
	 * @param data データ行列(IN)
	 * @param r rの計算値(IN)
	 * @param Nk \sum_{i=1}^N r_{i,k}の計算値(OUT)
	 * @param param パラメータ分布を格納する変数(OUT)
	 */
	void Mstep(NDEParameterDistribution initial_param, double data[], double r[][], double Nk[], NDEParameterDistribution param){
		log.info("Variational Bayesian M-step process started.");
		int N = data.length;
		int K = initial_param.m + initial_param.n;
		// \sum_{i=1}^N r_{i,k}を計算
		for(int i=0;i<K;i++){
			Nk[i] = 0;
			for(int j=0;j<N;j++) Nk[i] += r[j][i];
		}
		// \sum_{i=1}^N r_{i,k}x_i^2と\sum_{i=1}^N r_{i,k}|x_i|を計算
		double xk[] = new double[K];
		for(int k=0;k<param.m;k++){
			xk[k] = 0;
			for(int i=0;i<N;i++) xk[k] += r[i][k]*Math.pow(data[i], 2);
		}
		for(int k=param.m;k<K;k++){
			xk[k] = 0;
			for(int i=0;i<N;i++) xk[k] += r[i][k]*Math.abs(data[i]);
		}
		// パラメータの値の更新
		for(int k=0;k<param.m;k++){
			param.alpha[k] = initial_param.alpha[k] + Nk[k]; //式(6)
			param.a[k] = initial_param.a[k] + Nk[k]/2; //式(7)
			param.b[k] = initial_param.b[k] + xk[k]/2; //式(8)
		}
		for(int k=param.m;k<K;k++){
			param.alpha[k] = initial_param.alpha[k] + Nk[k]; //式(6)
			param.a[k] = initial_param.a[k] + Nk[k]; //式(7)
			param.b[k] = initial_param.b[k] + xk[k]; //式(8)
		}
		log.info("Variational Bayesian M-step process finished.");
	}
	
	/**
	 * 変分下限を計算。
	 * @param initial_param 事前分布(IN)
	 * @param r rの計算値(IN)
	 * @param Nk \sum_{i=1}^N r_{i,k}の計算値(IN)
	 * @param param パラメータ分布(IN)
	 * @return 変分下限
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
		double tmp = Gamma.logGamma(alpha0_hat) - Gamma.logGamma(alpha_hat); //式(9)のkによらない項
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
				- ((k<initial_param.m) ? Math.log(2*Math.PI)*Nk[k]/2 : Math.log(2)*Nk[k]); //式(9)
		}
		log.info("Variational Bayesian lower bound calculation process finished.");
		log.info("Variational lower bound: " + tmp);		
		return tmp;
	}
	
}
