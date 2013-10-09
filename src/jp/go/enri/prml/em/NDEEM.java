/**
 * 
 */
package jp.go.enri.prml.em;


import jp.go.enri.prml.dist.NDE;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * EM法によるガウス分布とラプラス分布の混合分布のパラメータ推定。
 * @author 藤田雅人（電子航法研究所）
 * @version 1.0.1　(Last update: 30/11/2011)
 *
 */
public class NDEEM {
	/**
	 * ログの取得
	 */
	public static Log log = LogFactory.getLog(NDEEM.class);
	/**
	 * 対数尤度用閾値
	 */
	public static double logLiklihood_threshold = 10E-3;
	/**
	 * EM法による計算値を格納
	 * @author 藤田雅人
	 *
	 */
	public static class Result{
		/**
		 * ガウス分布とラプラス分布の混合分布のパラメータを取得。
		 * @return ガウス分布とラプラス分布の混合分布のパラメータ
		 */
		public NDE getParam() {
			return param;
		}
		/**
		 * 対数尤度を取得
		 * @return 対数尤度
		 */
		public double getLogLiklihood() {
			return logLiklihood;
		}
		/**
		 * ガウス分布とラプラス分布の混合分布のパラメータ
		 */
		NDE param;
		/**
		 * 対数尤度
		 */
		double logLiklihood;
	}
	
	/**
	 * EM法によるパラメータ推定
	 * @param initial_params 初期パラメータ
	 * @param data 観測値
	 * @return 推定結果
	 */
	public Result estimate(NDE initial_params, double[] data){
		log.info("EM method estimation process started.");
		// 変数の初期化
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
			// 対数尤度の計算
			logLiklihood_pre = logLiklihood;
			logLiklihood = logLiklihood(params,data);
			// 停止条件
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
	 * E-stepを実装。負担率を計算。
	 * @param param 混合正規分布のパラメータ
	 * @param data 観測値
	 * @param gamma 負担率の計算値を格納する配列
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
	 * M-stepを実装。
	 * @param m 混合分布に含まれるガウス分布の数
	 * @param n 混合分布に含まれるラプラス分布の数
	 * @param data 観測値
	 * @param gamma 負担率の計算値
	 * @return ガウス分布とラプラス分布の混合分布のパラメータ
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
	 * 対数尤度を計算
	 * @param param ガウス分布とラプラス分布の混合分布のパラメータ
	 * @param data 観測値
	 * @return 対数尤度
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
