/**
 * 
 */
package jp.go.enri.prml.em;


import jp.go.enri.prml.dist.ONDE;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * EM法によるガウス分布とラプラス分布の混合オフセット分布のパラメータ推定。
 * @author 藤田雅人（電子航法研究所）
 * @version 1.0.1　(Last update: 06/12/2011)
 *
 */
public class ONDEEM {
	/**
	 * ログの取得
	 */
	public static Log log = LogFactory.getLog(ONDEEM.class);
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
		 * ガウス分布とラプラス分布の混合オフセット分布のパラメータを取得。
		 * @return ガウス分布とラプラス分布の混合オフセット分布のパラメータ
		 */
		public ONDE getParam() {
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
		 * ガウス分布とラプラス分布の混合オフセット分布のパラメータ
		 */
		ONDE param;
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
	public Result estimate(ONDE initial_params, double[] data){
		log.info("EM method estimation process started.");
		// 変数の初期化
		log.info("EM method initialization process started.");
		int N = data.length;
		int K = initial_params.getM() + initial_params.getN();
		int L = initial_params.getL();
		double gamma[][][] = new double[N][][];
		for(int i=0;i<gamma.length;i++){
			gamma[i] = new double[K][];
			for(int k=0;k<K;k++){
				gamma[i][k] = new double[L];
			}
		}
		double logLiklihood = Double.NEGATIVE_INFINITY;
		double logLiklihood_pre = Double.NEGATIVE_INFINITY;
		ONDE params = initial_params.clone();
		log.info("EM method initialization process finished.");
		boolean flag = true;
		do{
			// E step
			Estep(params,data,gamma);
			// M step
			params = Mstep(initial_params.getM(),initial_params.getN(),initial_params.getOffset(),data,gamma);
			// 対数尤度の計算
			logLiklihood_pre = logLiklihood;
			logLiklihood = logLiklihood(params,data);
//			System.out.println("Result");
//			for(int i=0;i<params.getL();i++){
//				System.out.println("offset[" + i + "]=" + params.getOffset()[i]);
//				System.out.println("omega[" + i + "]=" + params.getOmega()[i]);
//			}
//			for(int i=0;i<params.getM();i++){
//				System.out.println("pi[" + i + "]=" + params.getPi()[i]);
//				System.out.println("sigma[" + i + "]=" + params.getSigma()[i]);
//			}
//			for(int i=0;i<params.getN();i++){
//				System.out.println("pi[" + (i+params.getM()) + "]=" + params.getPi()[i+params.getM()]);
//				System.out.println("lambda[" + i + "]=" + params.getLambda()[i]);
//			}
//			System.out.println("-----------------------------");
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
	void Estep(ONDE param, double[] data, double gamma[][][]){
		log.info("EM method E-step process started.");
		int N = data.length;
		int K = param.getM() + param.getN();
		int L = param.getL();
		for(int i=0;i<N;i++){
			double sum = 0;
			double x = data[i];
			for(int k=0;k<param.getM();k++){
				for(int l=0;l<L;l++){
					gamma[i][k][l] = param.getOmega()[l]*param.getPi()[k]*Math.exp(-0.5*Math.pow(x-param.getOffset()[l], 2)/Math.pow(param.getSigma()[k],2))/(Math.sqrt(2*Math.PI)*param.getSigma()[k]);
					sum += gamma[i][k][l];
				}
			}
			for(int k=0;k<param.getN();k++){
				for(int l=0;l<L;l++){
					gamma[i][k+param.getM()][l] = param.getOmega()[l]*param.getPi()[k+param.getM()]*Math.exp(-Math.abs(x-param.getOffset()[l])/param.getLambda()[k])/(2*param.getLambda()[k]);
					sum += gamma[i][k+param.getM()][l];
				}
			}
			for(int k=0;k<K;k++){
				for(int l =0;l<L;l++) gamma[i][k][l] /= sum;
			}
		}
		log.info("EM method E-step process finished.");
	}
	
	/**
	 * M-stepを実装。
	 * @param m 混合分布に含まれるガウス分布の数
	 * @param n 混合分布に含まれるラプラス分布の数
	 * @param offset オフセット値
	 * @param data 観測値
	 * @param gamma 負担率の計算値
	 * @return ガウス分布とラプラス分布の混合分布のパラメータ
	 */
	ONDE Mstep(int m, int n, double offset[], double data[],double gamma[][][]){
		log.info("EM method M-step process started.");
		int N = data.length;
		int K = m+n;
		int L = offset.length;
		double omega[] = new double[L];
		for(int l=0;l<L;l++){
			for(int i=0;i<N;i++){
				double tmp_omega = 0;
				for(int k=0;k<K;k++) tmp_omega += gamma[i][k][l];
				omega[l] += tmp_omega;
			}
			omega[l] /= N;
		}
		double Nk[] = new double[K];
		double pi[] = new double[K];
		for(int i=0;i<K;i++){
			Nk[i] = 0;
			for(int j=0;j<N;j++){
				for(int l=0;l<L;l++) Nk[i] += gamma[j][i][l];
			}
			pi[i] = Nk[i]/N;
		}
		double sigma[] = new double[m];
		for(int k=0;k<m;k++){
			sigma[k] = 0;
			for(int i=0;i<N;i++){
				for(int l=0;l<L;l++){
					sigma[k] += gamma[i][k][l]*Math.pow(data[i]-offset[l], 2);
				}
			}
			sigma[k] /= Nk[k];
			sigma[k] = Math.sqrt(sigma[k]);
		}
		double lambda[] = new double[n];
		for(int k=0;k<n;k++){
			lambda[k] = 0;
			for(int i=0;i<N;i++){
				for(int l=0;l<L;l++){
					lambda[k] += gamma[i][m+k][l]*Math.abs(data[i]-offset[l]);
				}
			}
			lambda[k] /= Nk[m+k];
		}
		log.info("EM method M-step process finished.");
		return new ONDE(pi,sigma,lambda,offset,omega);
	}
	/**
	 * 対数尤度を計算
	 * @param param ガウス分布とラプラス分布の混合分布のパラメータ
	 * @param data 観測値
	 * @return 対数尤度
	 */
	double logLiklihood(ONDE param,double data[]){
		log.info("EM method log-liklihood evaluation process started.");
		int N = data.length;
		int L = param.getL();
		double tmp = 0;
		for(int i=0;i<N;i++){
			double tmp2 = 0;
			double x = data[i];
			for(int k=0;k<param.getM();k++){
				for(int l=0;l<L;l++){
					tmp2 += param.getOmega()[l]*param.getPi()[k]*Math.exp(-0.5*Math.pow(x-param.getOffset()[l], 2)/(param.getSigma()[k]*param.getSigma()[k]))/(Math.sqrt(2*Math.PI)*param.getSigma()[k]);
				}
			}
			for(int k=0;k<param.getN();k++){
				for(int l=0;l<L;l++){
					tmp2 += param.getOmega()[l]*param.getPi()[k+param.getM()]*Math.exp(-Math.abs(x-param.getOffset()[l])/param.getLambda()[k])/(2*param.getLambda()[k]);
				}
			}
			tmp += Math.log(tmp2);
		}
		log.info("EM method log-liklihood evaluation process finished.");
		return tmp;
	}

}
