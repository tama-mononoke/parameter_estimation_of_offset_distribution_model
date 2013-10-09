/**
 * 
 */
package jp.go.enri.prml.bf;


import jp.go.enri.prml.dist.NDE;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * 総当たり法によるガウス分布とラプラス分布の混合分布のパラメータ最尤推定。
 * @author 藤田雅人（電子航法研究所）
 * @version 1.0.1　(Last update: 30/11/2011)
 *
 */
public class NDEBF {
	/**
	 * ログの取得
	 */
	public static Log log = LogFactory.getLog(NDEBF.class);
	/**
	 * 対数尤度用閾値
	 */
	public static double logLiklihood_threshold = 10E-3;
	/**
	 * 総当たり法による計算値を格納
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
	 * パラメータの値を増加させるのに用いる。
	 * @author 藤田雅人
	 */
	static class Increment{
		/**
		 * ログの取得
		 */
		public static Log log = LogFactory.getLog(Increment.class);
		/**
		 * 分割数
		 */
		private int D;
		/**
		 * 混合分布に含まれるガウス分布の数
		 */
		private int m;
		/**
		 * 混合分布に含まれるラプラス分布の数
		 */
		private int n;
		/**
		 * alphaのとりうる範囲
		 */
		private double alpha_bound[][];
		/**
		 * sigmaのとりうる範囲
		 */
		private double sigma_bound[][];
		/**
		 * lambdaのとりうる範囲
		 */
		private double lambda_bound[][];
		/**
		 * ポインタ
		 */
		private int pointer[][];
		/**
		 * コンストラクタ
		 * @param d
		 * @param alpha_bound
		 * @param sigma_bound
		 * @param lambda_bound
		 */
		Increment(int d, double[][] alpha_bound, double[][] sigma_bound,
				double[][] lambda_bound) {
			super();
			D = d;
			this.alpha_bound = alpha_bound;
			this.sigma_bound = sigma_bound;
			this.lambda_bound = lambda_bound;
			m = sigma_bound.length;
			n = lambda_bound.length;
			pointer = new int[3][];
			pointer[0] = new int[m+n-1];
			pointer[1] = new int[m];
			pointer[2] = new int[n];
			for(int i=0;i<3;i++){
				for(int j=0;j<pointer[i].length;j++) pointer[i][j] = 0;
			}
		}
		/**
		 * ポインタの増加
		 * @return ポインタが終端に到達すればfalse。そうでなければtrue。
		 */
		boolean increment(){
			while(increment2()){
				double sum = 0;
				for(int i=0;i<alpha_bound.length;i++){
					sum += alpha_bound[i][0]*(1-(2*pointer[0][i]+1)/(2*D)) + alpha_bound[i][1]*(2*pointer[0][i]+1)/(2*D);
				}
				if(sum<=1){
					return true;
				}
			}
			return false;	
		}
		
		/**
		 * ポインタの増加
		 * @return ポインタが終端に到達すればfalse。そうでなければtrue。
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
		 * ポインタに対応するパラメータ値の取得
		 * @return パラメータ値
		 */
		NDE getCurrent(){
			double alpha[] = new double[m+n];
			double sum = 0;
			for(int i=0;i<m+n-1;i++){
				alpha[i] = alpha_bound[i][0]*(1-(2.0*pointer[0][i]+1)/(2*D)) + alpha_bound[i][1]*(2.0*pointer[0][i]+1)/(2*D);
				sum += alpha[i];
			}
			alpha[m+n-1] = 1-sum;
			double sigma[] = new double[m];
			for(int i=0;i<m;i++){
				sigma[i] = sigma_bound[i][0]*(1-(2.0*pointer[1][i]+1)/(2*D)) + sigma_bound[i][1]*(2.0*pointer[1][i]+1)/(2*D);
			}
			double lambda[] = new double[n];
			for(int i=0;i<n;i++){
				lambda[i] = lambda_bound[i][0]*(1-(2.0*pointer[2][i]+1)/(2*D)) + lambda_bound[i][1]*(2.0*pointer[2][i]+1)/(2*D);
			}
			return new NDE(alpha,sigma,lambda);
		}
		
		Increment getNewIncrement(){
			double tmp_alpha_bound[][] = new double[alpha_bound.length][];
			double tmp_sigma_bound[][] = new double[sigma_bound.length][];
			double tmp_lambda_bound[][] = new double[lambda_bound.length][];
			for(int i=0;i<alpha_bound.length;i++){
				tmp_alpha_bound[i] = new double[2];
				tmp_alpha_bound[i][0] = alpha_bound[i][0]*(1- (double) pointer[0][i]/D) + alpha_bound[i][1]*((double) pointer[0][i]/D);
				tmp_alpha_bound[i][1] = alpha_bound[i][0]*(1- (double) (pointer[0][i]+1)/D) + alpha_bound[i][1]*((double) (pointer[0][i]+1)/D);
			}
			for(int i=0;i<sigma_bound.length;i++){
				tmp_sigma_bound[i] = new double[2];
				tmp_sigma_bound[i][0] = sigma_bound[i][0]*(1- (double) pointer[1][i]/D) + sigma_bound[i][1]*((double) pointer[1][i]/D);
				tmp_sigma_bound[i][1] = sigma_bound[i][0]*(1- (double) (pointer[1][i]+1)/D) + sigma_bound[i][1]*((double) (pointer[1][i]+1)/D);
			}
			for(int i=0;i<lambda_bound.length;i++){
				tmp_lambda_bound[i] = new double[2];
				tmp_lambda_bound[i][0] = lambda_bound[i][0]*(1- (double) pointer[2][i]/D) + lambda_bound[i][1]*((double) pointer[2][i]/D);
				tmp_lambda_bound[i][1] = lambda_bound[i][0]*(1- (double) (pointer[2][i]+1)/D) + lambda_bound[i][1]*((double) (pointer[2][i]+1)/D);
			}
			return new Increment(D, tmp_alpha_bound, tmp_sigma_bound, tmp_lambda_bound);
		}
		
		void print(){
			if(log.isInfoEnabled()){
				StringBuilder sb = new StringBuilder();
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
	 * 総当たり法によるパラメータ推定
	 * @param D 
	 * @param initial_params 初期パラメータ
	 * @param data 観測値
	 * @return 推定結果
	 */
	public Result estimate(NDEParameter initial_params, double[] data){
		log.info("BF method estimation process started.");
		// 変数の初期化
		log.info("BF method initialization process started.");
		int N = data.length;
		double initial_alpha_bound[][] = new double[initial_params.initial_sigma_bound.length+initial_params.initial_lambda_bound.length-1][];
		for(int i=0;i<initial_alpha_bound.length;i++){
			initial_alpha_bound[i] = new double[2];
			for(int j=0;j<2;j++) initial_alpha_bound[i][j] = j;
		}
		double logLiklihood_pre = Double.NEGATIVE_INFINITY;
		Increment incre = new Increment(initial_params.D, initial_alpha_bound, initial_params.initial_sigma_bound, initial_params.initial_lambda_bound);
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
			// 停止条件
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
	 * 対数尤度を計算
	 * @param param ガウス分布とラプラス分布の混合分布のパラメータ
	 * @param data 観測値
	 * @return 対数尤度
	 */
	double logLiklihood(NDE param,double data[]){
		//log.info("BF method log-liklihood evaluation process started.");
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
		//log.info("BF method log-liklihood evaluation process finished.");
		return tmp;
	}

}
