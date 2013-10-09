/**
 * 
 */
package jp.go.enri.prml.vb;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.math.special.Gamma;

/**
 * 変分ベイズ法によるガウス分布とラプラス分布の混合オフセット分布のパラメータ推定。
 * アルゴリズムは「変分ベイズ法による正規分布と両側指数分布の混合オフセット分布のパラメータ推定」に記載。
 * @author 藤田雅人（電子航法研究所）
 * @version 1.0.1　(Last update: 06/12/2011)
 *
 */
public class ONDEVB {
	/**
	 * ログの取得
	 */
	public static Log log = LogFactory.getLog(ONDEVB.class);
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
		 * ガウス混合オフセット分布の変分ベイズ法パラメータを取得。
		 * @return ガウス混合オフセット分布の変分ベイズ法パラメータ
		 */
		public ONDEParameterDistribution getParam() {
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
		 * ガウス混合オフセット分布の変分ベイズ法パラメータ
		 */
		ONDEParameterDistribution param;
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
	public Result estimate(ONDEParameterDistribution initial_params, double data[]){
		log.info("Variational Bayesian estimation process started.");
		// 変数の初期化
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
			// 変分下限の計算
			lowerbound_pre = lowerbound;
			lowerbound = lowerbound(initial_params,params,r,Rk);
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
	
	void Estep(ONDEParameterDistribution param, double data[], double r[][][]){
		log.info("Variational Bayesian E-step process started.");
		int K = param.m + param.n;
		int L = param.offset.length;
		double tmp_t[] = new double[L];
		double tmp_r[] = new double[K];
		// 式(4)の値を計算
		// lのみに依存する項をあらかじめ計算
		for(int l=0;l<L;l++){
			tmp_t[l] = Gamma.digamma(param.p[l]); //式(4)
		}
		// kのみに依存する項をあらかじめ計算
		for(int k=0;k<param.m;k++){
			tmp_r[k] = Gamma.digamma(param.alpha[k]) + 0.5*(Gamma.digamma(param.a[k])-Math.log(param.b[k])-Math.log(2*Math.PI)); //式(4)
		}
		for(int k=param.m;k<K;k++){
			tmp_r[k] = Gamma.digamma(param.alpha[k]) + (Gamma.digamma(param.a[k])-Math.log(param.b[k])-Math.log(2)); //式(4)
		}
		
		int N = data.length;
		for(int i=0;i<N;i++){
			// rho_{ikl}の計算値が小さくなりすぎるので、log(rho_{ikl})の最大値の値を保持
			double rhoMax = Double.NEGATIVE_INFINITY;
			for(int l=0;l<L;l++){
				for(int k=0;k<param.m;k++){
					r[i][k][l] = tmp_t[l]+tmp_r[k]-0.5*param.a[k]*Math.pow(data[i]-param.offset[l],2)/param.b[k]; //式(4)
					rhoMax = Math.max(rhoMax, r[i][k][l]);
				}
				for(int k=param.m;k<K;k++){
					r[i][k][l] = tmp_t[l]+tmp_r[k]-param.a[k]*Math.abs(data[i]-param.offset[l])/param.b[k]; //式(4)
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
	 * M-stepを実装
	 * @param initial_param 事前分布(IN)
	 * @param data データ行列(IN)
	 * @param r rの計算値(IN)
	 * @param Rk \sum_{i=1}^N \sum_{l=1}^L r_{i,k,l}の計算値(OUT)
	 * @param param パラメータ分布を格納する変数(OUT)
	 */
	void Mstep(ONDEParameterDistribution initial_param, double data[], double r[][][], double Rk[], ONDEParameterDistribution param){
		log.info("Variational Bayesian M-step process started.");
		int N = data.length;
		int K = initial_param.m + initial_param.n;
		int L = initial_param.offset.length;
		// R_kを計算
		for(int k=0;k<K;k++){
			Rk[k] = 0;
			for(int l=0;l<L;l++){
				double tmp = 0;
				for(int i=0;i<N;i++) tmp += r[i][k][l];
				Rk[k] += tmp;
			}
		}
		double Tl[] = new double[L];
		// T_lを計算
		for(int l=0;l<L;l++){
			Tl[l] = 0;
			for(int k=0;k<K;k++){
				double tmp = 0;
				for(int i=0;i<N;i++) tmp += r[i][k][l];
				Tl[l] += tmp;
			}
		}
		//その他の項を計算
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
		// パラメータの値の更新
		for(int l=0;l<param.offset.length;l++){
			param.p[l] = initial_param.p[l] + Tl[l];//式(9)
		}
		for(int k=0;k<param.m;k++){
			param.alpha[k] = initial_param.alpha[k] + Rk[k]; //式(11)
			param.a[k] = initial_param.a[k] + Rk[k]/2; //式(12)
			param.b[k] = initial_param.b[k] + xk[k]/2; //式(13)
		}
		for(int k=param.m;k<K;k++){
			param.alpha[k] = initial_param.alpha[k] + Rk[k]; //式(11)
			param.a[k] = initial_param.a[k] + Rk[k]; //式(12)
			param.b[k] = initial_param.b[k] + xk[k]; //式(13)
		}
		log.info("Variational Bayesian M-step process finished.");
	}
	
	/**
	 * 変分下限を計算。
	 * @param initial_param 事前分布(IN)
	 * @param r rの計算値(IN)
	 * @param t tの計算値(IN)
	 * @param Rk \sum_{i=1}^N \sum_{l=1}^L r_{i,k,l}の計算値(IN)
	 * @param param パラメータ分布(IN)
	 * @return 変分下限
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
		double tmp = Gamma.logGamma(alpha0_hat) - Gamma.logGamma(alpha_hat) + Gamma.logGamma(p0_hat) - Gamma.logGamma(p_hat); //式(12)のlやkによらない項
		// kのみに依存する項
		for(int k=0;k<K;k++){
			tmp += - Gamma.logGamma(initial_param.alpha[k]) + Gamma.logGamma(param.alpha[k])
				+initial_param.a[k]*Math.log(initial_param.b[k]) - param.a[k]*Math.log(param.b[k])
				-Gamma.logGamma(initial_param.a[k]) + Gamma.logGamma(param.a[k])
				- ((k<initial_param.m) ? Math.log(2*Math.PI)*Rk[k]/2 : Math.log(2)*Rk[k]); //式(12)
		}
		// lのみに依存する項
		for(int l=0;l<L;l++){
			tmp += - Gamma.logGamma(initial_param.p[l]) + Gamma.logGamma(param.p[l]); //式(9)
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
