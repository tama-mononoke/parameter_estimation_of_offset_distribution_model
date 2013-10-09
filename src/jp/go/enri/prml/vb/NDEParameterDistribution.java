/**
 * 
 */
package jp.go.enri.prml.vb;

import java.util.Arrays;

import jp.go.enri.prml.dist.Dirichlet;
import jp.go.enri.prml.dist.Gamma;
import jp.go.enri.prml.dist.NDE;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * ガウス分布とラプラス分布の混合分布のパラメータの事前／事後分布のパラメータ
 * @author 藤田雅人（電子航法研究所）
 * @version 1.0.1　(Last update: 30/11/2011)
 */
public class NDEParameterDistribution {
	/**
	 * 混合分布に含まれるガウス分布の数の取得
	 * @return 混合分布に含まれるガウス分布の数
	 */
	public int getM() {
		return m;
	}

	/**
	 * 混合分布に含まれるラプラス分布の数の取得
	 * @return 混合分布に含まれるラプラス分布の数
	 */
	public int getN() {
		return n;
	}

	/**
	 * 混合比\piの事前／事後分布（ディリクレ分布）のパラメータの取得
	 * @return 混合比\piの事前／事後分布（ディリクレ分布）のパラメータ
	 */
	public double[] getAlpha() {
		return Arrays.copyOf(alpha, alpha.length);
	}

	/**
	 * 精度パラメータ\etaの事前／事後分布（ガンマ分布）のパラメータaの取得
	 * @return 精度パラメータ\etaの事前／事後分布（ガンマ分布）のパラメータa
	 */
	public double[] getA() {
		return Arrays.copyOf(a, a.length);
	}

	/**
	 * 精度パラメータ\etaの事前／事後分布（ガンマ分布）のパラメータbの取得
	 * @return 精度パラメータ\etaの事前／事後分布（ガンマ分布）のパラメータb
	 */
	public double[] getB() {
		return Arrays.copyOf(b, b.length);
	}

	/**
	 * ログの取得
	 */
	public static Log log = LogFactory.getLog(NDEParameterDistribution.class);
	/**
	 * 混合分布に含まれるガウス分布の数
	 */
	int m;
	/**
	 * 混合分布に含まれるラプラス分布の数
	 */
	int n;
	/**
	 * 混合比\piの事前／事後分布（ディリクレ分布）のパラメータ
	 */
	double alpha[];
	/**
	 * 精度パラメータ\etaの事前／事後分布（ガンマ分布）のパラメータa
	 */
	double a[];
	/**
	 * 精度パラメータ\etaの事前／事後分布（ガンマ分布）のパラメータb
	 */
	double b[];
	/**
	 * 乱数生成用ガンマ分布
	 */
	private Gamma gamma[];
	/**
	 * 乱数生成用ディリクレ分布
	 */
	private Dirichlet dirichlet;
	/**
	 * @see java.lang.Object#clone()
	 */
	@Override
	public NDEParameterDistribution clone(){
		double[] alpha_tmp = Arrays.copyOf(alpha, alpha.length);
		double[] a_tmp = Arrays.copyOf(a, a.length);
		double[] b_tmp = Arrays.copyOf(b, b.length);
		return new NDEParameterDistribution(m, n, alpha_tmp, a_tmp, b_tmp);
	}
	
	/**
	 * コンストラクタ
	 * @param m 混合分布に含まれるガウス分布の数
	 * @param n 混合分布に含まれるラプラス分布の数
	 * @param alpha 混合比\piの事前／事後分布（ディリクレ分布）のパラメータ ((m+n)次元配列）
	 * @param a 精度パラメータ\etaの事前／事後分布（ガンマ分布）のパラメータa ((m+n)次元配列）
	 * @param b 精度パラメータ\etaの事前／事後分布（ガンマ分布）のパラメータb ((m+n)次元配列）
	 */
	public NDEParameterDistribution(int m, int n, double[] alpha, double[] a, double[] b){
		initialize(m, n, alpha, a, b);
	}
	
	/**
	 * コンストラクタ
	 * @param m 混合分布に含まれるガウス分布の数
	 * @param n 混合分布に含まれるラプラス分布の数
	 * @param alpha0 混合比\piの事前／事後分布（ディリクレ分布）のパラメータ (共有）
	 * @param a 精度パラメータ\etaの事前／事後分布（ガンマ分布）のパラメータa ((m+n)次元配列）
	 * @param b 精度パラメータ\etaの事前／事後分布（ガンマ分布）のパラメータb ((m+n)次元配列）
	 */
	public NDEParameterDistribution(int m, int n, double alpha0, double[] a, double[] b){
		double tmpAlpha[] = new double[m+n];
		Arrays.fill(tmpAlpha, alpha0);
		initialize(m, n, tmpAlpha, a, b);
	}
	

	/**
	 * 初期化
	 * @param m 混合分布に含まれるガウス分布の数
	 * @param n 混合分布に含まれるラプラス分布の数
	 * @param alpha 混合比\piの事前／事後分布（ディリクレ分布）のパラメータ ((m+n)次元配列）
	 * @param a 精度パラメータ\etaの事前／事後分布（ガンマ分布）のパラメータa ((m+n)次元配列）
	 * @param b 精度パラメータ\etaの事前／事後分布（ガンマ分布）のパラメータb ((m+n)次元配列）
	 */
	private void initialize(int m, int n, double[] alpha, double[] a, double[] b){
		if(m<0 || n<0){
			log.error("The number m and n should not be negative.");
			throw new IllegalArgumentException();
		}
		this.m = m;
		this.n = n;
		int K = m+n;
		if(alpha==null || alpha.length != K || a==null || a.length != K || b==null || b.length != K){
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
		this.alpha = Arrays.copyOf(alpha, alpha.length);
		this.a = Arrays.copyOf(a, a.length);
		this.b = Arrays.copyOf(b, b.length);
		// 他のパラメータの更新
		updateOtherParameters();
	}
	/**
	 * 他のパラメータの更新
	 */
	void updateOtherParameters(){
		gamma = new Gamma[m+n];
		for(int i=0;i<m+n;i++) gamma[i] = new Gamma(a[i], b[i]);
		dirichlet = new Dirichlet(alpha);
	}
	
	/**
	 * パラメータの値をMAP推定
	 * @return MAP推定値
	 */
	public NDE MAPEstimation(){
		log.info("Variational Bayesian MAP estimation process started.");
		int K = m+n;
		double pi[] = new double[K];
		double sum = 0;
		for(int i=0;i<K;i++) sum += alpha[i];
		for(int i=0;i<K;i++) pi[i] = (alpha[i]-1)/(sum-K);
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
		NDE para = new NDE(pi,sigma,lambda);
		log.info("Variational Bayesian MAP estimation process finished.");
		if(log.isInfoEnabled()){
			StringBuilder sb = new StringBuilder();
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
	 * パラメータの値を平均推定
	 * @return 平均推定値
	 */
	public NDE ptEstimation(){
		log.info("Variational Bayesian point estimation process started.");
		int K = m+n;
		double pi[] = new double[K];
		double sum = 0;
		for(int i=0;i<K;i++) sum += alpha[i];
		for(int i=0;i<K;i++) pi[i] = alpha[i]/sum;
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
		NDE para = new NDE(pi,sigma,lambda);
		log.info("Variational Bayesian point estimation process finished.");
		if(log.isInfoEnabled()){
			StringBuilder sb = new StringBuilder();
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
	 * パラメータ分布からサンプリング
	 * @return
	 */
	public NDE sampling(){
		double sigma[] = new double[m];
		double lambda[] = new double[n];
		for(int i=0;i<m;i++) sigma[i] = 1/Math.sqrt(gamma[i].nextDouble());
		for(int i=0;i<n;i++) lambda[i] = 1/gamma[i+m].nextDouble();
		return new NDE(dirichlet.nextDouble(), sigma, lambda);
	}
	
}
