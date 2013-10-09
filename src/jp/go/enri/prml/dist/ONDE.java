/**
 * 
 */
package jp.go.enri.prml.dist;

import java.util.Arrays;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * ガウス分布とラプラス分布の混合オフセット分布のパラメータ
 * @author 藤田雅人（電子航法研究所）
 * @version 1.0.1　(Last update: 06/12/2011)
 *
 */
public class ONDE {
	/**
	 * オフセット値を取得
	 * @return オフセット値
	 */
	public double[] getOffset() {
		return offset;
	}
	/**
	 * オフセット混合比を取得
	 * @return オフセット混合比
	 */
	public double[] getOmega() {
		return omega;
	}

	/**
	 * ログの取得
	 */
	public static Log log = LogFactory.getLog(ONDE.class);
	/**
	 * 混合分布に含まれるガウス分布の数を取得
	 * @return 混合分布に含まれるガウス分布の数
	 */
	public int getM() {
		return m;
	}
	/**
	 * 混合分布に含まれるラプラス分布の数を取得
	 * @return 混合分布に含まれるラプラス分布の数
	 */
	public int getN() {
		return n;
	}
	/**
	 * 混合比の値を取得
	 * @return 混合比の値
	 */
	public double[] getPi() {
		return Arrays.copyOf(pi, pi.length);
	}
	/**
	 * ラプラス分布のスケールパラメータの値
	 * @return ラプラス分布のスケールパラメータの値を取得
	 */
	public double[] getLambda() {
		return Arrays.copyOf(lambda, lambda.length);
	}
	/**
	 * ガウス分布の標準偏差の値を取得
	 * @return sigma ガウス分布の標準偏差の値
	 */
	public double[] getSigma() {
		return Arrays.copyOf(sigma, sigma.length);
	}
	/**
	 * とりうるオフセットの数を返す。
	 * @return とりうるオフセットの数
	 */
	public int getL(){
		return offset.length;
	}
	/**
	 * 混合分布に含まれるガウス分布の数
	 */
	int m;
	/**
	 * 混合分布に含まれるラプラス分布の数
	 */
	int n;
	/**
	 * 混合比の値
	 */
	double pi[];
	/**
	 * ラプラス分布のスケールパラメータの値
	 */
	double lambda[];
	/**
	 * ガウス分布の標準偏差の値
	 */
	double sigma[];
	/**
	 * オフセット値
	 */
	double offset[];
	/**
	 * オフセット混合比
	 */
	double omega[];
	/**
	 * コンストラクタ
	 */
	private ONDE(){}
	/**
	 * コンストラクタ
	 * @param pi 混合比の値
	 * @param sigma ガウス分布の標準偏差の値
	 * @param lambda ラプラス分布のスケールパラメータの値
	 * @param offset オフセット値
	 * @param omega オフセット混合比
	 */
	public ONDE(double[] pi, double[] sigma, double[] lambda, double[] offset, double[] omega) {
		super();
		if(sigma==null){
			log.error("The length of array sigma is not correct.");
			throw new IllegalArgumentException();
		}
		this.m = sigma.length;
		if(lambda==null){
			log.error("The length of array lambda is not correct.");
			throw new IllegalArgumentException();
		}
		this.n = lambda.length;
		int K = m+n;
		if(pi==null || pi.length!=K){
			log.error("The length of array pi is not correct.");
			throw new IllegalArgumentException();
		}
		for(int i=0;i<pi.length;i++){
			if(pi[i]<0){
				log.error("The value pi[" + i + "] should not be negative. pi[" + i + "]=" + pi[i]);
				throw new IllegalArgumentException();
			}
		}
		this.pi = Arrays.copyOf(pi, pi.length);
		for(int i=0;i<sigma.length;i++){
			if(sigma[i]<=0){
				log.error("The value sigma[" + i + "] should be positive. sigma[" + i + "]=" + sigma[i]);
				throw new IllegalArgumentException();
			}
		}
		this.sigma = Arrays.copyOf(sigma, sigma.length);
		for(int i=0;i<lambda.length;i++){
			if(lambda[i]<=0){
				log.error("The value lambda[" + i + "] should be positive. lambda[" + i + "]=" + lambda[i]);
				throw new IllegalArgumentException();
			}
		}
		this.lambda = Arrays.copyOf(lambda, lambda.length);
		if(offset==null || offset.length<1){
			log.error("The length of array offfset is not correct.");
			throw new IllegalArgumentException();
		}
		this.offset = Arrays.copyOf(offset,offset.length);
		if(omega==null || omega.length!=offset.length){
			log.error("The length of array omega is not correct.");
			throw new IllegalArgumentException();
		}
		this.omega = Arrays.copyOf(omega, omega.length);
	}
	
	/**
	 * @see java.lang.Object#clone()
	 */
	@Override
	public ONDE clone(){
		ONDE tmp = new ONDE();
		tmp.m = m;
		tmp.n = n;
		tmp.pi = Arrays.copyOf(pi, pi.length);
		tmp.sigma = Arrays.copyOf(sigma, sigma.length);
		tmp.lambda = Arrays.copyOf(lambda, lambda.length);
		tmp.offset = Arrays.copyOf(offset,offset.length);
		tmp.omega = Arrays.copyOf(omega, omega.length);
		return tmp;
	}
	/**
	 * 各オフセット値が適用されている確率を推定
	 * @param x 観測値
	 * @return 各オフセット値が適用されている確率
	 */
	public double[] estimateOffsetProbability(double x){
		double ans[] = new double[offset.length];
		double sum = 0;
		for(int l=0;l<offset.length;l++){
			double tmp2 = 0;
			for(int k=0;k<m;k++){
				tmp2 += pi[k]*Math.exp(-0.5*Math.pow(x-offset[l], 2)/(sigma[k]*sigma[k]))/(Math.sqrt(2*Math.PI)*sigma[k]);
			}
			for(int k=0;k<n;k++){
				tmp2 += pi[k+m]*Math.exp(-Math.abs(x-offset[l])/lambda[k])/(2*lambda[k]);
			}
			tmp2 *= omega[l];
			ans[l] = tmp2;
			sum += tmp2;
		}
		for(int l=0;l<offset.length;l++) ans[l] /= sum;
		return ans;
	}
	
}
