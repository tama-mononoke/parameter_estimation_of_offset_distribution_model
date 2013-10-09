/**
 * 
 */
package jp.go.enri.prml.dist;

import java.util.Arrays;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.math.MathException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.special.Erf;

/**
 * ガウス分布とラプラス分布の混合分布のパラメータ
 * @author 藤田雅人（電子航法研究所）
 * @version 1.0.1　(Last update: 30/11/2011)
 *
 */
public class NDE {
	/**
	 * ログの取得
	 */
	public static Log log = LogFactory.getLog(NDE.class);
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
	 * コンストラクタ
	 */
	private NDE(){}
	/**
	 * コンストラクタ
	 * @param pi 混合比の値
	 * @param sigma ガウス分布の標準偏差の値
	 * @param lambda ラプラス分布のスケールパラメータの値
	 */
	public NDE(double[] pi, double[] sigma, double[] lambda) {
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
	}
	
	/**
	 * @see java.lang.Object#clone()
	 */
	@Override
	public NDE clone(){
		NDE tmp = new NDE();
		tmp.m = m;
		tmp.n = n;
		tmp.pi = Arrays.copyOf(pi, pi.length);
		tmp.sigma = Arrays.copyOf(sigma, sigma.length);
		tmp.lambda = Arrays.copyOf(lambda, lambda.length);
		return tmp;
	}
	/**
	 * -x以上x以下になる確率を計算
	 * @param x x>0となる値
	 * @return 確率
	 * @throws MathException
	 */
	public double cumulate(double x) throws MathException{
		if(x<0) throw new IllegalArgumentException();
		double sum = 0;
		for(int i=0;i<m;i++){
			double val = 1;
			double tmp = x/(Math.sqrt(2)*sigma[i]);
			if(tmp < 20){ //1-Erf(20)=5.40E-176でほとんど1と見なして構わない。
				try{
					val = Erf.erf(tmp);
				}
				catch(MaxIterationsExceededException ex){
					log.info("Failed to calculate Erf.erf(" + tmp + ")", ex);
					val = 1;
				}
			}
			sum += pi[i]*val;
		}
		for(int i=0;i<n;i++){
			sum += pi[i+m]*(1-Math.exp(-x/lambda[i]));
		}
		return sum;
	}
	
	private static final double incre = 1;
	/**
	 * -x以上x以下になる確率がgammaと等しくなるxを計算
	 * @param gamma
	 * @param epsilon 許容誤差
	 * @return
	 * @throws MathException
	 */
	public double containmentInterval(double gamma, double epsilon) throws MathException{
		if(gamma<=0 || gamma>=1 || epsilon <=0) throw new IllegalArgumentException();
		double x[] = new double[2];
		x[0] = 0;
		x[1] = incre;
		double val[] = new double[2];
		val[0] = 0;
		val[1] = cumulate(x[1]);
		// val[0] <= gamma <= val[1]を成立させる
		while(val[1]<gamma){
			for(int i=0;i<2;i++) x[i] += incre;
			val[0] = val[1];
			val[1] = cumulate(x[1]);
		}
		for(int i=0;i<2;i++)	if(Math.abs(val[i]-gamma)<epsilon) return x[i];
		while(true){
			double tmp = (x[0]+x[1])/2;
			double tmp2 = cumulate(tmp);
			if(Math.abs(tmp2-gamma)<epsilon) return tmp;
			int cnt = (tmp2<gamma) ? 0 :1 ;
			x[cnt] = tmp;
			val[cnt] = tmp2;
		}
	}
}
