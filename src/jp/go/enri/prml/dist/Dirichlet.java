/**
 * 
 */
package jp.go.enri.prml.dist;

import java.util.Arrays;

import cern.jet.random.AbstractDistribution;

/**
 * Dirichlet分布の乱数生成
 * @author masato
 *
 */
public class Dirichlet {
	/**
	 * パラメータ
	 */
	private double alpha[];

	public Dirichlet(double[] alpha) {
		super();
		this.alpha = Arrays.copyOf(alpha, alpha.length);
	}
	
	/**
	 * 乱数生成器
	 */
	private AbstractDistribution dist[] = null;
	/**
	 * 乱数の生成
	 * @return
	 */
	public double[] nextDouble(){
		if(dist==null){
			dist = new cern.jet.random.Gamma[alpha.length];
			for(int i=0;i<dist.length;i++){
				dist[i] = new cern.jet.random.Gamma(alpha[i],1,MyRandomEngine.getRandomEngine());
			}
		}
		double ans[] = new double[dist.length];
		double sum = 0;
		for(int i=0;i<dist.length;i++){
			ans[i]  = dist[i].nextDouble();
			sum += ans[i];
		}
		for(int i=0;i<dist.length;i++) ans[i] /= sum;
		return ans;
	}
}
