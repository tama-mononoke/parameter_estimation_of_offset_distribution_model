/**
 * 
 */
package jp.go.enri.prml.dist;

import cern.jet.random.AbstractDistribution;

/**
 * Gamma分布の乱数生成
 * @author masato
 *
 */
public class Gamma {
	/**
	 * パラメータ
	 */
	private double a;
	/**
	 * パラメータ
	 */
	private double b;
	/**
	 * 乱数生成器
	 */
	private AbstractDistribution dist = null;
	
	public Gamma(double a, double b){
		this.a = a;
		this.b = b;
	}
	/**
	 * 乱数の生成
	 * @return
	 */
	public double nextDouble(){
		if(dist==null){
			dist = new cern.jet.random.Gamma(a,b,MyRandomEngine.getRandomEngine());
		}
		return dist.nextDouble();
	}
}
