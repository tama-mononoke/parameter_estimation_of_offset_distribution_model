/**
 * 
 */
package jp.go.enri.prml.dist;

import cern.jet.random.AbstractDistribution;

/**
 * Random sample generator of Gamma distribution
 * @author Masato Fujita
 *
 */
public class Gamma {
	/**
	 * Parameter
	 */
	private double a;
	/**
	 * Parameter
	 */
	private double b;
	/**
	 * Random sample generator
	 */
	private AbstractDistribution dist = null;
	
	public Gamma(double a, double b){
		this.a = a;
		this.b = b;
	}
	/**
	 * Generate random samples.
	 * @return samples
	 */
	public double nextDouble(){
		if(dist==null){
			dist = new cern.jet.random.Gamma(a,b,MyRandomEngine.getRandomEngine());
		}
		return dist.nextDouble();
	}
}
