/**
 * 
 */
package jp.go.enri.prml.dist;

import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

/**
 * Random sample generator
 * @author Masato Fujita
 *
 */
public class MyRandomEngine {
	/**
	 * Random sample generator
	 */
	private static RandomEngine engine = null;
	private MyRandomEngine(){}
	/**
	 * Get the random sample generator.
	 * @return Random sample generator
	 */
	public static RandomEngine getRandomEngine(){
		if(engine==null) engine = new MersenneTwister();
		return engine;
	}
}
