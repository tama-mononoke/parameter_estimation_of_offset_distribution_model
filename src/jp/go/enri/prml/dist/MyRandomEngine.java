/**
 * 
 */
package jp.go.enri.prml.dist;

import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

/**
 * 乱数生成器
 * @author 藤田雅人
 *
 */
public class MyRandomEngine {
	/**
	 * 乱数生成器
	 */
	private static RandomEngine engine = null;
	private MyRandomEngine(){}
	/**
	 * 乱数生成器を取得
	 * @return 乱数生成器
	 */
	public static RandomEngine getRandomEngine(){
		if(engine==null) engine = new MersenneTwister();
		return engine;
	}
}
