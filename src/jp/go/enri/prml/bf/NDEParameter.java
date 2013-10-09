/**
 * 
 */
package jp.go.enri.prml.bf;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * 総当たり法によるガウス分布とラプラス分布の混合分布のパラメータ最尤推定に用いるパラメータ
 * @author 藤田雅人
 *
 */
public class NDEParameter {
	/**
	 * ログの取得
	 */
	public static Log log = LogFactory.getLog(NDEParameter.class);
	/**
	 * 区間分割数
	 */
	int D;
	/**
	 * シグマのとりうる範囲
	 */
	double initial_sigma_bound[][];
	/**
	 * ラムダのとりうる範囲
	 */
	double initial_lambda_bound[][];
	/**
	 * コンストラクタ
	 * @param d
	 * @param initial_sigma_bound
	 * @param initial_lambda_bound
	 */
	public NDEParameter(int d, double[][] initial_sigma_bound,
			double[][] initial_lambda_bound) {
		super();
		D = d;
		if(D<=0 || initial_lambda_bound==null || initial_sigma_bound==null){
			log.error("invalid argument");
			throw new IllegalArgumentException();
		}
		this.initial_sigma_bound = initial_sigma_bound;
		this.initial_lambda_bound = initial_lambda_bound;
	}
	
	
	
}
