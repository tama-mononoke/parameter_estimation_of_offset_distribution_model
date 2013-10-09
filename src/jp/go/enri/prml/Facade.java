/**
 * 
 */
package jp.go.enri.prml;

import java.util.Arrays;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.log4j.PropertyConfigurator;

/**
 * プログラムの実行
 * @author 藤田雅人
 *
 */
public class Facade {
	/**
	 * ログの取得
	 */
	public static Log log = LogFactory.getLog(Facade.class);
	/**
	 * 実行
	 * @param args
	 * arg[0]: "u" ベイズ更新
	 */
	public static void main(String[] args) {
		PropertyConfigurator.configure( "log4j.properties" );
		try{
			String tmp[] = Arrays.copyOfRange(args, 1, args.length);
			if(args[0].equals("u")) BayesianUpdate.run(tmp);
			else throw new IllegalArgumentException();
		}
		catch(Exception ex){
			log.error("Error!", ex);
			ex.printStackTrace();
		}
	}

}
