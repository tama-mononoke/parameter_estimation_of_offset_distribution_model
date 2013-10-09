/**
 * 
 */
package jp.go.enri.prml.factory;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.Properties;

import jp.go.enri.prml.dist.NDE;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.math.random.MersenneTwister;
import org.apache.commons.math.random.RandomGenerator;

/**
 * {@link VariationBayesMLM}用テストデータ生成器
 * @author 藤田雅人
 * 
 */
public class TestDataGeneratorNDE {
	/**
	 * ログの取得
	 */
	public static Log log = LogFactory.getLog(TestDataGeneratorNDE.class);
	/**
	 * 正規分布の個数
	 */
	private final int m;
	/**
	 * 両側指数分布の個数
	 */
	private final int n;
	/**
	 * alphaの値
	 */
	private double alpha[];
	/**
	 * lambdaの値
	 */
	private double lambda[];
	/**
	 * sigmaの値
	 */
	private double sigma[];
	/**
	 * 計算用変数
	 */
	private double tmp_alpha[];
	/**
	 * 乱数生成器
	 */
	private RandomGenerator random;
	/**
	 * コンストラクタ
	 * @param prop プロパティ
	 */
	public TestDataGeneratorNDE(Properties prop){
		m = Integer.parseInt(prop.getProperty("m"));
		n = Integer.parseInt(prop.getProperty("n"));
		if(m<0 || n<0) throw new IllegalArgumentException();
		alpha = new double[m+n];
		sigma = new double[m];
		lambda = new double[n];
		double sum = 0;
		for(int i=0;i<m;i++){
			alpha[i] = Double.parseDouble(prop.getProperty("pi_N" + Integer.toString(i+1)));
			sigma[i] = Double.parseDouble(prop.getProperty("sigma" + Integer.toString(i+1)));
			if(alpha[i]<0 || sigma[i]<=0) throw new IllegalArgumentException();
			sum += alpha[i];
		}
		for(int i=0;i<n;i++){
			alpha[m+i] = Double.parseDouble(prop.getProperty("pi_DE" + Integer.toString(i+1)));
			lambda[i] = Double.parseDouble(prop.getProperty("lambda" + Integer.toString(i+1)));
			if(alpha[m+i]<0 || lambda[i]<=0) throw new IllegalArgumentException();
			sum += alpha[m+i];
		}
		double tmp = 0;
		tmp_alpha = new double[m+n];
		for(int i=0;i<m+n;i++){
			tmp += alpha[i];
			tmp_alpha[i] = tmp/sum;
		}
		random = new MersenneTwister();
	}
	
	/**
	 * コンストラクタ
	 * @param nde NDEインスタンス
	 */
	public TestDataGeneratorNDE(NDE nde){
		m = nde.getM();
		n = nde.getN();
		alpha = Arrays.copyOf(nde.getPi(), nde.getPi().length);
		sigma = Arrays.copyOf(nde.getSigma(), nde.getSigma().length);
		lambda = Arrays.copyOf(nde.getLambda(), nde.getLambda().length);
		double sum = 0;
		for(int i=0;i<m+n;i++) sum += alpha[i];
		double tmp = 0;
		tmp_alpha = new double[m+n];
		for(int i=0;i<m+n;i++){
			tmp += alpha[i];
			tmp_alpha[i] = tmp/sum;
		}
		random = new MersenneTwister();
	}
	
	
	/**
	 * 乱数の生成
	 * @return 擬似乱数とそれを生成した分布の識別番号
	 */
	public double[] generate(){
		double tmp = random.nextDouble();
		for(int i=0;i<tmp_alpha.length;i++){
			if(tmp<=tmp_alpha[i]){
				double out[] = new double[2];
				out[1] = i;
				if(i<m){
					out[0] = sigma[i]*random.nextGaussian();
					return out;
				}
				else{
					double ans = 0.0;
					while(ans==0.0){
						ans = random.nextDouble();
					}
					ans = lambda[i-m]*Math.log(ans);
					if(random.nextBoolean()) ans *= -1;
					out[0] = ans;
					return out;
				}
			}
		}
		throw new IllegalArgumentException();
	}
	
	/**
	 * 本クラスの確認用プログラム。生成した乱数とそれを生成した分布モデルの番号を出力
	 * @param args args[0]:分布を定義したプロパティファイル, args[1]：出力ファイル, args[2]:データ出力個数
	 */
	public static void main(String[] args) {
		try{
			int num = Integer.parseInt(args[2]);
			if(num<=0) throw new IllegalArgumentException();
			FileInputStream inStream = new FileInputStream(args[0]);
			Properties propDist = new Properties();
			propDist.load(inStream);
			inStream.close();
			TestDataGeneratorNDE generator = new TestDataGeneratorNDE(propDist);
			BufferedWriter bw = new BufferedWriter(new FileWriter(args[1]));
			for(int i=0;i<num;i++){
				double out[] = generator.generate();
				bw.write(out[0] + "," + out[1] + "\n");
			}
			bw.close();
		}
		catch(Exception ex){
			ex.printStackTrace();
		}

	}
}
