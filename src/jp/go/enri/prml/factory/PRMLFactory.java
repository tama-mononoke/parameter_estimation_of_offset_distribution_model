/**
 * 
 */
package jp.go.enri.prml.factory;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Properties;
import java.util.StringTokenizer;

import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.RealVector;

import jp.go.enri.prml.bf.NDEParameter;
import jp.go.enri.prml.bf.ONDEParameter;
import jp.go.enri.prml.dist.NDE;
import jp.go.enri.prml.dist.ONDE;
import jp.go.enri.prml.vb.NDEParameterDistribution;
import jp.go.enri.prml.vb.ONDEParameterDistribution;



/**
 * 設定を読み取り、クラスを生成
 * @author 藤田雅人
 * @version 1.0.1　(Last update: 30/11/2011)
 */
public class PRMLFactory {
	private PRMLFactory(){}
	/**
	 * 設定より{@link ONDEParameterDistribution}クラスを構築
	 * @param prop 設定
	 * @return {@link ONDEParameterDistribution}クラス
	 */
	public static ONDEParameterDistribution getONDEInVariationalBayes(Properties prop){
		int L = Integer.parseInt(prop.getProperty("L"));
		int m = Integer.parseInt(prop.getProperty("m"));
		int n = Integer.parseInt(prop.getProperty("n"));
		if(m<0 || n<0 || L<=0) throw new IllegalArgumentException();
		int k = m + n;
		double alpha[] = new double[k];
		double a[] = new double[k];
		double b[] = new double[k];
		for(int i=0;i<k;i++){
			alpha[i] = Double.parseDouble(prop.getProperty("alpha" + (i+1)));
			a[i] = Double.parseDouble(prop.getProperty("a" + (i+1)));
			b[i] = Double.parseDouble(prop.getProperty("b" + (i+1)));
		}
		double offset[] = new double[L];
		double p[] = new double[L];
		for(int i=0;i<L;i++){
			offset[i] = Double.parseDouble(prop.getProperty("offset" + (i+1)));
			p[i] = Double.parseDouble(prop.getProperty("p" + (i+1)));
		}
		return new ONDEParameterDistribution(m,n,offset,p,alpha,a,b);
	}
	
	/**
	 * 設定より{@link ONDE}クラスを構築
	 * @param prop 設定
	 * @return {@link ONDE}クラス
	 */
	public static ONDE getONDE(Properties prop){
		int L = Integer.parseInt(prop.getProperty("L"));
		int m = Integer.parseInt(prop.getProperty("m"));
		int n = Integer.parseInt(prop.getProperty("n"));
		if(m<0 || n<0 || L<=0) throw new IllegalArgumentException();
		double alpha[] = new double[m+n];
		double sigma[] = new double[m];
		double lambda[] = new double[n];
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
		double offset[] = new double[L];
		double omega[] = new double[L];
		for(int i=0;i<L;i++){
			offset[i] = Double.parseDouble(prop.getProperty("offset" + (i+1)));
			omega[i] = Double.parseDouble(prop.getProperty("omega" + (i+1)));
		}
		return new ONDE(alpha,sigma,lambda,offset,omega);
	}
	
	/**
	 * 設定より{@link ONDEParameter}クラスを構築
	 * @param prop 設定
	 * @return {@link ONDEParameter}クラス
	 */
	public static ONDEParameter getONDEParameter(Properties prop){
		double offset[] = new double[Integer.parseInt(prop.getProperty("L"))];
		for(int i=0;i<offset.length;i++) offset[i] = Double.parseDouble(prop.getProperty("offset" + (i+1)));
		double initial_sigma_bound[][] = new double[Integer.parseInt(prop.getProperty("m"))][];
		double initial_lambda_bound[][] = new double[Integer.parseInt(prop.getProperty("n"))][];
		for(int i=0;i<initial_sigma_bound.length;i++){
			initial_sigma_bound[i] = new double[2];
			StringTokenizer st = new StringTokenizer(prop.getProperty("sigma" + (i+1) + "_bound"),",");
			for(int j=0;j<2;j++){
				initial_sigma_bound[i][j] = Double.parseDouble(st.nextToken());
			}
		}
		for(int i=0;i<initial_lambda_bound.length;i++){
			initial_lambda_bound[i] = new double[2];
			StringTokenizer st = new StringTokenizer(prop.getProperty("lambda" + (i+1) + "_bound"),",");
			for(int j=0;j<2;j++){
				initial_lambda_bound[i][j] = Double.parseDouble(st.nextToken());
			}
		}
		return new ONDEParameter(Integer.parseInt(prop.getProperty("D")), offset,initial_sigma_bound, initial_lambda_bound);
	}
	/**
	 * 設定より{@link NDEParameterDistribution}クラスを構築
	 * @param prop 設定
	 * @return {@link NDEParameterDistribution}クラス
	 */
	public static NDEParameterDistribution getNDEInVariationalBayes(Properties prop){
		int m = Integer.parseInt(prop.getProperty("m"));
		int n = Integer.parseInt(prop.getProperty("n"));
		if(m<0 || n<0) throw new IllegalArgumentException();
		int k = m + n;
		double alpha[] = new double[k];
		double a[] = new double[k];
		double b[] = new double[k];
		for(int i=0;i<k;i++){
			alpha[i] = Double.parseDouble(prop.getProperty("alpha" + (i+1)));
			a[i] = Double.parseDouble(prop.getProperty("a" + (i+1)));
			b[i] = Double.parseDouble(prop.getProperty("b" + (i+1)));
		}
		return new NDEParameterDistribution(m,n,alpha,a,b);
	}
	
	/**
	 * 設定より{@link NDE}クラスを構築
	 * @param prop 設定
	 * @return {@link NDE}クラス
	 */
	public static NDE getNDE(Properties prop){
		int m = Integer.parseInt(prop.getProperty("m"));
		int n = Integer.parseInt(prop.getProperty("n"));
		if(m<0 || n<0) throw new IllegalArgumentException();
		double alpha[] = new double[m+n];
		double sigma[] = new double[m];
		double lambda[] = new double[n];
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
		return new NDE(alpha,sigma,lambda);
	}
	
	/**
	 * 設定より{@link NDEParameter}クラスを構築
	 * @param prop 設定
	 * @return {@link NDEParameter}クラス
	 */
	public static NDEParameter getNDEParameter(Properties prop){
		double initial_sigma_bound[][] = new double[Integer.parseInt(prop.getProperty("m"))][];
		double initial_lambda_bound[][] = new double[Integer.parseInt(prop.getProperty("n"))][];
		for(int i=0;i<initial_sigma_bound.length;i++){
			initial_sigma_bound[i] = new double[2];
			StringTokenizer st = new StringTokenizer(prop.getProperty("sigma" + (i+1) + "_bound"),",");
			for(int j=0;j<2;j++){
				initial_sigma_bound[i][j] = Double.parseDouble(st.nextToken());
			}
		}
		for(int i=0;i<initial_lambda_bound.length;i++){
			initial_lambda_bound[i] = new double[2];
			StringTokenizer st = new StringTokenizer(prop.getProperty("lambda" + (i+1) + "_bound"),",");
			for(int j=0;j<2;j++){
				initial_lambda_bound[i][j] = Double.parseDouble(st.nextToken());
			}
		}
		return new NDEParameter(Integer.parseInt(prop.getProperty("D")), initial_sigma_bound, initial_lambda_bound);
	}
	/**
	 * ファイルから1次元数値データを読み込む
	 * @param fileName ファイル名
	 * @return データ
	 * @throws FileNotFoundException ファイルが見つからないとき
	 * @throws IOException 読み込みエラー
	 */
	public static double[] getData1(String fileName) throws FileNotFoundException, IOException{
		BufferedReader br = new BufferedReader(new FileReader(fileName));
		String line = br.readLine();
		RealVector list = new ArrayRealVector();
		while((line=br.readLine())!=null){
			StringTokenizer st = new StringTokenizer(line,",");
			list = list.append(Double.parseDouble(st.nextToken()));
		}
		br.close();
		return list.getData();
	}
}
