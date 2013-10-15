/**
 * 
 */
package jp.go.enri.prml;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.Date;
import java.util.Properties;
import java.util.StringTokenizer;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import jp.go.enri.prml.dist.NDE;
import jp.go.enri.prml.dist.ONDE;
import jp.go.enri.prml.em.NDEEM;
import jp.go.enri.prml.em.ONDEEM;
import jp.go.enri.prml.bf.NDEBF;
import jp.go.enri.prml.bf.ONDEBF;
import jp.go.enri.prml.bf.NDEParameter;
import jp.go.enri.prml.bf.ONDEParameter;
import jp.go.enri.prml.factory.PRMLFactory;
import jp.go.enri.prml.factory.TestDataGeneratorNDE;
import jp.go.enri.prml.factory.TestDataGeneratorONDE;
import jp.go.enri.prml.vb.NDEParameterDistribution;
import jp.go.enri.prml.vb.ONDEParameterDistribution;
import jp.go.enri.prml.vb.NDEVB;
import jp.go.enri.prml.vb.ONDEVB;

/**
 * Parameter estimation entry
 * @author Masato Fujita
 *
 */
public class BayesianUpdate {
	/**
	 * Log
	 */
	public static Log log = LogFactory.getLog(BayesianUpdate.class);
	/**
	 * Constructor
	 */
	private BayesianUpdate(){}
	
	/**
	 * Run the parameter estimation engine and related functions 
	 * @param args arg[0]: nde: non-offset onde: offset
	 * @throws Exception
	 */
	public static void run(String args[]) throws Exception{
		String tmp[] = Arrays.copyOfRange(args, 1, args.length);
		if(args[0].equals("nde")) NDE(tmp);
		else if(args[0].equals("onde")) ONDE(tmp);
		else throw new IllegalArgumentException();
	}
	
	/**
	 * ONDE型モデル
	 * @param args args[0]: rg: 乱数生成, vb:変分ベイズ, em:EM法, bf:総当たり法, oe:オフセット値推定
	 * @throws Exception
	 */
	private static void ONDE(String args[]) throws Exception{
		String tmp[] = Arrays.copyOfRange(args, 1, args.length);
		if(args[0].equals("rg")) ONDERandom(tmp);
		else if(args[0].equals("vb")) ONDEVB(tmp);
		else if(args[0].equals("em")) ONDEEM(tmp);
		else if(args[0].equals("bf")) ONDEBF(tmp);
		else if(args[0].equals("oe")) ONDEOffset(tmp);
		else throw new IllegalArgumentException();
	}
	/**
	 * オフセット推定
	 * @param args args[0]: d: 分布モデルを用いた推定, b: ベイズモデルパラメータ分布を用いた推定
	 * @throws Exception
	 */
	private static void ONDEOffset(String args[]) throws Exception{
		String tmp[] = Arrays.copyOfRange(args, 1, args.length);
		if(args[0].equals("d")) ONDEOffset1(tmp);
		else if(args[0].equals("b")) ONDEOffset2(tmp);
		else throw new IllegalArgumentException();
		
	}
	
	/**
	 * 分布モデルを用いたオフセット推定
	 * @param args args[0]:分布を定義したプロパティファイル, args[1]：データファイル, args[2]:出力ファイル
	 * @throws Exception
	 */
	private static void ONDEOffset1(String args[]) throws Exception{
		FileInputStream inStream = new FileInputStream(args[0]);
		Properties propDist = new Properties();
		propDist.load(inStream);
		inStream.close();
		ONDE dist = PRMLFactory.getONDE(propDist);
		BufferedReader br = new BufferedReader(new FileReader(args[1]));
		BufferedWriter bw = new BufferedWriter(new FileWriter(args[2]));
		String line = "";
		while((line=br.readLine())!=null){
			StringTokenizer st = new StringTokenizer(line,",");
			double est[] = dist.estimateOffsetProbability(Double.parseDouble(st.nextToken()));
			StringBuilder sb = new StringBuilder(line);
			for(double tmp : est){
				sb.append(",");
				sb.append(tmp);
			}
			sb.append("\n");
			bw.write(sb.toString());
		}
		bw.close();
		br.close();
	}
	
	/**
	 * ベイズモデルパラメータ分布を用いたオフセット推定
	 * @param args args[0]:分布を定義したプロパティファイル, args[1]：データファイル, args[2]:出力ファイル
	 * @throws Exception
	 */
	private static void ONDEOffset2(String args[]) throws Exception{
		FileInputStream inStream = new FileInputStream(args[0]);
		Properties propDist = new Properties();
		propDist.load(inStream);
		inStream.close();
		ONDEParameterDistribution dist = PRMLFactory.getONDEInVariationalBayes(propDist);
		BufferedReader br = new BufferedReader(new FileReader(args[1]));
		BufferedWriter bw = new BufferedWriter(new FileWriter(args[2]));
		String line = "";
		while((line=br.readLine())!=null){
			StringTokenizer st = new StringTokenizer(line,",");
			double est[] = dist.estimateOffsetProbability(Double.parseDouble(st.nextToken()));
			StringBuilder sb = new StringBuilder(line);
			for(double tmp : est){
				sb.append(",");
				sb.append(tmp);
			}
			sb.append("\n");
			bw.write(sb.toString());
		}
		bw.close();
		br.close();
	}
	
	/**
	 * ONDE型分布モデルサンプル生成
	 * @param args args[0]:分布を定義したプロパティファイル, args[1]：データ出力個数, args[2]:出力ファイル
	 */
	private static void ONDERandom(String args[]) throws Exception{
		for(int i=0;i<args.length;i++){
			log.info(args[i]);
		}
		int num = Integer.parseInt(args[1]);
		if(num<=0) throw new IllegalArgumentException();
		FileInputStream inStream = new FileInputStream(args[0]);
		Properties propDist = new Properties();
		propDist.load(inStream);
		inStream.close();
		TestDataGeneratorONDE generator = new TestDataGeneratorONDE(propDist);
		BufferedWriter bw = new BufferedWriter(new FileWriter(args[2]));
		for(int i=0;i<num;i++){
			double out[] = generator.generate();
			bw.write(out[0] + "," + out[1]  + "," + out[2] + "\n");
		}
		bw.close();
	}
	
	/**
	 * 総当たり法によるONDE型分布モデルのパラメータ最尤推定
	 * @param args args[0]:プロパティファイル, args[1]：データファイル, args[2]:出力ファイル
	 */
	private static void ONDEBF(String args[]) throws Exception{
		for(int i=0;i<args.length;i++){
			log.info(args[i]);
		}
		double data[] = PRMLFactory.getData1(args[1]);
		FileInputStream inStream2 = new FileInputStream(args[0]);
		Properties propDist2 = new Properties();
		propDist2.load(inStream2);
		ONDEParameter initial = PRMLFactory.getONDEParameter(propDist2);
		ONDEBF estimator = new ONDEBF();
		log.info("bf-estimation process start");
		ONDEBF.Result res = estimator.estimate(initial,data);
		log.info("bf-estimation process finish");
		ONDE para = res.getParam();
		
		Properties out = new Properties();
		out.setProperty("L", Integer.toString(para.getL()));
		out.setProperty("m", Integer.toString(para.getM()));
		out.setProperty("n", Integer.toString(para.getN()));
		log.info("BF Result");
		for(int i=0;i<para.getL();i++){
			log.info("offset[" + i + "]=" + para.getOffset()[i]);
			log.info("omega[" + i + "]=" + para.getOmega()[i]);
			out.setProperty("offset" + (i+1), Double.toString(para.getOffset()[i]));
			out.setProperty("omega" + (i+1), Double.toString(para.getOmega()[i]));
		}
		for(int i=0;i<para.getM();i++){
			log.info("pi[" + i + "]=" + para.getPi()[i]);
			log.info("sigma[" + i + "]=" + para.getSigma()[i]);
			out.setProperty("pi_N" + (i+1), Double.toString(para.getPi()[i]));
			out.setProperty("sigma" + (i+1), Double.toString(para.getSigma()[i]));
		}
		for(int i=0;i<para.getN();i++){
			log.info("pi[" + (i+para.getM()) + "]=" + para.getPi()[i+para.getM()]);
			log.info("lambda[" + i + "]=" + para.getLambda()[i]);
			out.setProperty("pi_DE" + (i+1), Double.toString(para.getPi()[i+para.getM()]));
			out.setProperty("lambda" + (i+1), Double.toString(para.getLambda()[i]));
		}
		log.info("logLiklihood=" + res.getLogLiklihood());
		out.setProperty("logLiklihood", Double.toString(res.getLogLiklihood()));
		BufferedWriter stream = new BufferedWriter(new FileWriter(args[2]));
		out.store(stream, "BaysienUpdate generated automatically in " + new Date(System.currentTimeMillis()));
		stream.close();
	}
	
	/**
	 * EM法によるONDE型分布モデルのパラメータ最尤推定
	 * @param args args[0]:プロパティファイル, args[1]：データファイル, args[2]:出力ファイル
	 */
	private static void ONDEEM(String args[]) throws Exception{
		for(int i=0;i<args.length;i++){
			log.info(args[i]);
		}
		double data[] = PRMLFactory.getData1(args[1]);

		FileInputStream inStream2 = new FileInputStream(args[0]);
		Properties propDist2 = new Properties();
		propDist2.load(inStream2);
		ONDE initial = PRMLFactory.getONDE(propDist2);
		ONDEEM estimator = new ONDEEM();
		log.info("em-estimation process start");
		ONDEEM.Result res = estimator.estimate(initial,data);
		log.info("em-estimation process finish");
		ONDE para = res.getParam();
		
		Properties out = new Properties();
		out.setProperty("L", Integer.toString(para.getL()));
		out.setProperty("m", Integer.toString(para.getM()));
		out.setProperty("n", Integer.toString(para.getN()));
		log.info("EM Result");
		for(int i=0;i<para.getL();i++){
			log.info("offset[" + i + "]=" + para.getOffset()[i]);
			log.info("omega[" + i + "]=" + para.getOmega()[i]);
			out.setProperty("offset" + (i+1), Double.toString(para.getOffset()[i]));
			out.setProperty("omega" + (i+1), Double.toString(para.getOmega()[i]));
		}
		for(int i=0;i<para.getM();i++){
			log.info("pi[" + i + "]=" + para.getPi()[i]);
			log.info("sigma[" + i + "]=" + para.getSigma()[i]);
			out.setProperty("pi_N" + (i+1), Double.toString(para.getPi()[i]));
			out.setProperty("sigma" + (i+1), Double.toString(para.getSigma()[i]));
		}
		for(int i=0;i<para.getN();i++){
			log.info("pi[" + (i+para.getM()) + "]=" + para.getPi()[i+para.getM()]);
			log.info("lambda[" + i + "]=" + para.getLambda()[i]);
			out.setProperty("pi_DE" + (i+1), Double.toString(para.getPi()[i+para.getM()]));
			out.setProperty("lambda" + (i+1), Double.toString(para.getLambda()[i]));
		}
		log.info("logLiklihood=" + res.getLogLiklihood());
		out.setProperty("logLiklihood", Double.toString(res.getLogLiklihood()));
		BufferedWriter stream = new BufferedWriter(new FileWriter(args[2]));
		out.store(stream, "BaysienUpdate generated automatically in " + new Date(System.currentTimeMillis()));
		stream.close();
	}
	
	/**
	 * 変分ベイズ法によるONDE型分布モデルのパラメータ最尤推定
	 * @param args args[0]:プロパティファイル, args[1]：データファイル, args[2]:出力ファイル
	 */
	private static void ONDEVB(String args[]) throws Exception{
		for(int i=0;i<args.length;i++){
			log.info(args[i]);
		}
		double data[] = PRMLFactory.getData1(args[1]);
		FileInputStream inStream2 = new FileInputStream(args[0]);
		Properties propDist2 = new Properties();
		propDist2.load(inStream2);
		ONDEParameterDistribution initial =  PRMLFactory.getONDEInVariationalBayes(propDist2);
		ONDEVB estimator = new ONDEVB();
		log.info("vb-estimation process start");
		ONDEVB.Result res = estimator.estimate(initial,data);
		log.info("vb-estimation process finish");
		ONDE para = res.getParam().MAPEstimation();
		
		Properties out = new Properties();
		log.info("Result");
		log.info("VB Posterior");
		log.info("m=" + res.getParam().getM());
		log.info("n=" + res.getParam().getN());
		out.setProperty("L", Integer.toString(para.getL()));
		out.setProperty("m", Integer.toString(res.getParam().getM()));
		out.setProperty("n", Integer.toString(res.getParam().getN()));
		double p[] = res.getParam().getP();
		for(int i=0;i<p.length;i++){
			log.info("p[" + i + "]=" + p[i]);
			out.setProperty("p" + (i+1), Double.toString(p[i]));
		}
		double alpha[] = res.getParam().getAlpha();
		double A[] = res.getParam().getA();
		double B[] = res.getParam().getB();
		for(int i=0;i<res.getParam().getM()+res.getParam().getN();i++){
			log.info("alpha[" + i + "]=" + alpha[i]);
			log.info("a[" + i + "]=" + A[i]);
			log.info("b[" + i + "]=" + B[i]);
			out.setProperty("alpha" + (i+1), Double.toString(alpha[i]));
			out.setProperty("a" + (i+1), Double.toString(A[i]));
			out.setProperty("b" + (i+1), Double.toString(B[i]));
		}
		log.info("VB Result");
		for(int i=0;i<para.getL();i++){
			log.info("offset[" + i + "]=" + para.getOffset()[i]);
			log.info("omega[" + i + "]=" + para.getOmega()[i]);
			out.setProperty("offset" + (i+1), Double.toString(para.getOffset()[i]));
			out.setProperty("omega" + (i+1), Double.toString(para.getOmega()[i]));
		}
		for(int i=0;i<para.getM();i++){
			log.info("pi[" + i + "]=" + para.getPi()[i]);
			log.info("sigma[" + i + "]=" + para.getSigma()[i]);
			out.setProperty("pi_N" + (i+1), Double.toString(para.getPi()[i]));
			out.setProperty("sigma" + (i+1), Double.toString(para.getSigma()[i]));
		}
		for(int i=0;i<para.getN();i++){
			log.info("pi[" + (i+para.getM()) + "]=" + para.getPi()[i+para.getM()]);
			log.info("lambda[" + i + "]=" + para.getLambda()[i]);
			out.setProperty("pi_DE" + (i+1), Double.toString(para.getPi()[i+para.getM()]));
			out.setProperty("lambda" + (i+1), Double.toString(para.getLambda()[i]));
		}
		log.info("lowerbound=" + res.getLowerbound());
		out.setProperty("lowerbound", Double.toString(res.getLowerbound()));
		BufferedWriter stream = new BufferedWriter(new FileWriter(args[2]));
		out.store(stream, "BaysienUpdate generated automatically in " + new Date(System.currentTimeMillis()));
		stream.close();
	}
	
	/**
	 * NDE型モデル
	 * @param args args[0]: rg: 乱数生成, vb:変分ベイズ, em:EM法, bf:総当たり法
	 * @throws Exception
	 */
	private static void NDE(String args[]) throws Exception{
		String tmp[] = Arrays.copyOfRange(args, 1, args.length);
		if(args[0].equals("rg")) NDERandom(tmp);
		else if(args[0].equals("vb")) NDEVB(tmp);
		else if(args[0].equals("em")) NDEEM(tmp);
		else if(args[0].equals("bf")) NDEBF(tmp);
		else throw new IllegalArgumentException();
	}
	
	/**
	 * NDE型分布モデルサンプル生成
	 * @param args args[0]:分布を定義したプロパティファイル, args[1]：データ出力個数, args[2]:出力ファイル
	 */
	private static void NDERandom(String args[]) throws Exception{
		for(int i=0;i<args.length;i++){
			log.info(args[i]);
		}
		int num = Integer.parseInt(args[1]);
		if(num<=0) throw new IllegalArgumentException();
		FileInputStream inStream = new FileInputStream(args[0]);
		Properties propDist = new Properties();
		propDist.load(inStream);
		inStream.close();
		TestDataGeneratorNDE generator = new TestDataGeneratorNDE(propDist);
		BufferedWriter bw = new BufferedWriter(new FileWriter(args[2]));
		for(int i=0;i<num;i++){
			double out[] = generator.generate();
			bw.write(out[0] + "," + out[1] + "\n");
		}
		bw.close();
	}
	
	/**
	 * 総当たり法によるNDE型分布モデルのパラメータ最尤推定
	 * @param args args[0]:プロパティファイル, args[1]：データファイル, args[2]:出力ファイル
	 */
	private static void NDEBF(String args[]) throws Exception{
		for(int i=0;i<args.length;i++){
			log.info(args[i]);
		}
		double data[] = PRMLFactory.getData1(args[1]);
		FileInputStream inStream2 = new FileInputStream(args[0]);
		Properties propDist2 = new Properties();
		propDist2.load(inStream2);
		NDEParameter initial = PRMLFactory.getNDEParameter(propDist2);
		NDEBF estimator = new NDEBF();
		log.info("bf-estimation process start");
		NDEBF.Result res = estimator.estimate(initial,data);
		log.info("bf-estimation process finish");
		NDE para = res.getParam();
		
		Properties out = new Properties();
		out.setProperty("m", Integer.toString(para.getM()));
		out.setProperty("n", Integer.toString(para.getN()));
		log.info("BF Result");
		for(int i=0;i<para.getM();i++){
			log.info("pi[" + i + "]=" + para.getPi()[i]);
			log.info("sigma[" + i + "]=" + para.getSigma()[i]);
			out.setProperty("pi_N" + (i+1), Double.toString(para.getPi()[i]));
			out.setProperty("sigma" + (i+1), Double.toString(para.getSigma()[i]));
		}
		for(int i=0;i<para.getN();i++){
			log.info("pi[" + (i+para.getM()) + "]=" + para.getPi()[i+para.getM()]);
			log.info("lambda[" + i + "]=" + para.getLambda()[i]);
			out.setProperty("pi_DE" + (i+1), Double.toString(para.getPi()[i+para.getM()]));
			out.setProperty("lambda" + (i+1), Double.toString(para.getLambda()[i]));
		}
		log.info("logLiklihood=" + res.getLogLiklihood());
		out.setProperty("logLiklihood", Double.toString(res.getLogLiklihood()));
		BufferedWriter stream = new BufferedWriter(new FileWriter(args[2]));
		out.store(stream, "BaysienUpdate generated automatically in " + new Date(System.currentTimeMillis()));
		stream.close();
	}
	
	/**
	 * EM法によるNDE型分布モデルのパラメータ最尤推定
	 * @param args args[0]:プロパティファイル, args[1]：データファイル, args[2]:出力ファイル
	 */
	private static void NDEEM(String args[]) throws Exception{
		for(int i=0;i<args.length;i++){
			log.info(args[i]);
		}
		double data[] = PRMLFactory.getData1(args[1]);

		FileInputStream inStream2 = new FileInputStream(args[0]);
		Properties propDist2 = new Properties();
		propDist2.load(inStream2);
		NDE initial = PRMLFactory.getNDE(propDist2);
		NDEEM estimator = new NDEEM();
		log.info("em-estimation process start");
		NDEEM.Result res = estimator.estimate(initial,data);
		log.info("em-estimation process finish");
		NDE para = res.getParam();
		
		Properties out = new Properties();
		out.setProperty("m", Integer.toString(para.getM()));
		out.setProperty("n", Integer.toString(para.getN()));
		log.info("EM Result");
		for(int i=0;i<para.getM();i++){
			log.info("pi[" + i + "]=" + para.getPi()[i]);
			log.info("sigma[" + i + "]=" + para.getSigma()[i]);
			out.setProperty("pi_N" + (i+1), Double.toString(para.getPi()[i]));
			out.setProperty("sigma" + (i+1), Double.toString(para.getSigma()[i]));
		}
		for(int i=0;i<para.getN();i++){
			log.info("pi[" + (i+para.getM()) + "]=" + para.getPi()[i+para.getM()]);
			log.info("lambda[" + i + "]=" + para.getLambda()[i]);
			out.setProperty("pi_DE" + (i+1), Double.toString(para.getPi()[i+para.getM()]));
			out.setProperty("lambda" + (i+1), Double.toString(para.getLambda()[i]));
		}
		log.info("logLiklihood=" + res.getLogLiklihood());
		out.setProperty("logLiklihood", Double.toString(res.getLogLiklihood()));
		BufferedWriter stream = new BufferedWriter(new FileWriter(args[2]));
		out.store(stream, "BaysienUpdate generated automatically in " + new Date(System.currentTimeMillis()));
		stream.close();
	}
	
	/**
	 * 変分ベイズ法によるNDE型分布モデルのパラメータ最尤推定
	 * @param args args[0]:プロパティファイル, args[1]：データファイル, args[2]:出力ファイル
	 */
	private static void NDEVB(String args[]) throws Exception{
		for(int i=0;i<args.length;i++){
			log.info(args[i]);
		}
		double data[] = PRMLFactory.getData1(args[1]);
		FileInputStream inStream2 = new FileInputStream(args[0]);
		Properties propDist2 = new Properties();
		propDist2.load(inStream2);
		NDEParameterDistribution initial =  PRMLFactory.getNDEInVariationalBayes(propDist2);
		NDEVB estimator = new NDEVB();
		log.info("vb-estimation process start");
		NDEVB.Result res = estimator.estimate(initial,data);
		log.info("vb-estimation process finish");
		NDE para = res.getParam().MAPEstimation();
		
		Properties out = new Properties();
		log.info("Result");
		log.info("VB Posterior");
		log.info("m=" + res.getParam().getM());
		log.info("n=" + res.getParam().getN());
		out.setProperty("m", Integer.toString(res.getParam().getM()));
		out.setProperty("n", Integer.toString(res.getParam().getN()));
		double alpha[] = res.getParam().getAlpha();
		double A[] = res.getParam().getA();
		double B[] = res.getParam().getB();
		for(int i=0;i<res.getParam().getM()+res.getParam().getN();i++){
			log.info("alpha[" + i + "]=" + alpha[i]);
			log.info("a[" + i + "]=" + A[i]);
			log.info("b[" + i + "]=" + B[i]);
			out.setProperty("alpha" + (i+1), Double.toString(alpha[i]));
			out.setProperty("a" + (i+1), Double.toString(A[i]));
			out.setProperty("b" + (i+1), Double.toString(B[i]));
		}
		log.info("VB Result");
		for(int i=0;i<para.getM();i++){
			log.info("pi[" + i + "]=" + para.getPi()[i]);
			log.info("sigma[" + i + "]=" + para.getSigma()[i]);
			out.setProperty("pi_N" + (i+1), Double.toString(para.getPi()[i]));
			out.setProperty("sigma" + (i+1), Double.toString(para.getSigma()[i]));
		}
		for(int i=0;i<para.getN();i++){
			log.info("pi[" + (i+para.getM()) + "]=" + para.getPi()[i+para.getM()]);
			log.info("lambda[" + i + "]=" + para.getLambda()[i]);
			out.setProperty("pi_DE" + (i+1), Double.toString(para.getPi()[i+para.getM()]));
			out.setProperty("lambda" + (i+1), Double.toString(para.getLambda()[i]));
		}
		log.info("lowerbound=" + res.getLowerbound());
		out.setProperty("lowerbound", Double.toString(res.getLowerbound()));
		BufferedWriter stream = new BufferedWriter(new FileWriter(args[2]));
		out.store(stream, "BaysienUpdate generated automatically in " + new Date(System.currentTimeMillis()));
		stream.close();
	}
	
	
}
