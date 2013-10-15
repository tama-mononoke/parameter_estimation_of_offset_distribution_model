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
 * {@link VariationBayesMLM} test data generator
 * @author Masato Fujita
 * 
 */
public class TestDataGeneratorNDE {
	/**
	 * Log
	 */
	public static Log log = LogFactory.getLog(TestDataGeneratorNDE.class);
	/**
	 * the number of Gaussian components in the mixture distribution.
	 */
	private final int m;
	/**
	 * the number of Laplace components in the mixture distribution.
	 */
	private final int n;
	/**
	 * the mixing coefficients
	 */
	private double alpha[];
	/**
	 * the scale parameter of Laplace distributions
	 */
	private double lambda[];
	/**
	 * the standard deviation of Gaussian distributions
	 */
	private double sigma[];
	/**
	 * temporary variable
	 */
	private double tmp_alpha[];
	/**
	 * Random sample generator
	 */
	private RandomGenerator random;
	/**
	 * Constructor
	 * @param prop configuration
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
	 * Constructor
	 * @param nde NDE instance
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
	 * Generate random samples.
	 * @return random samples and the identification number of the generating distribution component.
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
	 * Test program. Random samples and the identification number of the generating distribution component are recorded.
	 * @param args args[0]:Configuration file defining distribution, args[1]：output file, args[2]: NUmber of generated samples.
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
