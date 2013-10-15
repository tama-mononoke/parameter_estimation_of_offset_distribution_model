/**
 *
 */
package jp.go.enri.util.crm.integral;

import jp.go.enri.prml.dist.NDE;

import org.apache.commons.math.MathException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.special.Erf;

/**
 * Calculate the lateral overlap probability when the distribution type of TSE is specified.
 * All functions are implemented based on SASP-WG/WHL/11-WP/05 and SASP-WG/WHL/9-WP/14.
 * @author Masato Fujita
 * @version 1.1
 * 2009/04/21<br/>
 * 2011/11/22<br/>
 * 2013/10/15<br/>
 */
public class LateralOverlapProbability {
	
	/**
	 * Two numbers whose gap is smaller than this value are identified.
	 */
	private static final double epsilon = Math.pow(10,-5);
	/**
	 * Lateral overlap probability
	 * @param nde1 Distribution model
	 * @param nde2 Distribution model
	 * @param S Route spacing(NM)
	 * @param w Wing Span
	 * @return Lateral overlap probability
	 */
	public static double evaluateLateralOverlapProbability(NDE nde1, NDE nde2, double S, double w) throws MathException{
		int m1 = nde1.getM();
		int n1 = nde1.getN();
		int m2 = nde2.getM();
		int n2 = nde2.getN();
		double pi1[] = nde1.getPi();
		double sigma1[] = nde1.getSigma();
		double lambda1[] = nde1.getLambda();
		double pi2[] = nde1.getPi();
		double sigma2[] = nde1.getSigma();
		double lambda2[] = nde1.getLambda();
		double sum = 0;
		for(int i=0;i<m1;i++){
			for(int j=0;j<m2;j++){
				sum += pi1[i]*pi2[j]*PySyN(S,w,sigma1[i],sigma2[j]);
			}
		}
		for(int i=0;i<m1;i++){
			for(int j=0;j<n2;j++){
				sum += pi1[i]*pi2[m1+j]*PySyNVSDE(S,w,sigma1[i],lambda2[j]);
			}
		}
		for(int i=0;i<n1;i++){
			for(int j=0;j<m2;j++){
				sum += pi1[m1+i]*pi2[j]*PySyNVSDE(S,w,sigma2[j],lambda1[i]);
			}
		}
		for(int i=0;i<n1;i++){
			for(int j=0;j<n2;j++){
				sum += pi1[m1+i]*pi2[m1+j]*PySyDE(S,w,lambda1[i],lambda2[j]);
			}
		}
		return sum;
	}
	
	
	/**
	 * Return the Psi(x) defined in SASP-WG/WHL/11-WP/05.
	 * @param x input
	 * @return Psi(x)
	 */
	private static double Psi(double x) throws MathException{
		double val = 1;
		double tmp = Math.abs(x)/Math.sqrt(2);
		if(tmp < 20){ //1-Erf(20)=5.40E-176. It is almost 1.
			try{
				val = Erf.erf(tmp);
			}
			catch(MaxIterationsExceededException ex){
				val = 1;
			}
		}
		return 0.5+Math.signum(x)*0.5*val;
	}
	
	/**
	 * Calculate Py(S) in the case where the TSE distribution follows N-DE distribution.
	 * Analytic solution given in SASP-WG/WHL/11-WP/05 is utilized.
	 * @param S route spacing
	 * @param w average wing span
	 * @param alpha weight coefficient
	 * @param sigma scale parameter of Gaussian distribution
	 * @param lambda scale parameter of double exponential distribution
	 * @return Py(S)
	 */
	public static double PySyNDE(double S, double w, double alpha, double sigma, double lambda) throws MathException{
//		System.out.println("S\t" + S);
//		System.out.println("w\t" + w);
//		System.out.println("alpha\t" + alpha);
//		System.out.println("sigma\t" + sigma);
//		System.out.println("lambda\t" + lambda);
//		System.out.println(Math.pow(1-alpha,2)*(Psi((w-S)/(Math.sqrt(2)*sigma))-Psi((-w-S)/(Math.sqrt(2)*sigma))));
//		System.out.println(alpha*(1-alpha)*(
//					Math.exp(Math.pow(sigma, 2)/(2*Math.pow(lambda, 2)))
//					*(Math.exp((w-S)/lambda)*Psi(-(w-S)/sigma-sigma/lambda)-Math.exp(-(w+S)/lambda)*Psi((w+S)/sigma-sigma/lambda))
//					+Psi((w-S)/sigma)-Psi(-(w+S)/sigma)
//			));
//		System.out.println(alpha*(1-alpha)*(
//				Math.exp(Math.pow(sigma, 2)/(2*Math.pow(lambda, 2)))
//				*(Math.exp((S+w)/lambda)*Psi(-(S+w)/sigma-sigma/lambda)-Math.exp((S-w)/lambda)*Psi(-(S-w)/sigma-sigma/lambda))
//				+Psi((S+w)/sigma)-Psi((S-w)/sigma)
//		));
//		System.out.println( ((S<=w) ?
//				Math.pow(alpha, 2)+(Math.pow(alpha, 2)*Math.exp(-w/lambda)/(2*lambda))*(S*Math.sinh(S/lambda)-(w+2*lambda)*Math.cosh(S/lambda))
//				: (Math.pow(alpha, 2)*Math.exp(-S/lambda)/(2*lambda))*((S+2*lambda)*Math.sinh(w/lambda)-w*Math.cosh(w/lambda))
//			));
		return Math.pow(1-alpha,2)*(Psi((w-S)/(Math.sqrt(2)*sigma))-Psi((-w-S)/(Math.sqrt(2)*sigma)))
			+alpha*(1-alpha)*(
					Math.exp(Math.pow(sigma, 2)/(2*Math.pow(lambda, 2)))
					*(Math.exp((w-S)/lambda)*Psi(-(w-S)/sigma-sigma/lambda)-Math.exp(-(w+S)/lambda)*Psi((w+S)/sigma-sigma/lambda))
					+Psi((w-S)/sigma)-Psi(-(w+S)/sigma)
			)
			+alpha*(1-alpha)*(
					Math.exp(Math.pow(sigma, 2)/(2*Math.pow(lambda, 2)))
					*(Math.exp((S+w)/lambda)*Psi(-(S+w)/sigma-sigma/lambda)-Math.exp((S-w)/lambda)*Psi(-(S-w)/sigma-sigma/lambda))
					+Psi((S+w)/sigma)-Psi((S-w)/sigma)
			)
			+ ((S<=w) ?
				Math.pow(alpha, 2)+(Math.pow(alpha, 2)*Math.exp(-w/lambda)/(2*lambda))*(S*Math.sinh(S/lambda)-(w+2*lambda)*Math.cosh(S/lambda))
				: (Math.pow(alpha, 2)*Math.exp(-S/lambda)/(2*lambda))*((S+2*lambda)*Math.sinh(w/lambda)-w*Math.cosh(w/lambda))
			);
	}
	
	/**
	 * Calculate Py(S) in the case where the TSE distribution follows N-DE distribution.
	 * Approximation given in SASP-WG/WHL/11-WP/05 is utilized.
	 * @param S route spacing
	 * @param w average wing span
	 * @param alpha weight coefficient
	 * @param sigma scale parameter of Gaussian distribution
	 * @param lambda scale parameter of double exponential distribution
	 * @return Py(S)
	 */
	public static double PySyNDEApprox(double S, double w, double alpha, double sigma, double lambda) throws MathException{
		double tmp =  Math.pow(1-alpha,2)*Math.exp(-Math.pow(S, 2)/(4*Math.pow(sigma, 2)))/(2*sigma*Math.sqrt(Math.PI))
				+(alpha*(1-alpha)*Math.exp(Math.pow(sigma,2)/(2*Math.pow(lambda,2)))/lambda)
				*(Math.exp(-S/lambda)*Psi(S/sigma-sigma/lambda)+Math.exp(S/lambda)*Psi(-S/sigma-sigma/lambda))
				+Math.pow(alpha,2)*(lambda+S)*Math.exp(-S/lambda)/(4*Math.pow(lambda,2));
		return 2*w*tmp;
	}


	/**
	 * Calculate Py(S) in the case where the TSE distributions follow different Gaussian distributions.
	 * @param S route spacing
	 * @param w average wing span
	 * @param sigma1 scale parameter #1 of Gaussian distribution
	 * @param sigma2 scale parameter #2 of Gaussian distribution
	 * @return Py(S)
	 */
	public static double PySyN(double S, double w, double sigma1, double sigma2) throws MathException{
		return integralN(S,w,Math.sqrt(Math.pow(sigma1,2)+Math.pow(sigma2,2)));
	}

	/**
	 * Calculate Py(S) in the case where the TSE distribution follows NN distribution.
	 * @param S route spacing
	 * @param w average wing span
	 * @param alpha weight coefficient
	 * @param sigma1 scale parameter #1 of Gaussian distribution
	 * @param sigma2 scale parameter #2 of Gaussian distribution
	 * @return Py(S)
	 */
	public static double PySyNN(double S, double w, double alpha, double sigma1, double sigma2) throws MathException{
		return Math.pow(1-alpha, 2)*integralN(S,w,Math.sqrt(2)*sigma1)
			+2*alpha*(1-alpha)*integralN(S,w,Math.sqrt(Math.pow(sigma1,2)+Math.pow(sigma2,2)))
			+Math.pow(alpha, 2)*integralN(S,w,Math.sqrt(2)*sigma2);
	}

	/**
	 * Calculate Py(S) in the case where the TSE distribution follows a single Gaussian distribution.
	 * @param S route spacing
	 * @param w average wing span
	 * @param sigma scale parameter of Gaussian distribution
	 * @return Py(S)
	 */
	public static double PySyN(double S, double w, double sigma) throws MathException{
		return integralN(S,w,Math.sqrt(2)*sigma);
	}
	/**
	 * Calculate Py(S) in the case where the TSE distribution follows DDE distribution.
	 * (SASP-WG/WHL/9-WP/14)
	 * @param S route spacing
	 * @param w average wing span
	 * @param alpha weight coefficient
	 * @param lambda1 scale parameter #1 of DE distribution
	 * @param lambda2 scale parameter #2 of DE distribution
	 * @return Py(S)
	 */
	public static double PySyDEDE(double S, double w, double alpha, double lambda1, double lambda2){
		return Math.pow(1-alpha,2)*PySyDE(S,w,lambda1,lambda1)+2*alpha*(1-alpha)*PySyDE(S,w,lambda1,lambda2)
			+Math.pow(alpha,2)*PySyDE(S,w,lambda2,lambda2);
	}
	/**
	 * Calculate Py(S) in the case where the TSE distributions follow different DE distributions.
	 * (SASP-WG/WHL/9-WP/14)
	 * @param S route spacing
	 * @param w average wing span
	 * @param lambda1 scale parameter #1 of DE distribution
	 * @param lambda2 scale parameter #2 of DE distribution
	 * @return Py(S)
	 */
	public static double PySyDE(double S, double w, double lambda1, double lambda2){
		if(Math.abs(lambda1-lambda2)<epsilon){
			return 0.5*integralDE(S,w,lambda1,0)+0.5*integralDE(S,w,lambda1,1);
		}
		else{
			return Math.pow(lambda1,2)/(Math.pow(lambda1,2)-Math.pow(lambda2,2))*integralDE(S,w,lambda1,0)
				+ Math.pow(lambda2,2)/(Math.pow(lambda2,2)-Math.pow(lambda1,2))*integralDE(S,w,lambda2,0);
		}
	}

	/**
	 * Calculate Py(S) in the case where the TSE distribution follows a single DE distribution.
	 * (SASP-WG/WHL/9-WP/14)
	 * @param S route spacing
	 * @param w average wing span
	 * @param lambda scale parameter of DE distribution
	 * @return Py(S)
	 */
	public static double PySyDE(double S, double w, double lambda){
		return PySyDE(S,w,lambda,lambda);
	}
	/**
	 * Integrate DE^n(x) functions defined in SASP-WG/WHL/9-WP/14 on [S-w,S+w].
	 * (SASP-WG/WHL/9-WP/14)
	 * @param S route spacing
	 * @param w average wing span
	 * @param lambda scale parameter of DE^n distribution
	 * @param n degree of DE^n distribution
	 * @return integration
	 */
	private static double integralDE(double S, double w, double lambda, int n){
		if(S>=w){
			return integralDE(S+w,lambda,n)-integralDE(S-w,lambda,n);
		}
		else{
			return integralDE(S+w,lambda,n)+integralDE(w-S,lambda,n);
		}
	}
	/**
	 * Integrate DE^n(x) functions defined in SASP-WG/WHL/9-WP/14 on [0,x].
	 * It is assumed that x>0.
	 * (SASP-WG/WHL/9-WP/14)
	 * @param x input
	 * @param lambda scale parameter of DE^n distribution
	 * @param n degree of DE^n distribution
	 * @return integration
	 */
	private static double integralDE(double x, double lambda, int n){
		if(n>0){
			double den = 2*factorial(n)*Math.pow(lambda,n);
			return (Math.pow(0,n)*Math.exp(-0/lambda)/den-Math.pow(x,n)*Math.exp(-x/lambda)/den)
			+integralDE(x,lambda,n-1);
		}
		else{
			return 0.5-0.5*Math.exp(-x/lambda);
		}
	}
	/**
	 * Calculate the factorial of non-negative integers.
	 * @param n non-negative integer
	 * @return factorial of n
	 */
	private static double factorial(int n){
		double tmp = 1;
		for(int i=1;i<=n;i++) tmp*= i;
		return tmp;
	}

	/**
	 * Calculate Pr(S-w<Y<S+w) when Y follows a symmetric Gaussian distribution whose standard deviation is sigma.
	 * @param S route spacing
	 * @param w wing span
	 * @param sigma standard deviation
	 * @return Pr(S-w<Y<S+w)
	 */
	private static double integralN(double S, double w, double sigma) throws MathException{
		return Psi((S+w)/sigma) - Psi((S-w)/sigma);
	}

	// 2011/11/22 The following functions are added.
	/**
	 * Calculate Py(S) in the case where the TSE distribution follows N distribution and DE distribution.
	 * Analytic solution given in SASP-WG/WHL/11-WP/05 is utilized.
	 * @param S route spacing
	 * @param w average wing span
	 * @param sigma scale parameter of Gaussian distribution
	 * @param lambda scale parameter of double exponential distribution
	 * @return Py(S)
	 */
	public static double PySyNVSDE(double S, double w, double sigma, double lambda) throws MathException{
		// cut and past from PySyNDE
		double tmp1 = Math.exp(Math.pow(sigma, 2)/(2*Math.pow(lambda, 2)))
			*(Math.exp((w-S)/lambda)*Psi(-(w-S)/sigma-sigma/lambda)-Math.exp(-(w+S)/lambda)*Psi((w+S)/sigma-sigma/lambda))
			+Psi((w-S)/sigma)-Psi(-(w+S)/sigma);
		double tmp2 = Math.exp(Math.pow(sigma, 2)/(2*Math.pow(lambda, 2)))
			*(Math.exp((S+w)/lambda)*Psi(-(S+w)/sigma-sigma/lambda)-Math.exp((S-w)/lambda)*Psi(-(S-w)/sigma-sigma/lambda))
			+Psi((S+w)/sigma)-Psi((S-w)/sigma);
		return tmp1+tmp2;
	}
	
	/**
	 * Calculate Py(S) in the case where the TSE distribution follows N distribution and DE distribution.
	 * Analytic solution given in SASP-WG/WHL/11-WP/05 is utilized.
	 * @param S route spacing
	 * @param w average wing span
	 * @param sigma scale parameter of Gaussian distribution
	 * @param lambda scale parameter of double exponential distribution
	 * @return Py(S)
	 */
	public static double PySyNVSDEAppro(double S, double w, double sigma, double lambda) throws MathException{
		// cut and past from PySyNDE
		double tmp = (Math.exp(Math.pow(sigma,2)/(2*Math.pow(lambda,2)))/lambda)
			*(Math.exp(-S/lambda)*Psi(S/sigma-sigma/lambda)+Math.exp(S/lambda)*Psi(-S/sigma-sigma/lambda));
		return 2*w*tmp;
	}
	
}
