package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.distribution;
/*
 * @ClassName:   DirichletDistribution
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        9/3/21 10:14 AM
 */

import org.apache.commons.math3.distribution.AbstractMultivariateRealDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;


public class DirichletDistribution extends AbstractMultivariateRealDistribution {
    private double[] alpha;
    private RealDistribution[] gammas;

    public DirichletDistribution(double[] alpha) {
        this(new Well19937c(), alpha);

    }
    public DirichletDistribution(RandomGenerator rng, double[] alpha) {
        super(rng, alpha.length);
        this.alpha = alpha;
    }


    @Override
    public double[] sample() {
        initialize();
        double[] sample = new double[ alpha.length ];
        double sum = 0.0;
        for (int i = 0; i < sample.length; i++) {
            sample[i] = gammas[i].sample();
            sum += sample[i];
        }
        for (int i = 0; i < sample.length; i++) {
            sample[i] /= sum;
        }
        return sample;
    }

    public double logdensity(double[] x){
        if (x.length != this.alpha.length) {
            throw new IllegalArgumentException("Dirichlet Distribution: Length of observations != length of alpha");
        }
        double logP = 0.0;
        double alpha_sum = 0.0;
        for (int i = 0; i < x.length; i++){
            logP += (this.alpha[i] - 1) * Math.log(x[i]);
            logP += -Math.log(Gamma.gamma(this.alpha[i]));
            alpha_sum += this.alpha[i];
        }
//        System.out.println(alpha_sum);
//        System.out.println("gamma(sum alpha)="+Gamma.gamma(alpha_sum));
//        System.out.println(Gamma.logGamma(alpha_sum));
        logP += Gamma.logGamma(alpha_sum);
        return logP;
    }

    @Override
    public double density(double[] x){
        return FastMath.exp(logdensity(x));
    }

    private void initialize() {
        if (this.gammas == null) {
            this.gammas = new RealDistribution[alpha.length];
            for (int i = 0; i < alpha.length; i++) {
                gammas[i] = new GammaDistribution(alpha[i], 1);
            }
        }
    }

    public static void main(String[] args) {

        double[] alpha = {2, 2, 2};
        double[] quantities = {0.2, 0.2, 0.6};
        DirichletDistribution dirichletDistribution = new DirichletDistribution(alpha);
        double pdf = dirichletDistribution.density(quantities);
        double logpdf = dirichletDistribution.logdensity(quantities);
        System.out.println(pdf);
        System.out.println(Math.log(pdf));
        System.out.println(logpdf);
//        double[] samples = dirichletDistribution.sample();
//        for (int i = 0; i < alpha.length; i++){
//            System.out.println(samples[i]);
//
//        }
    }
}