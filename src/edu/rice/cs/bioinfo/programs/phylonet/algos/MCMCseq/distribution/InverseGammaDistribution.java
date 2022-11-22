package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.distribution;
/*
 * @ClassName:   InverseGammaDistribution
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        11/10/22 9:35 AM
 */


/*
 * InverseGammaDistribution.java
 *
 * BEAST: Bayesian Evolutionary Analysis by Sampling Trees
 * Copyright (C) 2014 BEAST Developers
 *
 * BEAST is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * BEAST is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with BEAST.  If not, see <http://www.gnu.org/licenses/>.
 */


//        import beast.math.GammaFunction;
//        import beast.math.UnivariateFunction;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.distribution.GammaDistribution;

public class InverseGammaDistribution extends AbstractRealDistribution {



    private double shape, scale;

    private final double factor;
    private final double logFactor;

    public InverseGammaDistribution(double shape, double scale) throws NotStrictlyPositiveException {
        this.shape = shape;
        this.scale = scale;
        this.factor = 1.0; // Math.pow(scale, shape) / Math.exp(GammaFunction.lnGamma(shape));
        this.logFactor = 0.0; // shape * Math.log(scale) - GammaFunction.lnGamma(shape);
    }

    public double getShape() {
        return shape;
    }

    public void setShape(double value) {
        shape = value;
    }

    public double getScale() {
        return scale;
    }

    public void setScale(double value) {
        scale = value;
    }

    public double density(double x) {
        return pdf(x, shape, scale, factor);
    }

    public double logPdf(double x) {
        return logPdf(x, shape, scale, logFactor);
    }

    public double cumulativeProbability(double x) {
        return cdf(x, shape, scale);
    }


    public double sample() {
        return sample(shape, scale);
    }

    /**
     * probability density function of the Gamma distribution
     *
     * @param x     argument
     * @param shape shape parameter
     * @param scale scale parameter
     * @return pdf value
     */
    public static double pdf(double x, double shape, double scale, double factor) {
        if (x <= 0)
            return 0.0;

        final double a = Math.exp(logPdf(x, shape, scale, 0.0));

        return factor * a;
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x     argument
     * @param shape shape parameter
     * @param scale scale parameter
     * @return log pdf value
     */
    public static double logPdf(double x, double shape, double scale, double factor) {
        if (x <= 0)
            return Double.NEGATIVE_INFINITY;

        return factor - (scale / x) - (shape+1) * Math.log(x) + shape * Math.log(scale) - Gamma.logGamma(shape);
//        return  factor + shape*Math.log(scale) - (shape + 1)*Math.log(x) - (scale/x) - GammaFunction.lnGamma(shape);
    }

    /**
     * cumulative density function of the Gamma distribution
     *
     * @param x     argument
     * @param shape shape parameter
     * @param scale scale parameter
     * @return cdf value
     */
    public static double cdf(double x, double shape, double scale) {
        if (x <= 0.0 || shape <= 0.0) {
            return 0.0;
        }
//        return GammaFunction.incompleteGammaQ(shape, scale/x);
        return Gamma.regularizedGammaQ(shape, scale/x);
    }


    /**
     * mean of the Gamma distribution
     *
     * @param shape shape parameter
     * @param scale scale parameter
     * @return mean
     */
    public static double mean(double shape, double scale) {
        if( shape > 1 ) {
            return scale / (shape - 1);
        }
        return Double.POSITIVE_INFINITY;
    }

    /**
     * variance of the Gamma distribution.
     *
     * @param shape shape parameter
     * @param scale scale parameter
     * @return variance
     */
    public static double variance(double shape, double scale) {
        if( shape > 2 ) {
            return scale*scale / ((shape - 1)*(scale-1)*(scale-2));
        }
        return Double.POSITIVE_INFINITY;
    }

    /**
     * sample from the Gamma distribution.
     *
     * @param shape shape parameter
     * @param scale scale parameter
     * @return sample
     */
    public static double sample(double shape, double scale) {
        return 1.0 / new GammaDistribution(shape, 1/scale).sample();
    }

    public double getNumericalMean() {
        return mean(shape, scale);

    }

    public double getNumericalVariance() {
        return variance(shape, scale);
    }

    public double getSupportLowerBound() {
        return 0.0D;
    }

    public double getSupportUpperBound() {
        return 1.0D / 0.0;
    }

    public boolean isSupportLowerBoundInclusive() {
        return true;
    }

    public boolean isSupportUpperBoundInclusive() {
        return false;
    }

    public boolean isSupportConnected() {
        return true;
    }

    public static void main(String[] args) {
        InverseGammaDistribution nodeheight = new InverseGammaDistribution(3, 0.2);
        System.out.println(nodeheight.logPdf(0.00001));
//        System.out.println(nodeheight.cumulativeProbability(0.1));
//        System.out.println(nodeheight.getNumericalMean());
//        System.out.println(nodeheight.getNumericalVariance());

    }

}
