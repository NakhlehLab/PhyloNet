package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variationalMulti;


import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.Arrays;

/**
 * Wrapper of the multivariate normal distribution used as a variational posterior of the entire model. i.e. not factorized.
 *
 * Created by Xinhao Liu on 7/27/20.
 */
public class MultiVariateGaussian {
    private MultivariateNormalDistribution distribution;
    private double[] muVector;
    private double[][] halfSigmaMatrix; // upper triangular matrix

    /**
     * Complete an upper triangular matrix
     * @param halfMatrix upper triangular matrix
     * @return a symmetric matrix
     */
    public static double[][] full(double[][] halfMatrix) {
        double[][] fullMatrix = Arrays.stream(halfMatrix).map(double[]::clone).toArray(double[][]::new);
        for (int i = 0; i < halfMatrix.length; i++) {
            for (int j = 0; j < halfMatrix.length; j++) {
                if (i > j) {
                    fullMatrix[i][j] = halfMatrix[j][i];
                }
            }
        }
        return fullMatrix;
    }

    /**
     * @param muVector
     * @param halfSigmaMatrix upper triangular matrix
     */
    public MultiVariateGaussian (double[] muVector, double[][] halfSigmaMatrix) {
        this.muVector = muVector;
        this.halfSigmaMatrix = halfSigmaMatrix;
        distribution = new MultivariateNormalDistribution(muVector, full(halfSigmaMatrix));
    }

    public void setMu(RealVector newMu) {
        muVector = newMu.toArray();
        distribution = new MultivariateNormalDistribution(muVector, full(halfSigmaMatrix));
    }

    public void setSigma(RealMatrix newSigma) {
        halfSigmaMatrix = newSigma.getData();
        System.out.println("In setSigma -- current mean is: " + Arrays.toString(muVector));
        System.out.println("In setSigma -- current sigma is: " + Arrays.deepToString(full(halfSigmaMatrix)));
        distribution = new MultivariateNormalDistribution(muVector, full(halfSigmaMatrix));
    }

    public double logDensity(double[] x) {
        return StrictMath.log(distribution.density(x));
    }

    public double[] sample() {
        return distribution.sample();
    }

    public RealVector scoreMu(double[] x) {
        RealMatrix sigma = MatrixUtils.createRealMatrix(full(halfSigmaMatrix));
        RealMatrix sigmaInverse = MatrixUtils.inverse(sigma);
        RealVector xVec = MatrixUtils.createRealVector(x);
        RealVector muVec = MatrixUtils.createRealVector(muVector);
        return sigmaInverse.operate(xVec.subtract(muVec));
    }

    public RealMatrix scoreSigma(double[] x) {
        RealVector xVec = MatrixUtils.createRealVector(x);
        RealVector muVec = MatrixUtils.createRealVector(muVector);

        RealMatrix sigmaInverseTranspose = MatrixUtils.inverse(MatrixUtils.createRealMatrix(full(halfSigmaMatrix))).transpose();
        RealMatrix firstTerm = sigmaInverseTranspose.scalarMultiply(-0.5);
        RealMatrix secondTerm = sigmaInverseTranspose.multiply(xVec.subtract(muVec).outerProduct(xVec.subtract(muVec))).multiply(sigmaInverseTranspose).scalarMultiply(0.5);
        return firstTerm.add(secondTerm);
    }

    public double[] getMuVector() {
        return muVector;
    }

    public double[][] getHalfSigmaMatrix() {
        return halfSigmaMatrix;
    }

    public static void main(String[] args) {
        double[][] halfMatrix = new double[][] { { 0,  2, 3, 4},
                                                 { -1, 0, 7, 8},
                                                 { -1, -1, 0, 12},
                                                 { -1, -1, -1, 0}};
        double[][] completed = full(halfMatrix);
        for (int i = 0; i < completed.length; i++) {
            for (int j = 0; j < completed.length; j++) {
                if (completed[i][j] != completed[j][i]) {
                    System.out.println("ERROR!!!!!!");
                }
            }
        }
        System.out.println(Arrays.deepToString(completed));
    }
}
