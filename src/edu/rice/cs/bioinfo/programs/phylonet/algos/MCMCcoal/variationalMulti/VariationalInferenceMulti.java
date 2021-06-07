package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variationalMulti;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HmmBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HmmCore;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.Prior;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.VariationalVariable;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import javax.security.sasl.RealmCallback;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Variational inference procedure for multivariate gaussian posterior.
 *
 * Created by Xinhao Liu on 7/30/20.
 */
public class VariationalInferenceMulti {
    private VariationalModelMulti variationalPosterior;
    private Prior prior;

    public VariationalInferenceMulti(ModelTree model, Prior prior) {
        this.variationalPosterior = new VariationalModelMulti(model);
        this.prior = prior;
    }

    public void run() {
        for (int i = 0; i < Utils.nIterations; i++) {
            System.out.println("=====================================ITERATION " + i +"=====================================");

            List<double[]> samples = new ArrayList<>();
            for (int s = 0; s < Utils.nSamples; s++) {
                samples.add(variationalPosterior.sample());
            }

            System.out.println("Finished sampling");

            // evaluate coalhmm likelihood for each sample
            List<Double> logJoints = new ArrayList<>();
            for (int s = 0; s < samples.size(); s++) {
                double[] sample = samples.get(s);
                variationalPosterior.setTreeBySample(sample);
                ModelTree sampledVariationalPosterior = variationalPosterior.getModel();
                // code below is to circumvent a strange bug causing likelihood to return positive value.
                // not sure why. maybe beagle screwed up.
                double likelihood = 1;
                while (likelihood >= 0) {
                    HmmBuilder builder = new HmmBuilder(sampledVariationalPosterior.getTree(), sampledVariationalPosterior.getRecombRate());
                    // log time used by builder.build(), add up to time used by building HMM by simulation
                    long buildingStartTime = System.currentTimeMillis();
                    HmmCore hmm = builder.build();
                    Utils.buildingTime += System.currentTimeMillis() - buildingStartTime;
                    // log time used by hmm.logLikelihood(), add up to time used by forward algorithm
                    long likelihoodStartTime = System.currentTimeMillis();
                    likelihood = hmm.logLikelihood();
                    Utils.likelihoodTime += System.currentTimeMillis() - likelihoodStartTime;
                }
                double prior = this.prior.logPrior(sampledVariationalPosterior);
                logJoints.add(likelihood + prior);
            }

            System.out.println("Finished evaluating likelihood");

            // Mu
            System.out.println("Current mu is: " + Arrays.toString(variationalPosterior.getDistribution().getMuVector()));
            List<RealVector> fMu = new ArrayList<>();
            List<RealVector> hMu = new ArrayList<>();
            for (int s = 0; s < samples.size(); s++) {
                double[] sample = samples.get(s);
                RealVector fMu_s = variationalPosterior.scoreMu(sample).mapMultiply(logJoints.get(s) - variationalPosterior.logDensity(sample));
                RealVector hMu_s = variationalPosterior.scoreMu(sample);
                fMu.add(fMu_s);
                hMu.add(hMu_s);
            }
            RealVector aMu = estimateVectorCovariance(fMu, hMu).ebeDivide(estimateVectorVariance(hMu));

            RealVector sumMu = MatrixUtils.createRealVector(new double[variationalPosterior.parameterCount]);
            for (int s = 0; s < samples.size(); s++) {
               sumMu = sumMu.add(fMu.get(s).subtract(aMu.ebeMultiply(hMu.get(s))));
            }
            RealVector gradMu = sumMu.mapDivide(samples.size());
            System.out.println("gradMu:");
            System.out.println(Arrays.toString(gradMu.toArray()));

            // Sigma
            System.out.println("Current sigma is: " + Arrays.deepToString(variationalPosterior.getDistribution().getHalfSigmaMatrix()));
            List<RealMatrix> fSigma = new ArrayList<>();
            List<RealMatrix> hSigma = new ArrayList<>();
            for (int s = 0; s < samples.size(); s++) {
                double[] sample = samples.get(s);
                RealMatrix fSigma_s = variationalPosterior.scoreSigma(sample).scalarMultiply(logJoints.get(s) - variationalPosterior.logDensity(sample));
                RealMatrix hSigma_s = variationalPosterior.scoreSigma(sample);
                fSigma.add(fSigma_s);
                hSigma.add(hSigma_s);
            }
            RealMatrix matrixCov = estimateMatrixCovariance(fSigma, hSigma);
            RealMatrix matrixVar = estimateMatrixVariance(hSigma);
            for (int row = 0; row < matrixCov.getRowDimension(); row++) {
                for (int col = 0; col < matrixCov.getColumnDimension(); col++) {
                    matrixCov.multiplyEntry(row, col, 1.0 / matrixVar.getEntry(row, col));
                }
            }
            RealMatrix aSigma = matrixCov;

            RealMatrix sumSigma = MatrixUtils.createRealMatrix(variationalPosterior.parameterCount, variationalPosterior.parameterCount);
            for (int s = 0; s < samples.size(); s++) {
                sumSigma = sumSigma.add(fSigma.get(s).subtract(Utils.hadamardProduct(aSigma, hSigma.get(s))));
            }
            RealMatrix gradSigma = sumSigma.scalarMultiply(1.0 / samples.size());
            System.out.println("gradSigma:");
            System.out.println(Arrays.deepToString(gradSigma.getData()));

            if (Utils.ILLEGAL_SAMPLE_GENERATED) {
                System.out.println("!!!!!!!!!!!!!!!!!!!!!!ILLEGAL SAMPLE GENERATED!!!!!!!!!!!!!!!!!!!!!!");
            }
            variationalPosterior.muGradientUpdate(gradMu, i);
            variationalPosterior.sigmaGradientUpdate(gradSigma, i);
            // flip Utils.ILLEGAL_SAMPLE_GENERATED
            if (Utils.ILLEGAL_SAMPLE_GENERATED) {
                Utils.ILLEGAL_SAMPLE_GENERATED = false;
            }
        }
    }

    /**
     * Calculate sample covariance for each row of the vector separately.
     */
    private static RealVector estimateVectorCovariance(List<RealVector> x, List<RealVector> y) {
        List<Tuple<List<Double>, List<Double>>> rows = new ArrayList<>();
        for (int rIdx = 0; rIdx < x.get(0).getDimension(); rIdx++) {
            List<Double> xAtThisRow = new ArrayList<>();
            List<Double> yAtThisRow = new ArrayList<>();
            for (int cIdx = 0; cIdx < x.size(); cIdx++) {
                xAtThisRow.add(x.get(cIdx).getEntry(rIdx));
                yAtThisRow.add(y.get(cIdx).getEntry(rIdx));
            }
            rows.add(new Tuple<>(xAtThisRow, yAtThisRow));
        }

        RealVector rowWiseCovVector = MatrixUtils.createRealVector(new double[]{});
        for (Tuple<List<Double>, List<Double>> row : rows) {
            rowWiseCovVector = rowWiseCovVector.append(estimateCovariance(row.Item1, row.Item2));
        }
        return rowWiseCovVector;
    }

    /**
     * Calculate sample covariance for each entry of the matrix separately.
     */
    private static RealMatrix estimateMatrixCovariance(List<RealMatrix> xMatrices, List<RealMatrix> yMatrices) {
        assert xMatrices.get(0).getRowDimension() == xMatrices.get(0).getColumnDimension() && yMatrices.get(0).getRowDimension() == yMatrices.get(0).getColumnDimension() && xMatrices.get(0).getRowDimension() == yMatrices.get(0).getColumnDimension();
        int matrixSize = xMatrices.get(0).getRowDimension();
        Tuple<ArrayList<Double>, ArrayList<Double>>[][] aggregateMatrix = new Tuple[matrixSize][matrixSize];
        for (int i = 0; i < matrixSize; i++) {
            for (int j = 0; j < matrixSize; j++) {
                aggregateMatrix[i][j] = new Tuple<>(new ArrayList<>(), new ArrayList<>());
            }
        }

        for (int idx = 0; idx < xMatrices.size(); idx++) {
            for (int i = 0; i < matrixSize; i++) {
                for (int j = 0; j < matrixSize; j++) {
                    aggregateMatrix[i][j].Item1.add(xMatrices.get(idx).getEntry(i, j));
                    aggregateMatrix[i][j].Item2.add(yMatrices.get(idx).getEntry(i, j));
                }
            }
        }

        RealMatrix entryWiseCovMatrix = MatrixUtils.createRealMatrix(new double[matrixSize][matrixSize]);
        for (int i = 0; i < matrixSize; i++) {
            for (int j = 0; j < matrixSize; j++) {
                entryWiseCovMatrix.setEntry(i, j, estimateCovariance(aggregateMatrix[i][j].Item1, aggregateMatrix[i][j].Item2));
            }
        }
        return entryWiseCovMatrix;
    }

    /**
     * Calculate sample covariance given two list of numbers.
     */
    private static double estimateCovariance(List<Double> x, List<Double> y) {
        double meanX = x.stream().mapToDouble(a -> a).average().getAsDouble();
        double meanY = y.stream().mapToDouble(a -> a).average().getAsDouble();

        double sum = 0.0;
        for (int i = 0; i < x.size(); i++) {
            sum += (x.get(i) - meanX) * (y.get(i) - meanY);
        }

        return sum / (x.size() - 1);
    }

    /**
     * Calculate sample variance for each row of the vector separately.
     */
    private static RealVector estimateVectorVariance(List<RealVector> x) {
        List<List<Double>> rows = new ArrayList<>();
        for (int rIdx = 0; rIdx < x.get(0).getDimension(); rIdx++) {
            List<Double> row = new ArrayList<>();
            for (int cIdx = 0; cIdx < x.size(); cIdx++) {
                row.add(x.get(cIdx).getEntry(rIdx));
            }
            rows.add(row);
        }

        RealVector rowWiseVarVector = MatrixUtils.createRealVector(new double[]{});
        for (List<Double> row : rows) {
            rowWiseVarVector = rowWiseVarVector.append(estimateVariance(row));
        }
        return rowWiseVarVector;
    }

    /**
     * Calculate sample variance for each entry of the matrix separately.
     */
    private static RealMatrix estimateMatrixVariance(List<RealMatrix> matrices) {
        assert matrices.get(0).getRowDimension() == matrices.get(0).getColumnDimension();
        int matrixSize = matrices.get(0).getRowDimension();
        ArrayList<Double>[][] aggregateMatrix = new ArrayList[matrixSize][matrixSize];
        for (int i = 0; i < matrixSize; i++) {
            for (int j = 0; j < matrixSize; j++) {
                aggregateMatrix[i][j] = new ArrayList<>();
            }
        }
        for (RealMatrix matrix : matrices) {
            for (int i = 0; i < matrixSize; i++) {
                for (int j = 0; j < matrixSize; j++) {
                    aggregateMatrix[i][j].add(matrix.getEntry(i, j));
                }
            }
        }

        RealMatrix entryWiseVarMatrix = MatrixUtils.createRealMatrix(new double[matrixSize][matrixSize]);
        for (int i = 0; i < matrixSize; i++) {
            for (int j = 0; j < matrixSize; j++) {
                entryWiseVarMatrix.setEntry(i, j, estimateVariance(aggregateMatrix[i][j]));
            }
        }
        return entryWiseVarMatrix;
    }

    /**
     * Calculate sample variance for a list of numbers
     */
    private static double estimateVariance(List<Double> x) {
        double meanX = x.stream().mapToDouble(a -> a).average().getAsDouble();

        double sum = 0.0;
        for (double num:x) {
            sum += (num - meanX) * (num - meanX);
        }

        return sum / (x.size() - 1);
    }

    public VariationalModelMulti getVariationalPosterior() {
        return variationalPosterior;
    }

}
