package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variationalMulti;

import Jama.Matrix;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.TreeNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.Arrays;


/**
 * Variational approximation to the posterior of the entire model, where the variational distribution is a multivariate gaussian.
 *
 * Created by Xinhao Liu on 7/29/20.
 */
public class VariationalModelMulti {
    private ModelTree model;
    private MultiVariateGaussian variationalDistribution;
    public int parameterCount;

    private RealVector rMu;
    private RealMatrix rSigma;
    private RealVector lrMu;
    private RealMatrix lrSigma;

    private final double deltaNumericalStability = 1E-7;

    /**
     * In the mu vector and sigma matrix, node height parameters are the former indices, population size parameters are the latter indices. They are not interleaved.
     * All parameters are in postorder.
     * That is, first n elements of mu are node height means in postorder. Last n elements of mu are population size means in postorder.
     */
    public VariationalModelMulti(ModelTree model) {
        this.model = model;
        STITree<TreeNodeInfo> tree = model.getTree();
        parameterCount = (tree.getNodeCount() - tree.getLeafCount()) * 2;

        double[] initMu = new double[parameterCount];
        double[][] initHalfSigma = new double[parameterCount][parameterCount];
        rMu = MatrixUtils.createRealVector(new double[parameterCount]);
        rSigma = MatrixUtils.createRealMatrix(parameterCount, parameterCount);
        lrMu = MatrixUtils.createRealVector(new double[parameterCount]);
        lrSigma = MatrixUtils.createRealMatrix(parameterCount, parameterCount);


        // Initialize node height means and variances
        int parameterIdx = 0;
        for (TNode n:tree.postTraverse()) {
            if (n.isLeaf()) {
                continue;
            }

            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) n;
            initMu[parameterIdx] = node.getNodeHeight();
            initHalfSigma[parameterIdx][parameterIdx] = Utils.NODE_HEIGHT_INIT_STDDEV * Utils.NODE_HEIGHT_INIT_STDDEV;
            lrMu.setEntry(parameterIdx, Utils.NODE_HEIGHT_MEAN_LEARNING_RATE);
            parameterIdx++;
        }
        // Initialize pop size means and variances
        for (TNode n:tree.postTraverse()) {
            if (n.isLeaf()) {
                continue;
            }

            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) n;
            initMu[parameterIdx] = node.getData().getPopSize();
            initHalfSigma[parameterIdx][parameterIdx] = Utils.POP_SIZE_INIT_STDDEV * Utils.POP_SIZE_INIT_STDDEV;
            lrMu.setEntry(parameterIdx, Utils.POP_SIZE_MEAN_LEARNING_RATE);
            parameterIdx++;
        }

        // Initialize all covariances to be zero. Initialize sigma learning rate as an upper triangular matrix.
        for (int i = 0; i < parameterCount; i++) {
            for (int j = 0; j < parameterCount; j++) {
                if (i < j) {
                    initHalfSigma[i][j] = 0;
                }
                if (i <= j) {
                    lrSigma.setEntry(i, j, Utils.SIGMA_LEARNING_RATE);
                }
            }
        }

        this.variationalDistribution = new MultiVariateGaussian(initMu, initHalfSigma);
    }

    public void muGradientUpdate(RealVector gradMu, int iter) {
        // decrease learning rate in event of illegal sample
        RealVector lrMuThisIter;
        if (Utils.ILLEGAL_SAMPLE_GENERATED || iter == 0) {
            lrMuThisIter = lrMu.mapDivide(1000);
        } else {
            lrMuThisIter = lrMu;
        }

        if (!Utils.ILLEGAL_SAMPLE_GENERATED) {
            rMu = rMu.add(gradMu.ebeMultiply(gradMu));
        }
        // calculate adapted learning rate
        RealVector rMuSqrt = MatrixUtils.createRealVector(new double[]{});
        for (int i = 0; i < parameterCount; i++) {
            rMuSqrt = rMuSqrt.append(Math.sqrt(rMu.getEntry(i)));
        }
//        System.out.println("rMuSqrt:");
//        System.out.println(Arrays.toString(rMuSqrt.toArray()));
        RealVector denominator = rMuSqrt.mapAdd(deltaNumericalStability);  // denominator of AdaGrad adapted learning rate
//        System.out.println("Denominator:");
//        System.out.println(Arrays.toString(denominator.toArray()));
        RealVector adaptedLR = lrMuThisIter.ebeDivide(denominator);
        RealVector update = adaptedLR.ebeMultiply(gradMu); // Update to Mu under AdaGrad
        System.out.println("Mu update: " + Arrays.toString(update.toArray()));
        variationalDistribution.setMu(MatrixUtils.createRealVector(variationalDistribution.getMuVector()).add(update));
    }

    public void sigmaGradientUpdate(RealMatrix gradSigma, int iter) {
        // decrease learning rate in event of illegal sample
        RealMatrix lrSigmaThisIter;
        if (Utils.ILLEGAL_SAMPLE_GENERATED || iter == 0) {
            lrSigmaThisIter = lrSigma.scalarMultiply(1.0 / 1000000); // 1000^2?
        } else {
            lrSigmaThisIter = lrSigma;
        }

        if (!Utils.ILLEGAL_SAMPLE_GENERATED) {
            rSigma = rSigma.add(Utils.hadamardProduct(gradSigma, gradSigma));
        }
        // calculate adapted learning rate
        RealMatrix rSigmaSqrt = MatrixUtils.createRealMatrix(parameterCount, parameterCount);
        for (int i = 0; i < parameterCount; i++) {
            for (int j = 0; j < parameterCount; j++) {
                rSigmaSqrt.setEntry(i, j, Math.sqrt(rSigma.getEntry(i, j)));
            }
        }
        RealMatrix denominator = rSigmaSqrt.scalarAdd(deltaNumericalStability);  // denominator of AdaGrad adapted learning rate
        RealMatrix adaptedLR = MatrixUtils.createRealMatrix(parameterCount, parameterCount);
        for (int i = 0; i < parameterCount; i++) {
            for (int j = 0; j < parameterCount; j++) {
                adaptedLR.setEntry(i, j, lrSigmaThisIter.getEntry(i, j) / denominator.getEntry(i, j)); // calculate per-entry adapted learning rate
            }
        }
        RealMatrix update = Utils.hadamardProduct(adaptedLR, gradSigma);
        System.out.println("Sigma update: " + Arrays.deepToString(update.getData()));
//        System.out.println("New Sigma is: ");
//        System.out.println(Arrays.deepToString(MatrixUtils.createRealMatrix(variationalDistribution.getHalfSigmaMatrix()).add(update).getData()));
        variationalDistribution.setSigma(MatrixUtils.createRealMatrix(variationalDistribution.getHalfSigmaMatrix()).add(update));
    }

    public double[] sample() {
        return variationalDistribution.sample();
    }

    /**
     * Set tree parameters to sample values
     */
    public void setTreeBySample(double[] sample) {
        STITree<TreeNodeInfo> tree = this.model.getTree();

        int parameterIdx = 0;
        for (TNode n:tree.postTraverse()) {
            if (n.isLeaf()) {
                continue;
            }

            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) n;
            node.setNodeHeight(sample[parameterIdx]);
            parameterIdx++;
        }
        for (TNode n:tree.postTraverse()) {
            if (n.isLeaf()) {
                continue;
            }

            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) n;
            node.getData().setPopSize((int) sample[parameterIdx]);
            parameterIdx++;
        }

        this.model.refresh();
    }

    public ModelTree getModel() {
        return model;
    }

    public MultiVariateGaussian getDistribution() {
        return variationalDistribution;
    }

    public RealVector scoreMu(double[] x) {
        return variationalDistribution.scoreMu(x);
    }

    public RealMatrix scoreSigma(double[] x) {
        return variationalDistribution.scoreSigma(x);
    }

    public double logDensity(double[] x) {
        return variationalDistribution.logDensity(x);
    }
}
