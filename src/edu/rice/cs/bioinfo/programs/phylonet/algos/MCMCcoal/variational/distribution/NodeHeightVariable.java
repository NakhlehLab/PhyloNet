package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.TreeNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * Variational distribution on node height variables for root and each internal node.
 *
 * Created by Xinhao Liu on 3/14/20.
 */
public class NodeHeightVariable extends VariationalVariable {

    public STINode<TreeNodeInfo> node;

    private double rMean; // accumulated squared gradient for mean
    private double rStdDev; // accumulated squared gradient for standard deviation

    private final double deltaNumericalStability = 1E-7;

    public NodeHeightVariable(double mean, double standardDeviation, STINode<TreeNodeInfo> node) {
        super(mean, standardDeviation);
        this.node = node;
        this.rMean = 0;
        this.rStdDev = 0;
    }

    // for rerunning purpose
    public NodeHeightVariable(double mean, double standardDeviation, STINode<TreeNodeInfo> node, double rMean, double rStdDev) {
        super(mean, standardDeviation);
        this.node = node;
        this.rMean = rMean;
        this.rStdDev = rStdDev;
    }

    @Override
    public void setVariableValue(double value) {
        node.setNodeHeight(value);
    }

    @Override
    public void meanGradientUpdate(double gradient) {
        if (!Utils.ILLEGAL_SAMPLE_GENERATED)
            rMean += gradient * gradient;
        System.out.println("Gradient of node height mean is: " + (gradient / (deltaNumericalStability + Math.sqrt(rMean))));
        double update = (Utils.NODE_HEIGHT_MEAN_LEARNING_RATE / (deltaNumericalStability + Math.sqrt(rMean))) * gradient;
        System.out.println("Update of node height mean is: " + update);
        //setMean(getMean() + Utils.NODE_HEIGHT_MEAN_LEARNING_RATE * gradient);
        setMean(getMean() + update);

//        if (Utils.lastiter) {
//            System.out.println("rMean for node height mean is: " + rMean);
//        }

    }

    @Override
    public void standardDeviationGradientUpdate(double gradient) {
        if (!Utils.ILLEGAL_SAMPLE_GENERATED)
            rStdDev += gradient * gradient;
        System.out.println("Gradient of node height standard deviation is: " + (gradient / (deltaNumericalStability + Math.sqrt(rStdDev))));
        double update = (Utils.NODE_HEIGHT_STDDEV_LEARNING_RATE / (deltaNumericalStability + Math.sqrt(rStdDev))) * gradient;
        System.out.println("Update of node height standard deviation is: " + update);
        double newStandardDeviation = getStandardDeviation() + update;

//        if (Utils.lastiter) {
//            System.out.println("rStdDev for node height stddev is: " + rStdDev);
//        }

//        double newStandardDeviation = getStandardDeviation() + Utils.NODE_HEIGHT_STDDEV_LEARNING_RATE * gradient;
        if (newStandardDeviation < Utils.NODE_HEIGHT_MIN_STDDEV) {
            return;
        }
//        setStandardDeviation(getStandardDeviation() + Utils.NODE_HEIGHT_STDDEV_LEARNING_RATE * gradient);
        setStandardDeviation(newStandardDeviation);
    }
}
