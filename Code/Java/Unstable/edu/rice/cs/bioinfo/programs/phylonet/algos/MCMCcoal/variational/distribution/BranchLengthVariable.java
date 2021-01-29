package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution;


import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.TreeNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;

/**
 * Variational distribution on branch length variables for each necessary branch.
 *
 * Created by Xinhao Liu on 1/24/21.
 */
public class BranchLengthVariable extends VariationalVariable {

    public STINode<TreeNodeInfo> source;
    public STINode<TreeNodeInfo> dest;

    private double rMean; // accumulated squared gradient for mean
    private double rStdDev; // accumulated squared gradient for standard deviation

    private final double deltaNumericalStability = 1E-7;

    public BranchLengthVariable(double mean, double standardDeviation, STINode<TreeNodeInfo> source, STINode<TreeNodeInfo> dest) {
        super(mean, standardDeviation);
        this.source = source;
        this.dest = dest;
        this.rMean = 0;
        this.rStdDev = 0;
    }

    @Override
    public void setVariableValue(double value) {
        // assert dest.getParent() == source; // should be the same object
        dest.setParentDistance(value);
        // TODO: need to set source's node height? the order of setting node heights is a problem, but postorder should be good
    }

    @Override
    public void meanGradientUpdate(double gradient) {
        if (!Utils.ILLEGAL_SAMPLE_GENERATED)
            rMean += gradient * gradient;
        System.out.println("Gradient of branch length mean is: " + (gradient / (deltaNumericalStability + Math.sqrt(rMean))));
        double update = (Utils.BRANCH_LENGTH_MEAN_LEARNING_RATE / (deltaNumericalStability + Math.sqrt(rMean))) * gradient;
        System.out.println("Update of branch length mean is: " + update);
        //setMean(getMean() + Utils.BRANCH_LENGTH_MEAN_LEARNING_RATE * gradient);
        setMean(getMean() + update);
    }

    @Override
    public void standardDeviationGradientUpdate(double gradient) {
        if (!Utils.ILLEGAL_SAMPLE_GENERATED)
            rStdDev += gradient * gradient;
        System.out.println("Gradient of branch length standard deviation is: " + (gradient / (deltaNumericalStability + Math.sqrt(rStdDev))));
        double update = (Utils.BRANCH_LENGTH_STDDEV_LEARNING_RATE / (deltaNumericalStability + Math.sqrt(rStdDev))) * gradient;
        System.out.println("Update of branch length standard deviation is: " + update);
        double newStandardDeviation = getStandardDeviation() + update;

        if (newStandardDeviation < Utils.BRANCH_LENGTH_MIN_STDDEV) {
            return;
        }
        setStandardDeviation(newStandardDeviation);
    }
}
