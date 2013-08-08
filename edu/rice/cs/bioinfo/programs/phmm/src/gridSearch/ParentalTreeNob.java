package gridSearch;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

public class ParentalTreeNob extends Nob {
    // an annoying detail about using Network class to represent trees
    protected NetNode<Double> parentNode;
    private NetNode<Double> childNode;

    public ParentalTreeNob(int gIn, double minIn, double maxIn, NetNode<Double> childNodeIn)  {
        super (gIn, minIn, maxIn);
        this.childNode = childNodeIn;

	if (childNode.getParentNumber() != 1) {
	    throw (new RuntimeException("ERROR: ParentalTreeNob accepts only tree nodes with exactly one parent, even though they're represented using NetNode objects. This excludes root nodes."));
	}

	parentNode = childNode.getParents().iterator().next();
    }

    @Override
    public void set_param(double value) {
        backupParam = childNode.getParentDistance(parentNode);
        childNode.setParentDistance(parentNode, value);
    }

    @Override
    public double get_param() {
        return childNode.getParentDistance(parentNode);
    }

    public void restoreParameterValue() {
        childNode.setParentDistance(parentNode, backupParam);
    }

}
