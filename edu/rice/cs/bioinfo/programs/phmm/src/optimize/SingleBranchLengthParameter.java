/**
 * Add in storage for specific nodes.
 * No ability to multiplex between a single length parameter
 * and multiple branches in this case.
 *
 * Only use this for PhyloNet TNode objects.
 */

package optimize;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;

public class SingleBranchLengthParameter extends LengthParameter {
    protected TNode node;

    /**
     * User must explicitly request whether or not an update is requested 
     * during construction.
     */
    public SingleBranchLengthParameter (String inName, double inValue, TNode inNode, boolean updateFlag) {
	// order forced by Java language constraints
	super(inName, inValue);

	this.node = inNode;

	if (updateFlag) {
	    updateModelState();
	}
    }

    /**
     * WARNING - doesn't update model state by default!
     */
    public void setValue (double inValue) {
	setValue(inValue, false);
    }

    public void setValue (double inValue, boolean updateFlag) {
	super.setValue(inValue);
	if (updateFlag) {
	    updateModelState();
	}
    }

    public void updateModelState () {
	this.node.setParentDistance(getValue());
    }
}
