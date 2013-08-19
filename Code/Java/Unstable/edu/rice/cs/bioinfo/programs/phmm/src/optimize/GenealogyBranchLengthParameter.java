/**
 * Add in storage for specific nodes.
 * No ability to multiplex between a single length parameter
 * and multiple branches in this case.
 *
 * Only use this for PhyloNet TNode objects.
 */

package optimize;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;

public class GenealogyBranchLengthParameter extends Parameter {
    public static final double DEFAULT_MINIMUM_BRANCH_LENGTH = MultivariateOptimizer.DEFAULT_MINIMUM_BRANCH_LENGTH;
    public static final double DEFAULT_INITIAL_BRANCH_LENGTH = MultivariateOptimizer.DEFAULT_INITIAL_BRANCH_LENGTH;
    public static final double DEFAULT_MAXIMUM_BRANCH_LENGTH = MultivariateOptimizer.DEFAULT_MAXIMUM_BRANCH_LENGTH;

    protected TNode node;
    protected CalculationCache calculationCache;

    /**
     * User must explicitly request whether or not an update is requested 
     * during construction.
     */
    public GenealogyBranchLengthParameter (String inName, 
					   double inValue, 
					   TNode inNode, 
					   CalculationCache inCalculationCache,
					   boolean checkValueMinimumConstraintFlag,
					   boolean checkValueMaximumConstraintFlag,
					   boolean updateModelStateFlag) {
	// order forced by Java language constraints
	super(inName, inValue, checkValueMinimumConstraintFlag, checkValueMinimumConstraintFlag, false);

	this.node = inNode;
	this.calculationCache = inCalculationCache;

	if (updateModelStateFlag) {
	    updateModelState();
	}
    }

    public void updateModelState () {
	// kliu - is this getting set properly?
	this.node.setParentDistance(getValue());

	// emission probabilities calculated on-the-fly in HMM object
	// for that reason, no need to reference HMM object here
	// 
	// later, however, may want to consider caching here???
	// or just within HiddenState???

	// invalidate caches for associated gene genealogy node and tree
	calculationCache.cacheSubstitutionProbabilityMatrix.remove(node);
	calculationCache.cacheSubstitutionProbability.clear(node.getTree());
    }

    public double getMinimumValue () {
	return (DEFAULT_MINIMUM_BRANCH_LENGTH);
    }

    public double getDefaultInitialValue () {
	return (DEFAULT_INITIAL_BRANCH_LENGTH);
    }

    public double getMaximumValue () {
	return (DEFAULT_MAXIMUM_BRANCH_LENGTH);
    }
}



    /**
     * WARNING - doesn't update model state by default!
     */
    // public void setValue (double inValue) {
    // 	setValue(inValue, false);
    // }

    // public void setValue (double inValue, boolean updateFlag) {
    // 	super.setValue(inValue);
    // 	if (updateFlag) {
    // 	    updateModelState();
    // 	}
    // }
