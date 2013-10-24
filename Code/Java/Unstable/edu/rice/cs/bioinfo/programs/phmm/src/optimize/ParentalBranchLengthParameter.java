package optimize;

import runHmm.runHmm;
import util.Constants;
import edu.rice.cs.bioinfo.library.programming.BijectiveHashtable;
import edu.rice.cs.bioinfo.library.programming.BidirectionalMultimap;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.library.programming.Tuple;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF.CoalescePattern;
import java.util.Hashtable;

/**
 * Up to caller to associate this ParentalBranchLengthParameter class
 * with a specific branch.
 *
 * No internal enforcement that the first member of a length-parameter-constraint-set is always set to 1.0.
 * Enforcement is performed externally in MultivariateOptimizer.
 */

public class ParentalBranchLengthParameter extends Parameter {
    public static final double DEFAULT_MINIMUM_BRANCH_LENGTH = MultivariateOptimizer.DEFAULT_MINIMUM_BRANCH_LENGTH;
    public static final double DEFAULT_INITIAL_BRANCH_LENGTH = MultivariateOptimizer.DEFAULT_INITIAL_BRANCH_LENGTH;
    public static final double DEFAULT_MAXIMUM_BRANCH_LENGTH = MultivariateOptimizer.DEFAULT_MAXIMUM_BRANCH_LENGTH;

    public static final double EPSILON = Constants.ZERO_DELTA;

    // Use external BidirectionalMap to
    // switch between ParentalBranchLengthParameter object and
    // edge.
    // Nice it explicitly ties between objects -
    // if multiple HiddenState objects share a Tree object,
    // we can reflect that in the map.
    //
    // Similar to t1 parameter shared across two network edges in
    // OptimizeContinuousNetworkModelParameter,
    // except that we can also share across multiple gene trees.

    protected runHmm runHmmObject; // for pushing updated transition probabilities
    protected ParentalTreesDecoration parentalTreesDecoration; // works with CalculationCache directly
    protected CalculationCache calculationCache;

    protected double minimumValue;
    protected double maximumValue;

    public ParentalBranchLengthParameter (String inName, 
					  double inValue, 
					  runHmm inRunHmmObject,
					  ParentalTreesDecoration inParentalTreesDecoration,
					  CalculationCache inCalculationCache,
					  boolean checkValueMinimumConstraintFlag,
					  boolean checkValueMaximumConstraintFlag,
					  boolean updateModelStateFlag) {
	// only sets name
	super(inName);
	
	this.runHmmObject = inRunHmmObject;
	this.parentalTreesDecoration = inParentalTreesDecoration;
	this.calculationCache = inCalculationCache;

	setMinimumValue(DEFAULT_MINIMUM_BRANCH_LENGTH, false);
	setMaximumValue(DEFAULT_MAXIMUM_BRANCH_LENGTH, false);
	
	// need to delay inequality constraint stuff until after construction
	//setValue(inValue, checkValueMinimumConstraintFlag, checkValueMaximumConstraintFlag, false);
	this.value = inValue;

	if (updateModelStateFlag) {
	    updateModelState();
	}
    }

    public void updateModelState () {
	// absolutely no Parameter.setValue() stuff happening in here
	// only modify PhyloNet-HMM state external to Parameter objects

	// propagate length-parameter values (under length-parameter-constraint-set constraints, if applicable)
	// to parental tree branch lengths
	//
	// also clear out caches based on parental branch changes
	// ...
	parentalTreesDecoration.updateBranchesBasedOnParentalBranchLengthParameters(this);

	// ...
	// and finally propagate parental tree branch lengths on to HMM transition probability matrix
	//
	// this will also re-fill and use the caches appropriately
	runHmmObject.updateTransitionProbabilities();
    }

    public double getMinimumValue () {
	return (minimumValue);
    }

    public double getDefaultInitialValue () {
	return (DEFAULT_INITIAL_BRANCH_LENGTH);
    }

    public double getMaximumValue () {
	return (maximumValue);
    }

    // need to call superclass method appropriately
    public void setValue (double inValue, 
			  boolean checkValueMinimumConstraintFlag,
			  boolean checkValueMaximumConstraintFlag,
			  boolean updateModelStateFlag) {
	super.setValue(inValue, checkValueMinimumConstraintFlag, checkValueMaximumConstraintFlag, updateModelStateFlag);

	// setValue on this ParentalBranchLengthParameter object 
	// will affect at most one other ParentalBranchLengthParameter object	
	// according to constraints in ParentalTreeDecoration parsing
	if (parentalTreesDecoration.getLpInequalitiesMap().containsKey(this)) {
	    // this parameter is the lesser parameter
	    // greater parameter must be lower bounded by this object's value
	    parentalTreesDecoration.getLpInequalitiesMap().get(this).setMinimumValue(this.getValue() + EPSILON, checkValueMinimumConstraintFlag);
	}
	else if (parentalTreesDecoration.getLpInequalitiesMap().containsValue(this)) {
	    // this parameter is the greater parameter
	    // lesser parameter must be upper bounded by this object's value
	    parentalTreesDecoration.getLpInequalitiesMap().rget(this).setMaximumValue(this.getValue() - EPSILON, checkValueMaximumConstraintFlag);
	}
    }

    // also add ability to change min/max values
    // needed for dependency constraints
    // 
    // WARNING - checkConstraintsFlag is strict!
    public void setMinimumValue (double inMinimumValue, boolean checkConstraintsFlag) {
	this.minimumValue = inMinimumValue;
	if (checkConstraintsFlag && (getValue() < getMinimumValue())) {
	    throw (new RuntimeException("ERROR: Parameter.setMinimumValue(...) called with value larger than getValue() amount. " + inMinimumValue + " " + getValue()));
	}
    }

    public void setMaximumValue (double inMaximumValue, boolean checkConstraintsFlag) {
	this.maximumValue = inMaximumValue;
	if (checkConstraintsFlag && (getValue() > getMaximumValue())) {
	    throw (new RuntimeException("ERROR: Parameter.setMaximumValue(...) called with value smaller than getValue() amount. " + inMaximumValue + " " + getValue()));
	}
    }
}




    // kliu - crap - need to also constrain max based on rest of HMM state
    // no closed form solution for this???

    // Arghhhh...
    // Four maps needed to properly update model state for a length-parameter, depending
    // on if it belongs to a length-parameter-constraint-set.
    // Annoying detail of layers of indirection between length-parameters and parental tree branches.
    // protected BijectiveHashtable<NetNode<CoalescePattern[]>,String> parentalNodeLabelMap;
    // protected BidirectionalMultimap<ParentalBranchLengthParameter,Tuple<String,String>> lpEidMap;
    // protected BidirectionalMultimap<ParameterConstraintSet,ParentalBranchLengthParameter> setLpMap;
    // protected Hashtable<ParameterConstraintSet,ParentalBranchLengthParameter> setFirstLpMap;

	// this.parentalNodeLabelMap = inParentalNodeLabelMap;
	// this.lpEidMap = inLpEidMap;
	// this.setLpMap = inSetLpMap;
	// this.setFirstLpMap = inSetFirstLpMap;
