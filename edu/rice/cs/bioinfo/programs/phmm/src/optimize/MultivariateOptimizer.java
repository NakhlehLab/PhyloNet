/**
 * Performs expectation-maximization (E-M) of
 * continuous model parameters in PhyloNet-HMM.
 * E-M iterations use a multivariate optimization heuristic that
 * incorporates Brent's method for univariate optimization.
 * Use this approach in lieu of Baum-Welch.
 */

// TODO:
// 1. Cache/restore HMM parameter values.

// 2. Organize HMM parameters, enter them into a queue. 
//    Single iteration == exhaust queue once.

// 3. Add in llh evaluator.

// 4. Alternative initializations: process queue
//    and give either initial, default, or random choices.

// 5. Pull in rest: search interval calculator,
//    brent optimization calls. For both distances and probabilities.

package optimize;

import java.util.List;
import java.util.Set;
import java.util.Map;
import phylogeny.EvoTree;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
import be.ac.ulg.montefiore.run.jahmm.phmm.OpdfMap;
import be.ac.ulg.montefiore.run.jahmm.phmm.TransitionProbabilityParameters;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;
import org.apache.commons.math3.optimization.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optimization.GoalType;

public class MultivariateOptimizer {
    public static final double RELATIVE_ACCURACY = 1e-12;
    public static final double ABSOLUTE_ACCURACY = 1e-8;
    public static final double SEARCH_INTERVAL_MINIMUM_WIDTH = 1e-4;
    public static final int BRENT_METHOD_SINGLE_ROUND_MAXIMUM_ITERATIONS = 100;
    public static final int MAXIMUM_NUM_ROUNDS = 100;
    public static final double MINIMUM_LOG_LIKELIHOOD_DELTA_FOR_CONVERGENCE = ABSOLUTE_ACCURACY;
    public static final double MINIMUM_FACTOR_WIDTH_FOR_MINIMUM_MAXIMUM_INTERVAL = 4.0;

    // search defaults
    // don't allow zero branch lengths
    // 1e-10 too small
    public static final double DEFAULT_MINIMUM_BRANCH_LENGTH = 1e-3;
    public static final double DEFAULT_INITIAL_BRANCH_LENGTH = 1e-1;
    public static final double DEFAULT_MAXIMUM_BRANCH_LENGTH = 1e1; // scientific notation - e for decimal exponent

    // default weights for members of branch-length-constraint-sets
    public static final double DEFAULT_INITIAL_BRANCH_LENGTH_CONSTRAINT_SET_MEMBER_WEIGHT = 1.0;

    // only for randomization purposes
    public static final double MAXIMUM_BRANCH_LENGTH_FOR_RANDOMIZATION_PURPOSES = 2.0;

    // hmm... is GeneTreeProbability able to handle probability == 0 or 1?
    // yes, it handles this fine
    public static final double DEFAULT_MINIMUM_PROBABILITY = 0.0;
    // initial probability always initialized to uniform
    public static final double DEFAULT_MAXIMUM_PROBABILITY = 1.0;

    // hack - to support set branch length constraints
    // PhyloNet's Rich Newick support has some craziness about '_' underscore letter.
    // use dash '-' instead.
    public static final String IDENTIFIER_SET_BRANCH_LENGTH_CONSTRAINTS = "-CONSTRAINT-SET-";

    // to support fixed network probabilities
    // over-parameterization issue with optimization of model of Yu et al. 2012
    //public static final boolean ENABLE_BRANCH_PROBABILITY_OPTIMIZATION_FLAG = true;

    public enum InitialSearchSettings { CURRENT, RANDOM, DEFAULT }





    protected Hmm<ObservationMap> hmm;
    // should really just wrap all into custom HMM class
    protected List<HiddenState> hiddenStates;
    protected TransitionProbabilityParameters transitionProbabilityParameters;
    protected Map<EvoTree,Set<HiddenState>> parentalTreeClasses;

    public MultivariateOptimizer (Hmm<ObservationMap> inHmm,
				  List<HiddenState> inHiddenStates,
				  TransitionProbabilityParameters inTransitionProbabilityParameters,
				  Map<EvoTree,Set<HiddenState>> inParentalTreeClasses
				  ) {
	this.hmm = inHmm;
	this.hiddenStates = inHiddenStates;
	this.transitionProbabilityParameters = inTransitionProbabilityParameters;
	this.parentalTreeClasses = inParentalTreeClasses;
    }



    /**
     * Update associated model values associated with a single parameter
     * e.g., for LengthParameter objects that belong to a LengthParameterConstraintSet
     * or multiple branches share a LengthParameter
     * or multiple parental trees have branches that share a LengthParameter
     * etc.
     * 
     * don't push this into Parameter
     * too complicated
     */
    protected void updateHMM (Parameter parameter) {
	System.err.println ("TODO");
	System.exit(1);
    }

    protected double computeHMMLikelihood () {
	System.err.println ("TODO");
	return (-1.0);
    }


    /**
     * Other than min/max values, no difference between LengthParameter and FrequencyParameter.
     */
    protected class ParameterUnivariateFunction implements MutableUnivariateFunction<Parameter> {
	protected Parameter parameter;

	public ParameterUnivariateFunction
	    (Parameter inParameter) {
	    setParameter(inParameter);
	}

	public void setParameter (Parameter inParameter) {
	    this.parameter = inParameter;
	}

	public Parameter getParameter () {
	    return (parameter);
	}

	/**
	 * Evaluate likelihood function f(x) for a particular x'.
	 * WARNING: f(x) is always a log likelihood, since computeGTProb() always returns a log likelihood!
	 *
	 * Meaning of x differs for nodes that are subject to set branch length constraints:
	 * for these nodes, x is weighted ratio.
	 */
	public double value (double x) {
	    // cache original setting, set new setting
	    double originalSetting = parameter.getValue();
	   
	    // update
	    parameter.setValue(x);
	    // update associated model values associated with a single parameter
	    // e.g., for LengthParameter objects that belong to a LengthParameterConstraintSet
	    // or multiple branches share a LengthParameter
	    // or multiple parental trees have branches that share a LengthParameter
	    // etc.
	    // 
	    // don't push this into Parameter
	    // too complicated
	    updateHMM(parameter);

	    // evaluate f(x) using forward/backwards algorithm
	    double result = computeHMMLikelihood();

	    // restore original setting
	    parameter.setValue(originalSetting);
	    updateHMM(parameter);

	    // above restore original setting op
	    // not really necessary if always accept new branch length even if no likelihood improvement
	    // since Brent's method is guaranteed to never make the likelihood worse
	    // meh, no update if no likelihood improvement
	    // -> above is necessary

	    return (result);
	}	
    }
    
}
