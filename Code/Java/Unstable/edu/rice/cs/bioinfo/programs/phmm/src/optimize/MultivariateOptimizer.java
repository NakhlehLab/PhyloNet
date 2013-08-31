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

// Argh - need to optimize branch lengths for both parental tree and
// gene genealogy. But parental tree is Network<Double> object and
// gene genealogy is Tree<Double> object.
// Just care about getParentDistance() and setParentDistance(...) methods.
// Push into LengthParameter object?
// How to specify in config file?

// Hmm... don't permit any sort of constraints on gene genealogies.
// Simpler that way.

// Only allow parameter sharing and parameter constraints on
// parental trees.
// Force different internal node names on parental trees
// 
//
// Easiest to do it this way.
// Later, if need to have more careful parameterization,
// can change it.

package optimize;

import java.util.List;
import java.util.Set;
import java.util.Map;
import java.util.Hashtable;
import java.util.Vector;
import java.util.StringTokenizer;
import java.io.StringWriter;
import java.io.File;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.IOException;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
import be.ac.ulg.montefiore.run.jahmm.phmm.OpdfMap;
import be.ac.ulg.montefiore.run.jahmm.phmm.TransitionProbabilityParameters;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;
import org.apache.commons.math3.optimization.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optimization.GoalType;
import runHmm.runHmm;
import util.Constants;
import substitutionModel.GTRSubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.library.programming.BijectiveHashtable;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;

public class MultivariateOptimizer {
    public static final double RELATIVE_ACCURACY = 1e-12;
    public static final double ABSOLUTE_ACCURACY = 1e-8;
    public static final double SEARCH_INTERVAL_MINIMUM_WIDTH = 1e-4;
    public static final int BRENT_METHOD_SINGLE_ROUND_MAXIMUM_ITERATIONS = 100;
    public static final int MAXIMUM_NUM_ROUNDS = 10000;
    public static final double MINIMUM_LOG_LIKELIHOOD_DELTA_FOR_CONVERGENCE = 1e-1;
    public static final double MINIMUM_FACTOR_WIDTH_FOR_MINIMUM_MAXIMUM_INTERVAL = 4.0;

    // search defaults
    // don't allow zero branch lengths
    // 1e-10 too small
    public static final double DEFAULT_MINIMUM_BRANCH_LENGTH = 1e-3;
    public static final double DEFAULT_INITIAL_BRANCH_LENGTH = 1e-1;
    public static final double DEFAULT_MAXIMUM_BRANCH_LENGTH = 1e1; // scientific notation - e for decimal exponent

    public static final double DEFAULT_MINIMUM_RATE = 1e-2;
    public static final double DEFAULT_INITIAL_RATE = 1.0;
    public static final double DEFAULT_MAXIMUM_RATE = 1e1; // scientific notation - e for decimal exponent

    // only for randomization purposes
    //
    // only go to 80% of max for randomization purposes
    public static final double FRACTION_OF_MAXIMUM_FOR_RANDOMIZATION_PURPOSES = 0.8;

    // only for initializing ParentalBranchLengthParameter objects that are subject to 
    // inequality constraints
    public static final double INITIALIZATION_RELATIVE_WEIGHT_FOR_INEQUALITY_CONSTRAINTS = 2.0;

    // hmm... is GeneTreeProbability able to handle probability == 0 or 1?
    // yes, it handles this fine
    public static final double DEFAULT_MINIMUM_PROBABILITY = 0.0;
    public static final double DEFAULT_INITIAL_PROBABILITY = 1e-2;
    // initial probability always initialized to uniform
    // don't let switching probability get to close to half
    public static final double DEFAULT_MAXIMUM_PROBABILITY = 0.25;
    // public static final double DEFAULT_MAXIMUM_PROBABILITY = 1.0;

    // hack - to support set branch length constraints
    // PhyloNet's Rich Newick support has some craziness about '_' underscore letter.
    // use dash '-' instead.
    public static final String IDENTIFIER_SET_BRANCH_LENGTH_CONSTRAINTS = "-CONSTRAINT-SET-";

    public static final String CHECKPOINT_PARAMETER_STATE_FILENAME = "CHECKPOINT-PARAMETER-STATE";
    public static final String CHECKPOINT_LIKELIHOOD_FILENAME = "CHECKPOINT-LIKELIHOOD";
    public static final String CHECKPOINT_VITERBI_SEQUENCE_FILENAME = "CHECKPOINT-VITERBI-SEQUENCE";
    public static final String CHECKPOINT_ENTRY_PAIR_DELIMITER_CHARACTER = "\t";
    //public static final String CHECKPOINT_LATEST_SUFFIX = "LATEST";

    public static final String FILENAME_SUFFIX_DELIMITER = ".";

    // to support fixed network probabilities
    // over-parameterization issue with optimization of model of Yu et al. 2012
    //public static final boolean ENABLE_BRANCH_PROBABILITY_OPTIMIZATION_FLAG = true;

    public enum InitialSearchSettings { CURRENT, RANDOM, DEFAULT }




    // carry over references to runHmm/etc. state
    protected Hmm<ObservationMap> hmm;
    // need a few helper routines in this object
    // lots of shared state between this class and runHmm class
    // that's fine
    protected runHmm runHmmObject;
    // should really just wrap all into custom HMM class
    protected List<HiddenState> hiddenStates;
    protected TransitionProbabilityParameters transitionProbabilityParameters;
    // bleh - need to assume a GTR model here due to specific parameterization
    protected GTRSubstitutionModel gtrSubstitutionModel;
    //protected Map<Network<Double>,Set<HiddenState>> parentalTreeClasses;
    protected BijectiveHashtable<String,Network<Double>> parentalTreeNameMap;
    protected BijectiveHashtable<String,Tree> geneGenealogyNameMap;
    protected List<ObservationMap> observation;
    protected CalculationCache calculationCache;

    // create parental tree decorations/parameter maps
    protected ParentalTreesDecoration parentalTreesDecoration;
    
    // kliu - implement this later
    protected boolean enableParentalTreeOptimizationFlag = true;
    protected boolean enableGeneGenealogyOptimizationFlag = true;
    protected boolean enableSwitchingFrequencyOptimizationFlag = true;
    protected boolean enableSubstitutionModelOptimizationFlag = true;

    // heuristic for univariate optimization
    protected BrentOptimizer brentOptimizer;
    
    protected Vector<ParentalBranchLengthParameter> parentalBranchLengthParameters;
    protected Vector<GenealogyBranchLengthParameter> genealogyBranchLengthParameters;
    protected Vector<SwitchingFrequencyParameter> switchingFrequencyParameters;
    protected Vector<GTRRateParameter> gtrRateParameters;
    protected Vector<GTRBaseFrequencyParameter> gtrBaseFrequencyParameters;

    // for convenience, keep list of Parameters
    // and bijective map between parameters and their names
    //
    // force parameters to have unique names and no tabs, otherwise barf
    protected BijectiveHashtable<String,Parameter> parameterNameMap;

    // meh - just keep above for clarity
    // use lpEidMap.keys() to get a list of all parental tree branch length parameters (ParentalBranchLengthParameter objects)

    public MultivariateOptimizer (Hmm<ObservationMap> inHmm,
				  // HMM update routines are located in this object
				  runHmm inRunHmm,
				  List<HiddenState> inHiddenStates,
				  TransitionProbabilityParameters inTransitionProbabilityParameters,
				  GTRSubstitutionModel inGTRSubstitutionModel,
				  //Map<Network<Double>,Set<HiddenState>> inParentalTreeClasses,
				  BijectiveHashtable<String,Network<Double>> inParentalTreeNameMap,
				  BijectiveHashtable<String,Tree> inGeneGenealogyNameMap,
				  List<ObservationMap> inObservation,
				  String inputParentalBranchLengthParameterToEdgeMapFilename,
				  String inputParentalBranchLengthParameterInequalitiesFilename,
				  String inputParentalBranchLengthParameterSetConstraintsFilename,
				  CalculationCache inCalculationCache,
				  boolean inEnableParentalTreeOptimizationFlag,
				  boolean inEnableGeneGenealogyOptimizationFlag,
				  boolean inEnableSwitchingFrequencyOptimizationFlag,
				  boolean inEnableSubstitutionModelOptimizationFlag,
				  String inputRestoreCheckpointFilename
				  ) {
	this.hmm = inHmm;
	this.runHmmObject = inRunHmm;
	this.hiddenStates = inHiddenStates;
	this.transitionProbabilityParameters = inTransitionProbabilityParameters;
	this.gtrSubstitutionModel = inGTRSubstitutionModel;
	//this.parentalTreeClasses = inParentalTreeClasses;
	this.parentalTreeNameMap = inParentalTreeNameMap;
	this.geneGenealogyNameMap = inGeneGenealogyNameMap;
	this.observation = inObservation;
	this.calculationCache = inCalculationCache;
	this.enableParentalTreeOptimizationFlag = inEnableParentalTreeOptimizationFlag;
	this.enableGeneGenealogyOptimizationFlag = inEnableGeneGenealogyOptimizationFlag;
	this.enableSwitchingFrequencyOptimizationFlag = inEnableSwitchingFrequencyOptimizationFlag;
	this.enableSubstitutionModelOptimizationFlag = inEnableSubstitutionModelOptimizationFlag;

	brentOptimizer = new BrentOptimizer(RELATIVE_ACCURACY, ABSOLUTE_ACCURACY);

	verifySearchSettings();

	// ParentalBranchLengthParameter objects created in here
	parentalTreesDecoration = new ParentalTreesDecoration(parentalTreeNameMap,
							      inputParentalBranchLengthParameterToEdgeMapFilename,
							      inputParentalBranchLengthParameterInequalitiesFilename,
							      inputParentalBranchLengthParameterSetConstraintsFilename,
							      runHmmObject,
							      calculationCache);
	parentalBranchLengthParameters = parentalTreesDecoration.getParentalBranchLengthParameters();

	// no model update
	createGenealogyBranchLengthParameters();
	createSwitchingFrequencyParameters();
	createSubstitutionModelParameters();	

	// for convenience only
	createNameMapOfAllParameters();

	// do restore if necessary
	if ((inputRestoreCheckpointFilename != null) &&
	    (!inputRestoreCheckpointFilename.trim().equals(""))) {
	    if (Constants.WARNLEVEL > 1) { System.out.println ("Restoring from checkpoint file " + inputRestoreCheckpointFilename + "."); }
	    restoreParameterValuesFromCheckpointWithNoUpdate(inputRestoreCheckpointFilename);
	    if (Constants.WARNLEVEL > 1) { System.out.println ("Restoring from checkpoint file " + inputRestoreCheckpointFilename + " DONE."); }
	}

	// finally do update
	updateModelStateForEveryParameter();
    }

    protected void createSubstitutionModelParameters () {
	gtrBaseFrequencyParameters = new Vector<GTRBaseFrequencyParameter>();
	double[] currentStationaryProbabilities = gtrSubstitutionModel.getStationaryProbabilities();
	for (int i = 0; i < gtrSubstitutionModel.getAlphabet().length(); i++) {
	    gtrBaseFrequencyParameters.add(new GTRBaseFrequencyParameter(GTRBaseFrequencyParameter.class.getName() + HiddenState.HIDDEN_STATE_NAME_DELIMITER + Character.toString(gtrSubstitutionModel.getAlphabet().getAlphabet().charAt(i)),
									 // by convention, first parameter in a constraint-set
									 // gets set to 1.0
									 //
									 // need to normalizes into relative weights
									 // \pi_A weight canonically set to 1.0
									 (i == 0) ? 1.0 : currentStationaryProbabilities[i] / currentStationaryProbabilities[0], // is this syntax ok?
									 gtrSubstitutionModel,
									 i,
									 gtrBaseFrequencyParameters,
									 calculationCache,
									 true,
									 true,
									 false // no need to update
									 ));
	}

	gtrRateParameters = new Vector<GTRRateParameter>();
	double[] currentSubstitutionRateParameters = gtrSubstitutionModel.getOriginalRateParameters();
	for (int i = 0; i < gtrSubstitutionModel.getRateParameterCount(); i++) {
	    gtrRateParameters.add(new GTRRateParameter(GTRRateParameter.class.getName() + HiddenState.HIDDEN_STATE_NAME_DELIMITER + Integer.toString(i),
						       currentSubstitutionRateParameters[i], // is this syntax ok?
						       gtrSubstitutionModel,
						       i,
						       calculationCache,
						       true,
						       true,
						       false // no need to update
						       ));
	}
    }

    protected void createSwitchingFrequencyParameters () {
	switchingFrequencyParameters = new Vector<SwitchingFrequencyParameter>();
	for (TransitionProbabilityParameters.ParameterChoice parameterChoice : TransitionProbabilityParameters.ParameterChoice.values()) {
	    SwitchingFrequencyParameter sfp = new SwitchingFrequencyParameter(SwitchingFrequencyParameter.class.getName() + HiddenState.HIDDEN_STATE_NAME_DELIMITER + parameterChoice.toString(),
									      transitionProbabilityParameters.get(parameterChoice),
									      runHmmObject,
									      transitionProbabilityParameters,
									      parameterChoice,
									      true,
									      true,
									      false // no need to update
									      );
	    //runHmmObject,
	    switchingFrequencyParameters.add(sfp);
	}
    }

    /**
     * In current model, hidden states share both parental trees and
     * gene genealogies. Gene genealogy parameterization done
     * per object only.
     */
    protected void createGenealogyBranchLengthParameters () {
	genealogyBranchLengthParameters = new Vector<GenealogyBranchLengthParameter>();
	for (Tree geneGenealogy : geneGenealogyNameMap.values()) {
	    for (TNode node : geneGenealogy.postTraverse()) {
		if (node.isRoot() || node.getParent() == null) {
		    continue;
		}
		
		GenealogyBranchLengthParameter gblp = new GenealogyBranchLengthParameter(GenealogyBranchLengthParameter.class.getName() + HiddenState.HIDDEN_STATE_NAME_DELIMITER + geneGenealogyNameMap.rget(geneGenealogy) + HiddenState.HIDDEN_STATE_NAME_DELIMITER + node.getName(),
											 node.getParentDistance(),
											 node,
											 calculationCache,
											 true,
											 true,
											 // no need for update
											 false);
		genealogyBranchLengthParameters.add(gblp);
	    }
	}
    }

    protected void createNameMapOfAllParameters () {
	Vector<Parameter> parameters = new Vector<Parameter>();
	parameters.addAll(parentalBranchLengthParameters); // parental tree branch length parameters
	parameters.addAll(genealogyBranchLengthParameters);
	parameters.addAll(switchingFrequencyParameters);
	parameters.addAll(gtrRateParameters);
	parameters.addAll(gtrBaseFrequencyParameters);

	parameterNameMap = new BijectiveHashtable<String,Parameter>();
	for (Parameter parameter : parameters) {
	    // strict!
	    if (parameterNameMap.containsKey(parameter.getName())) {
		// barf
		throw (new RuntimeException("ERROR: parameters must have unique names. Check for duplicate parameter name " + parameter.getName() + "."));
	    }

	    if (parameter.getName().contains(CHECKPOINT_ENTRY_PAIR_DELIMITER_CHARACTER)) {
		throw (new RuntimeException("ERROR: parameter names cannot contain the character used to delimit checkpoint pair items. Check parameter name " + parameter.getName() + "."));
	    }

	    parameterNameMap.put(parameter.getName(), parameter);
	}
    }

    protected boolean verifyProbability (double p) {
	if ((p >= 0.0) && (p <= 1.0)) {
	    return (true);
	}
	else {
	    return (false);
	}
    }

    /**
     * Paranoid!
     * Warning - this method throws unchecked exceptions.
     * Make sure that search settings make sense.
     */
    protected void verifySearchSettings () {
	// disallow zero probabilities
	// doesn't make sense anyways - don't need to collapse model parameters this way
	// protect against degenerate likelihood calculations/parameter optimizations
	if (!verifyProbability(DEFAULT_MINIMUM_PROBABILITY) || 
	    !verifyProbability(DEFAULT_MAXIMUM_PROBABILITY)) {
	    throw(new IllegalSearchIntervalException("ERROR: probability minimum/maximum settings are not proper probabilities."));
	}

	if ((DEFAULT_MINIMUM_PROBABILITY > DEFAULT_MAXIMUM_PROBABILITY)
	    ) {
	    throw(new IllegalSearchIntervalException("ERROR: probability minimum/maximum settings are out of order."));
	}

	// ditto with branch lengths
	if ((DEFAULT_MINIMUM_BRANCH_LENGTH < 0.0) ||
	    (DEFAULT_INITIAL_BRANCH_LENGTH < 0.0) ||
	    (DEFAULT_MAXIMUM_BRANCH_LENGTH < 0.0)) {
	    throw(new IllegalSearchIntervalException("ERROR: negative branch length search settings are not allowed."));
	}

	if ((DEFAULT_MINIMUM_BRANCH_LENGTH > DEFAULT_INITIAL_BRANCH_LENGTH) ||
	    (DEFAULT_INITIAL_BRANCH_LENGTH > DEFAULT_MAXIMUM_BRANCH_LENGTH)) {
	    throw(new IllegalSearchIntervalException("ERROR: branch length minimum/initial/maximum settings are out of order."));
	}

	if (DEFAULT_MINIMUM_BRANCH_LENGTH * MINIMUM_FACTOR_WIDTH_FOR_MINIMUM_MAXIMUM_INTERVAL > DEFAULT_MAXIMUM_BRANCH_LENGTH) {
	    throw(new IllegalSearchIntervalException("ERROR: branch length minimum/maximum bounds width is too small."));
	}
    }

    /**
     * Convenience function. Use this if you called Parameter.setValue(...) for every Parameter object
     * with updateModelStateFlag set to false.
     */
    protected void updateModelStateForEveryParameter () {
	// duplicated work here for each member of a length-parameter-constraint-set
	// oh well
	// ok since each call below just rebalances weight to satisfy constraint on a constraint-set
	// repeated re-balancings don't change anything, so long as individual
	// parameter values are unchanged
	for (Parameter p : parameterNameMap.values()) {
	    // manually update Jahmm's probability matrices/vectors *once* at the end
	    p.updateModelState();
	}
    }

	// see above comment
    //runHmmObject.updateTransitionProbabilities();


    protected boolean checkConvergence (int round, double prevLogLikelihood, double currLogLikelihood) {
	if (round >= MAXIMUM_NUM_ROUNDS) {
	    if (Constants.WARNLEVEL > 4) { System.out.println ("Converged in checkConvergence(): Maximum number of optimization rounds reached."); }
	    return (true);
	}

	// at least one round required before likelihood comparisons can occur
	if (round > 0) {
	    if (Constants.WARNLEVEL > 1) { System.out.println ("Checking convergence for previous log likelihood " + prevLogLikelihood + " and current log likelihood " + currLogLikelihood); }

	    // paranoid!
	    if (prevLogLikelihood > currLogLikelihood) {
		if (Constants.WARNLEVEL > 1) { System.out.println ("ERROR: current round's likelihood is worse than previous round's likelihood! Check optimize() guarantees. Proceeding anyways."); }
	    }

	    double absdiff = Math.abs(currLogLikelihood - prevLogLikelihood);
	    if (absdiff < MINIMUM_LOG_LIKELIHOOD_DELTA_FOR_CONVERGENCE) {
		if (Constants.WARNLEVEL > 1) { System.out.println ("Converged in checkConvergence(): log likelihood delta smaller than cutoff. " + absdiff + " " + MINIMUM_LOG_LIKELIHOOD_DELTA_FOR_CONVERGENCE); }
		return (true);
	    }
	}

	return (false);
    }

    /**
     * Need to satisfy four constraints:
     * 0. min, max make sense - see verifySearchSettings()
     * 1. x \in [l, u]
     * 2. x \in [min, max] and similarly for l and u
     * 3. f(x) is better than both f(l) and f(u)
     *
     * No log use here. Consistent with GeneTreeProbability library.
     *
     * Output is (l', x', u') search interval that satisfies above four constraints.
     *
     * Warning - these method only works for nonnegative x, l, u, min, max.
     */
    protected Tuple3<Double,Double,Double> getSearchInterval (UnivariateFunction f,
							      //double l, // not used
							      double x, 
							      //double u, // not used
							      double min, 
							      double max, 
							      String debugMessage // for debugging purposes
							      ) {
	// To proceed with search-by-halves,
	// max must be at least 4X the magnitude of min.
	// If this constraint fails, just do a simple search bounded by min/max, retaining current x.
	if (max < MINIMUM_FACTOR_WIDTH_FOR_MINIMUM_MAXIMUM_INTERVAL * min) {
	    // paranoid
	    if (min > max) {
		throw (new RuntimeException("ERROR: invalid input min/max in getSearchInterval(...). " + min + " " + max));
	    }

	    if ((x < min) || (x > max)) {
		throw (new RuntimeException("ERROR: x out of min/max range in getSearchInterval(...). " + min + " " + x + " " + max));
	    }

	    return (new Tuple3<Double,Double,Double>(min, x, max));
	}

	// tweak these and x -> output
	double l;
	double u;

	if (Constants.WARNLEVEL > 4) { System.out.println (debugMessage + " search interval min and max: |" + min + " " + max + "|"); }

	// constraint #0 guaranteed in verifySearchSettings()
	// in case some set() calls have happened after construction

	// constraint #1 satisfied throughout this function
	
	// constraint #2:
	// check first steps
	double lFirst = x / 2.0;
	double uFirst = x * 2.0;
	if (min > lFirst) { // first step of search-by-halves results in search interval below lower bound
	    // search at interval adjoining min
	    l = min;
	    x = min * 2.0;
	    u = min * 4.0;
	}
	else if (uFirst > max) {
	    // search at interval adjoining max
	    u = max;
	    x = max / 2.0;
	    l = max / 4.0;
	}
	else {
	    // default search halves around probe x
	    l = lFirst;
	    u = uFirst;
	}

	// looks reasonable
	//if (Constants.WARNLEVEL > 4) { System.out.println (debugMessage + " search setting after satisfying constraint #2 in getSearchInterval(): |" + l + " " + x + " " + u + "|"); }

	// testing
	if (Constants.WARNLEVEL > 4) { System.out.println (debugMessage + " l x u fl fx fu log likelihoods: |" + l + " " + x + " " + u + " " + f.value(l) + " " + f.value(x) + " " + f.value(u) + "|"); }

	// This doesn't really work that well for probabilities.
	// With x = 0.5 and maximum probability of 1, doubling doesn't work well.
	// Also, clunky to keep halving until loss of precision -> 0.
	// Smallest positive double = 2^-1074 = 8^358 ~ 10^324 or so.
	// Shouldn't rely on underflow to get to minimum.

	// TODO:
	// - If difference between l and min is smaller than ABSOLUTE_ACCURACY, then stop and set l = min.
	// - Instead of always doubling u towards max, if doubling u -> greater than max, then 
	//   advance half the distance between u and max (similarly to always halving distance between l and min).
	//   Need to also watch out for above underflow issues.

	// Wait a sec - why is this only getting called on branch with probability less than half??
	
	// Search might also be stuck in a local optimum. Try alternative initial search settings.
	// Try this before code changes above.

	// constraint #3:
	// expand search interval endpoints by half/double until satisfied
	while (!(f.value(x) > f.value(l))) {
	    // don't just walk over the edge
	    // if next move takes us over the edge
	    // then just approach the edge by halves
	    // until "close enough"

	    // search width is big enough
	    if (Math.abs(l - min) > SEARCH_INTERVAL_MINIMUM_WIDTH) {
		if (l / 2.0 > min) {
		    // next step won't take us over the edge
		    l /= 2.0;
		}
		else {
		    // otherwise, approach the edge by halves
		    l -= (Math.abs(l - min) / 2.0);
		}
	    }
	    else {
		// go to the edge
		l = min;
		break;
	    }
	    
	    // if l == 0.0 exactly,
	    // then divide-by-two op does nothing
	    // need to set up min to guard against this case

	    if (Constants.WARNLEVEL > 4) { System.out.println (debugMessage + " l x u fx fl fu log likelihoods searching l: |" + l + " " + x + " " + u + " " + f.value(l) + " " + f.value(x) + " " + f.value(u) + "|"); }
	}

	while (!(f.value(x) > f.value(u))) {
	    // process u same as l

	    // search width is big enough
	    if (Math.abs(max - u) > SEARCH_INTERVAL_MINIMUM_WIDTH) {
		if (u * 2.0 < max) {
		    // next step won't take us over the edge
		    u *= 2.0;
		}
		else {
		    // approach the edge by halves
		    u += (Math.abs(max - u) / 2.0);
		}
	    }
	    else {
		// go to the edge
		u = max;
		break;
	    }

	    // if u == 0.0 exactly,
	    // then multiply-by-two op does nothing
	    // need to set up max to guard against this case

	    if (Constants.WARNLEVEL > 4) { System.out.println (debugMessage + " l x u fx fl fu log likelihoods searching u: |" + l + " " + x + " " + u + " " + f.value(l) + " " + f.value(x) + " " + f.value(u) + "|"); }
	}

	if (Constants.WARNLEVEL > 4) { System.out.println (debugMessage + " l x u fx fl fu log likelihoods after satisfying constraint #3 in getSearchInterval() : |" + l + " " + x + " " + u + " " + f.value(l) + " " + f.value(x) + " " + f.value(u) + "|"); }

	if (Constants.WARNLEVEL > 2) { System.out.println (debugMessage + " l x u fx fl fu log likelihoods at end of getSearchInterval() : |" + l + " " + x + " " + u + " " + f.value(l) + " " + f.value(x) + " " + f.value(u) + "|"); }

	Tuple3<Double,Double,Double> result = new Tuple3<Double,Double,Double>
	    (new Double(l), new Double(x), new Double(u));
	return (result);
    }



    /**
     * Helper function.
     * Only update branch length if likelihood is strictly better.
     */
    protected double optimizeSingleParameter (Parameter p, double logLikelihood, int round) {
	ParameterUnivariateFunction f = new ParameterUnivariateFunction(p);
	// l, x, u
	Tuple3<Double,Double,Double> searchInterval = getSearchInterval
	    (f, 
	     p.getValue(), 
	     p.getMinimumValue(), 
	     p.getMaximumValue(),
	     "optimizeSingleParameter: ");
	UnivariatePointValuePair upvp = brentOptimizer.optimize(BRENT_METHOD_SINGLE_ROUND_MAXIMUM_ITERATIONS,
								f,
								GoalType.MAXIMIZE,
								searchInterval.Item1.doubleValue(), // l
								searchInterval.Item3.doubleValue(), // u
							        searchInterval.Item2.doubleValue()); // x

	if (Constants.WARNLEVEL > 2) { System.out.println ("Brent optimization point: |" + upvp.getPoint() + "| likelihood: |" + upvp.getValue() + "|"); }

	double brentOptimizedLogLikelihood = upvp.getValue();
	// see function comments
	if (brentOptimizedLogLikelihood > logLikelihood) {
	    if (Constants.WARNLEVEL > 2) { System.out.println ("INFO: Brent's method resulted in strict improvement in likelihood. Updating."); }

	    // update 
	    p.setValue(upvp.getPoint());
	    //updateHMM(p); // migrated to setValue
	    logLikelihood = brentOptimizedLogLikelihood;
	}
	else {
	    // no update - info instead
	    if (Constants.WARNLEVEL > 2) { System.out.println ("INFO: Round " + round + " optimized point and log likelihood for length parameter " + p.getName() + " resulted in log likelihood " + brentOptimizedLogLikelihood + " which isn't better than current log likelihood " + logLikelihood + ". Not updating branch length nor best round log likelihood."); }
	}

	return (logLikelihood);
    }

    /**
     * For diagnostic purposes.
     */
    // protected void debugParameters () {
    // 	for (Parameter p : parameterNameMap.values()) {
    // 	    System.out.println (p.toString());
    // 	}
    // }

    /**
     * For diagnostic purposes.
     */
    // protected void debugModel () {
    // 	System.out.println ("============================================");
    // 	System.out.println ("Current PhyloNet-HMM state: ");
    // 	System.out.println ();
    // 	System.out.println ("Hidden states: ");
    // 	for (int i = 0; i < hiddenStates.size(); i++) {
    // 	    System.out.println("Hidden state " + i + ":");
    // 	    System.out.println(hiddenStates.get(i).toString());
    // 	}
    // 	System.out.println ();
    // 	System.out.println ("Parameters: ");
    // 	for (Parameter parameter : parameterNameMap.values()) {
    // 	    System.out.println (parameter.toString());
    // 	}
    // 	System.out.println ();
    // 	System.out.println ("HMM transition and emission probabilities: ");
    // 	System.out.println (hmm.toString());
    // 	System.out.println ("============================================");
    // }

    // (int pass, InitialSearchSettings initialSearchSettings)

    /**
     * Perform one optimization pass from a single starting point.
     * Either starting point uses current settings, random settings, or default initial rate/uniform probabilities.
     */
    protected double singlePassOptimization (int pass, InitialSearchSettings initialSearchSettings) {
	double inputLogLikelihood = computeHMMLikelihood();

	if (Constants.WARNLEVEL > 1) { 
	    System.out.println ("Processing pass " + pass + " with initial search setting " + initialSearchSettings.toString() + ".");
	}

	if (Constants.WARNLEVEL > 1) { 
	    //debugModel();
	    System.out.println ("Input log likelihood: |" + inputLogLikelihood + "|");
	}

	// paranoid
	if (Double.isNaN(inputLogLikelihood)) {
	    System.err.println ("ERROR: unable to evaluate input log likelihood. Aborting and returning NaN.");
	    return (Double.NaN);
	}
	
	// Cuong's network representation:
	// Nodes contain name, one datum, and links to children and parents (also "hidden"
	// members with distances, support, and probabilities - unnecessary).
	// Distances always associated with child->parent reference.
	// Man, children and parent lists are LinkedLists. Fine if node degrees are small.
	// Wrapper class just adds in root designation for a particular node
	// and utility functions.
	// Nice and simple.
	switch (initialSearchSettings) {
	case CURRENT:
	    // NOOP
	    break;
	case RANDOM:
	    initializeRandom();
	    break;
	default:
	    initializeDefault();
	    break;
	}

	double initialLogLikelihood = computeHMMLikelihood();

	// paranoid
	if (Double.isNaN(initialLogLikelihood)) {
	    System.err.println ("ERROR: unable to evaluate initial log likelihood. Aborting and returning NaN.");
	    return (Double.NaN);
	}

	// looks good
	//
	// child.inDeg >= 2 -> only then parent probabilities printed
	// defaults to probability zero
	if (Constants.WARNLEVEL > 1) { 
	    //debugModel();
	    System.out.println ("Initial log likelihood: |" + initialLogLikelihood + "|");
	    // // Also print out length-parameters.
	    // System.out.println ("Initial parameters:");
	    // debugParameters();
	}

	// kliu - run iterative heuristic
	// in each iteration, run Brent's optimization method
	// to set each model parameter to local optimum
	// repeat in rounds until convergence
	int round = 0;
	// previous round's log likelihood
	double prevLogLikelihood = initialLogLikelihood; // checkConvergence() doesn't look at this during first round
	double currLogLikelihood = initialLogLikelihood; // checkConvergence() doesn't look at this during first round

	while (!checkConvergence(round, prevLogLikelihood, currLogLikelihood)) {
	    if (Constants.WARNLEVEL > 1) { System.out.println ("Processing round " + round + "."); }

	    double roundLogLikelihood = initialLogLikelihood;

	    // iterate through parameters
	    // may want to make order random??
	    int pCount = 0;
	    for (Parameter p : parameterNameMap.values()) {
		// skip parental tree branch length optimization if disabled
		if (!enableParentalTreeOptimizationFlag && 
		    (p instanceof ParentalBranchLengthParameter)
		    ) {
		    continue;
		}

		// skip gene genealogy branch length optimization if disabled
		if (!enableGeneGenealogyOptimizationFlag &&
		    (p instanceof GenealogyBranchLengthParameter)) {
		    continue;
		}

		// skip frequency optimization if disabled
		if (!enableSwitchingFrequencyOptimizationFlag &&
		    (p instanceof SwitchingFrequencyParameter)) {
		    continue;
		}

		// skip substitution model optimization if disabled
		if (!enableSubstitutionModelOptimizationFlag && 
		    ((p instanceof GTRRateParameter) ||
		     (p instanceof GTRBaseFrequencyParameter))) {
		    continue;
		}

		// by convention, skip the first lengthParameter
		// in a length-parameter-constraint-set
		if ((p instanceof ParentalBranchLengthParameter) && 
		    parentalTreesDecoration.checkFirstInParameterConstraintSet((ParentalBranchLengthParameter) p)) {
		    continue;
		}

		// skip first substitution model first base frequency parameter
		if ((p instanceof GTRBaseFrequencyParameter) && 
		    (gtrBaseFrequencyParameters.get(0) == (GTRBaseFrequencyParameter) p)) {
		    continue;
		}

		if (Constants.WARNLEVEL > 1) { System.out.println ("Processing " + p.getClass().getName() + " parameter " + p.getName() + " count " + pCount + "."); }

		// single branch length optimization
		// does update appropriately
		roundLogLikelihood = optimizeSingleParameter(p, roundLogLikelihood, round);

		// paranoid
		if (Double.isNaN(roundLogLikelihood)) {
		    System.err.println ("ERROR: unable to evaluate round " + round + " parameter " + p.getName() + " log likelihood. Aborting and returning NaN.");
		    return (Double.NaN);
		}

		if (Constants.WARNLEVEL > 1) { System.out.println ("Processing " + p.getClass().getName() + " parameter " + p.getName() + " count " + pCount + " DONE."); }
		pCount++;
	    }

	    // end of round
	    // update likelihoods
	    // if we did the swap at the start of each iteration, wouldn't need the third variable - meh
	    prevLogLikelihood = currLogLikelihood;
	    currLogLikelihood = roundLogLikelihood;

	    // kliu - cache state to disk also
	    writeCheckpoint(pass, round, currLogLikelihood);

	    if (Constants.WARNLEVEL > 1) { System.out.println ("Processing round " + round + " DONE."); }

	    round++;
	}

	// recompute for debugging output
	double resultLikelihood = computeHMMLikelihood();

	if (Constants.WARNLEVEL > 1) { 
	    //debugModel();
	    System.out.println ("Pass log likelihood: |" + resultLikelihood + "|");	    
	    // // Also print out length-parameters.
	    // System.out.println ("Pass parameters:");
	    // debugParameters();
	    System.out.println ("Processing pass " + pass + " with initial search setting " + initialSearchSettings.toString() + " DONE.");
	}

	return (resultLikelihood);
    }

    /**
     * Caller's responsibility to make sure that HMM is up to date prior to calling 
     * this method.
     * Hmm... this is the only place that we actually compute a model likelihood,
     * and it's a log likelihood.
     */
    protected double computeHMMLikelihood () {
	return (hmm.lnProbability(observation));
    }

    protected void initializeDefault () {
	if (enableParentalTreeOptimizationFlag) {
	    for (ParentalBranchLengthParameter parentalLengthParameter : parentalBranchLengthParameters) {
		if (parentalTreesDecoration.checkFirstInParameterConstraintSet(parentalLengthParameter)) {
		    // canonical!
		    parentalLengthParameter.setValue(1.0, true, true, false);
		}
		else {
		    // skip parameters subject to inequality constraints, 
		    // handle this separately below
		    if (parentalTreesDecoration.getLpInequalitiesMap().containsKey(parentalLengthParameter)) {
			parentalLengthParameter.setValue(parentalLengthParameter.getDefaultInitialValue() / INITIALIZATION_RELATIVE_WEIGHT_FOR_INEQUALITY_CONSTRAINTS,
							 false,
							 false,
							 false);
		    } 
		    else if (parentalTreesDecoration.getLpInequalitiesMap().containsValue(parentalLengthParameter)) {
			parentalLengthParameter.setValue(parentalLengthParameter.getDefaultInitialValue(),
							 false,
							 false,
							 false);
		    }
		    else {
			// no inequality constraint
			parentalLengthParameter.setValue(parentalLengthParameter.getDefaultInitialValue(), true, true, false);
		    }
		}
	    }
	}

	if (enableGeneGenealogyOptimizationFlag) {
	    for (GenealogyBranchLengthParameter gblp : genealogyBranchLengthParameters) {
		gblp.setValue(gblp.getDefaultInitialValue(), true, true, false);
	    }
	}	    

	if (enableSwitchingFrequencyOptimizationFlag) {
	    for (SwitchingFrequencyParameter sfp : switchingFrequencyParameters) {
		sfp.setValue(sfp.getDefaultInitialValue(), true, true, false);
	    }
	}

	if (enableSubstitutionModelOptimizationFlag) {
	    for (GTRRateParameter grp : gtrRateParameters) {
		grp.setValue(grp.getDefaultInitialValue(), true, true, false);
	    }
	    
	    for (int i = 0; i < gtrBaseFrequencyParameters.size(); i++) {
		GTRBaseFrequencyParameter gbfp = gtrBaseFrequencyParameters.get(i);
		if (i == 0) {
		    // by convention
		    gbfp.setValue(1.0, true, true, false);
		}
		else {
		    gbfp.setValue(gbfp.getDefaultInitialValue(), true, true, false);
		}
	    }
	}

	updateModelStateForEveryParameter();
    }

    /**
     * Can't let all branch lengths approach max.
     * Can cause model likelihood calculation to barf.
     * Cap the maximum randomized branch length.
     */
    protected void initializeRandom () {
	if (enableParentalTreeOptimizationFlag) {
	    for (ParentalBranchLengthParameter parentalLengthParameter : parentalBranchLengthParameters) {
		if (parentalTreesDecoration.checkFirstInParameterConstraintSet(parentalLengthParameter)) {
		    // canonical!
		    parentalLengthParameter.setValue(1.0, true, true, false);
		}
		else {
		    double parameterMinimumValue = parentalLengthParameter.getMinimumValue();
		    double parameterMaximumValue = FRACTION_OF_MAXIMUM_FOR_RANDOMIZATION_PURPOSES * parentalLengthParameter.getMaximumValue();
		    // uninitialized
		    double randomDistance1 = -1.0;
		    double randomDistance2 = -1.0;
		    do {
			randomDistance1 = Math.random() * (parameterMaximumValue - parameterMinimumValue) + parameterMinimumValue;
			randomDistance2 = Math.random() * (parameterMaximumValue - parameterMinimumValue) + parameterMinimumValue;
		    } while ((randomDistance1 < 0.0) || (randomDistance2 < 0.0) || (randomDistance1 == randomDistance2));
		    
		    if (randomDistance2 < randomDistance1) {
			double tr = randomDistance2;
			randomDistance2 = randomDistance1;
			randomDistance1 = tr;
		    }

		    // meh, just set it twice
		    // parentalLengthParameter is lesser
		    if (parentalTreesDecoration.getLpInequalitiesMap().containsKey(parentalLengthParameter)) {
			parentalLengthParameter.setValue(randomDistance1, false, false, false);
			parentalTreesDecoration.getLpInequalitiesMap().get(parentalLengthParameter).setValue(randomDistance2, false, false, false);
		    }
		    else if (parentalTreesDecoration.getLpInequalitiesMap().containsValue(parentalLengthParameter)) {
			parentalTreesDecoration.getLpInequalitiesMap().rget(parentalLengthParameter).setValue(randomDistance1, false, false, false);
			parentalLengthParameter.setValue(randomDistance2, false, false, false);

		    }
		    else {
			// double parameterMinimumValue = parentalLengthParameter.getMinimumValue();
			// double parameterMaximumValue = FRACTION_OF_MAXIMUM_FOR_RANDOMIZATION_PURPOSES * parentalLengthParameter.getMaximumValue();
			// double randomDistance = Math.random() * (parameterMaximumValue - parameterMinimumValue) + parameterMinimumValue;
			parentalLengthParameter.setValue(randomDistance1, true, true, false);
		    }
		}
	    }
	}

	if (enableGeneGenealogyOptimizationFlag) {
	    for (GenealogyBranchLengthParameter gblp : genealogyBranchLengthParameters) {
		double parameterMinimumValue = gblp.getMinimumValue();
		double parameterMaximumValue = FRACTION_OF_MAXIMUM_FOR_RANDOMIZATION_PURPOSES * gblp.getMaximumValue();
		double randomDistance = Math.random() * (parameterMaximumValue - parameterMinimumValue) + parameterMinimumValue;
		gblp.setValue(randomDistance, true, true, false);
	    }
	}	    

	if (enableSwitchingFrequencyOptimizationFlag) {
	    for (SwitchingFrequencyParameter sfp : switchingFrequencyParameters) {
		double parameterMinimumValue = sfp.getMinimumValue();
		double parameterMaximumValue = FRACTION_OF_MAXIMUM_FOR_RANDOMIZATION_PURPOSES * sfp.getMaximumValue();
		double randomProbability = Math.random() * (parameterMaximumValue - parameterMinimumValue) + parameterMinimumValue;
		sfp.setValue(randomProbability, true, true, false);
	    }
	}

	if (enableSubstitutionModelOptimizationFlag) {
	    for (GTRRateParameter grp : gtrRateParameters) {
		double parameterMinimumValue = grp.getMinimumValue();
		double parameterMaximumValue = FRACTION_OF_MAXIMUM_FOR_RANDOMIZATION_PURPOSES * grp.getMaximumValue();
		double randomDistance = Math.random() * (parameterMaximumValue - parameterMinimumValue) + parameterMinimumValue;
		grp.setValue(randomDistance, true, true, false);
	    }
	    
	    for (int i = 0; i < gtrBaseFrequencyParameters.size(); i++) {
		GTRBaseFrequencyParameter gbfp = gtrBaseFrequencyParameters.get(i);

		if (i == 0) {
		    // by convention
		    gbfp.setValue(1.0, true, true, false);
		}
		else {
		    double parameterMinimumValue = gbfp.getMinimumValue();
		    double parameterMaximumValue = FRACTION_OF_MAXIMUM_FOR_RANDOMIZATION_PURPOSES * gbfp.getMaximumValue();
		    double randomWeight = Math.random() * (parameterMaximumValue - parameterMinimumValue) + parameterMinimumValue;
		    gbfp.setValue(randomWeight, true, true, false);
		}
	    }
	}
	
	updateModelStateForEveryParameter();
    }

    /**
     * f 
     * No clone constructors for any of the network data structures.
     * Just cache and restore all branch lengths and probabilities by hand.
     * For child_node <- parent_node edge, cache key is the string "label(child_node) label(parent_node)".
     *
     * Need to pass in two empty maps.
     */
    protected void pushCacheValues () {
	for (Parameter parameter : parameterNameMap.values()) {
	    parameter.pushCacheValue();
	}
    }

    /**
     * Restore from a checkpoint file filename.
     * WARNING - doesn't perform an update after writing all Parameter object values! 
     * Caller must call updateModelStateForEveryParameter() after calling this function!
     */
    protected void restoreParameterValuesFromCheckpointWithNoUpdate (String filename) {
	Hashtable<String,Double> map = new Hashtable<String,Double>();

	// read checkpoint first, make sure everything is legit
	try {
	    BufferedReader br = new BufferedReader(new FileReader(filename));
	    String line;
	    while ((line = br.readLine()) != null) {
		StringTokenizer st = new StringTokenizer(line);
		// strict!
		if (st.countTokens() != 2) {
		    System.err.println ("Invalid line in checkpoint file " + filename + ". Not restoring from checkpoint.");
		    return;
		}

		try {
		    String name = st.nextToken();
		    Double value = Double.valueOf(st.nextToken());

		    // make sure that a Parameter object
		    // with the specified name exists.
		    if (!parameterNameMap.containsKey(name)) {
			System.err.println ("Invalid parameter name " + name + " in checkpoint file " + filename + " line " + line + ". Not restoring from checkpoint.");
			return;
		    }

		    // make sure no duplicate names in checkpoint file
		    if (map.containsKey(name)) {
			System.err.println ("Duplicate parameter name " + name + " in checkpoint file " + filename + " line " + line + ". Not restoring from checkpoint.");
			return;
		    }

		    map.put(name, value);
		}
		catch (NumberFormatException nfe) {
		    System.err.println (nfe);
		    nfe.printStackTrace();
		    System.err.println ("Invalid line in checkpoint file " + filename + " line " + line + ". Not restoring from checkpoint.");
		    return;
		}
	    }
	}
	catch (IOException ioe) {
	    System.err.println (ioe);
	    ioe.printStackTrace();
	    // no update
	    System.err.println ("ERROR: restoreCheckpoint(...) failed on checkpoint file " + filename + ". Not restoring from checkpoint.");
	    return;
	}

	// then do atomic update
	for (String name : map.keySet()) {
	    // kliu - no delay update until the end
	    parameterNameMap.get(name).setValue(map.get(name).doubleValue(), true, true, false);
	}
    }

    /**
     * Simple map. 
     * Everything goes through Parameter objects.
     *
     * Use tabs to delimit Parameter name/value pairs.
     *
     * Also need to cache pass/round counts. 
     */
    protected void writeCheckpoint (int pass, int round, double likelihood) {
	try {
	    // write out Parameter-based state
	    String filename = runHmmObject.getWorkingDirectory() + File.separator + CHECKPOINT_PARAMETER_STATE_FILENAME + FILENAME_SUFFIX_DELIMITER + Integer.toString(pass) + FILENAME_SUFFIX_DELIMITER + Integer.toString(round);
	    if (Constants.WARNLEVEL > 4) { System.out.println ("Writing checkpoint  " + filename + "."); }
	    BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
	    // safe due to guards in createNameMapOfAllParameters()
	    for (Parameter parameter : parameterNameMap.values()) {
		bw.write(parameter.getName() + CHECKPOINT_ENTRY_PAIR_DELIMITER_CHARACTER + parameter.getValue()); bw.newLine();
	    }
	    bw.flush();
	    bw.close();
	    if (Constants.WARNLEVEL > 4) { System.out.println ("Writing checkpoint  " + filename + " DONE."); }

	    // write Viterbi-optimal hidden state sequence
	    filename = runHmmObject.getWorkingDirectory() + File.separator + CHECKPOINT_VITERBI_SEQUENCE_FILENAME + FILENAME_SUFFIX_DELIMITER + Integer.toString(pass) + FILENAME_SUFFIX_DELIMITER + Integer.toString(round);
	    Tuple<int[],Double> viterbiResult = hmm.viterbiStateSequence(observation);
	    bw = new BufferedWriter(new FileWriter(filename));
	    double viterbiLLH = viterbiResult.Item2;
	    for (int index : viterbiResult.Item1) {
		bw.write(hiddenStates.get(index).getName()); bw.newLine();
	    }
	    bw.flush();
	    bw.close();

	    // write likelihood
	    filename = runHmmObject.getWorkingDirectory() + File.separator + CHECKPOINT_LIKELIHOOD_FILENAME + FILENAME_SUFFIX_DELIMITER + Integer.toString(pass) + FILENAME_SUFFIX_DELIMITER + Integer.toString(round);
	    if (Constants.WARNLEVEL > 4) { System.out.println ("Writing checkpoint  " + filename + "."); }
	    bw = new BufferedWriter(new FileWriter(filename));
	    // safe due to guards in createNameMapOfAllParameters()
	    bw.write(Double.toString(likelihood)); bw.newLine();
	    bw.write(Double.toString(viterbiLLH)); bw.newLine();
	    bw.flush();
	    bw.close();
	    
	    if (Constants.WARNLEVEL > 4) { System.out.println ("Writing checkpoint  " + filename + " DONE."); }
	}
	catch (IOException ioe) {
	    System.err.println (ioe);
	    ioe.printStackTrace();
	    // ok to return without writing checkpoint
	    // just lose a bit of progress - that's fine
	}
    }

    /**
     * WARNING: don't call this prior to a corresponding
     * MultivariateOptimizer.pushCacheValues() call!
     * Otherwise wonky behavior will result.
     */
    protected void popCacheValues (boolean setValueFlag) {
	for (Parameter parameter : parameterNameMap.values()) {
	    // delay Parameter.updateModelState() calls until all parameter values updated
	    // no need to check min/max values since all guards satisfied during full cache operation
	    parameter.popCacheValue(false, false, false, setValueFlag);
	}

	// finally, update internal state of HMM and propagate through to Jahmm's probability matrices/vectors
	if (setValueFlag) {
	    updateModelStateForEveryParameter();
	}
    }

    /**
     * WARNING - modifies HMM object and associated objects.
     * Returns likelihood of optimized model given observations.
     * 
     * WARNING - only InitialSearchSettings.CURRENT will use checkpoint state to begin search.
     */
    public double optimize (MultivariateOptimizer.InitialSearchSettings[] initialSearchSettings) {
	// initial search settings matter quite a bit
	// try fewer rounds, but a couple of random restarts
	// InitialSearchSettings[] initialSearchSettings = new InitialSearchSettings[] 
	//     { InitialSearchSettings.CURRENT, // think about how to get parameter settings in line with inputs???
	//       // // kliu - hold off on this now for testing purposes
	//       InitialSearchSettings.DEFAULT,
	//       //InitialSearchSettings.RANDOM,
	//       //InitialSearchSettings.RANDOM,
	//       InitialSearchSettings.RANDOM
	//     };

	// testing
	// InitialSearchSettings[] initialSearchSettings = new InitialSearchSettings[] 
	//     { InitialSearchSettings.RANDOM,
	//       InitialSearchSettings.RANDOM
	//     };

	double finalLikelihood = 1.0; // impossible log likelihood - why won't the compiler allow this to be uninitialized?

	for (int pass = 0; pass < initialSearchSettings.length; pass++) {
	    // cache all continuous model parameters
	    // could save a cache if we just did a restore - meh
	    //
	    // Maintain the invariant that branch-length-constraint-set constraints are
	    // satisfied prior to each cache operation.
	    // Thus, no need for cache to worry about branch-length-constraint-set weights.
	    //
	    // Crap - this makes checkpointing more difficult than it needs to be.
	    // A simpler way would just be to keep track of complete Parameter object state
	    // for solution with best likelihod so far.
	    // If current solution not as good as best, then revert.
	    // Just need to constantly checkpoint best and current solution (Parameter state),
	    // as well as iteration state.
	    // Meh - still too complicated.
	    // Max 1 day for 1 pass isn't too much of a restriction.
	    // DMTCP easiest anyways.
	    pushCacheValues();

	    double passLikelihood = singlePassOptimization(pass, initialSearchSettings[pass]);

	    // update
	    if ((finalLikelihood > 0.0) || // uninitialized 
            (passLikelihood > finalLikelihood)) { // improvement
		finalLikelihood = passLikelihood;
		// no need for restore
		popCacheValues(false);

		if (Constants.WARNLEVEL > 1) { 
		    System.out.println ("Updating with likelihood " + passLikelihood + ".");
		}
	    }
	    else {
		// need to restore
		popCacheValues(true);

		if (Constants.WARNLEVEL > 1) { 
		    System.out.println ("Not updating.");
		}
	    }
	}

	// String resultString = 
	//     "Final network: |" + getNetworkString() + "|\n" +
	//     "Final log likelihood: |" + finalLikelihood + "|\n";

	// return output string later
	return (finalLikelihood);
    }




    /**
     * Other than min/max, no difference between branch lengths and frequencies.
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
	    parameter.pushCacheValue();

	    //double originalSetting = parameter.getValue();

	    // update
	    parameter.setValue(x);

	    // kliu - hmm... updateHmm is *WAYYYY* too slow
	    // need to cache computations as much as possible - only 
	    // recompute when necessary
	    // *especially* for probability of gene tree given parental tree P[g|T]

	    // update associated model values associated with a single parameter
	    // e.g., for LengthParameter objects that belong to a ParameterConstraintSet
	    // or multiple branches share a LengthParameter
	    // or multiple parental trees have branches that share a LengthParameter
	    // etc.
	    // 
	    // don't push this into Parameter
	    // too complicated
	    // updateHMM(parameter); // migrated to setValue

	    // evaluate f(x) using forward/backwards algorithm
	    double result = computeHMMLikelihood();

	    // restore original setting
	    // use all guards, all set ops
	    //
	    // fine to check min/max - only change one Parameter at a time,
	    // and this Parameter's min/max depends on at most one *other* Parameter
	    // object's value.
	    parameter.popCacheValue(true, true, true, true);

	    //parameter.setValue(originalSetting);
	    //updateHMM(parameter); migrated to setValue
	    

	    // above restore original setting op
	    // not really necessary if always accept new branch length even if no likelihood improvement
	    // since Brent's method is guaranteed to never make the likelihood worse
	    // meh, no update if no likelihood improvement
	    // -> above is necessary

	    return (result);
	}	
    }
    

    public class IllegalSearchIntervalException extends RuntimeException {
	public IllegalSearchIntervalException (String message) {
	    super(message);
	}
    }

    /**
     * Signal error from f() evaluation, since no nice way to return meaningful messages from UnivariateFunction.value().
     */
    public class UnivariateFunctionEvaluationException extends RuntimeException {
	public UnivariateFunctionEvaluationException (String message) {
	    super(message);
	}
    }

}




// Shoot - need to name everything.
// Since might need to fix gene genealogy branch lengths.
// Also neater this way.
// This means that all trees must be named.


    // protected void updateHMM (Parameter parameter) {
    // 	updateHMM(parameter, false);
    // }

    // /**
    //  * Update associated model values associated with a single parameter
    //  * e.g., for LengthParameter objects that belong to a ParameterConstraintSet
    //  * or multiple branches share a LengthParameter
    //  * or multiple parental trees have branches that share a LengthParameter
    //  * etc.
    //  * 
    //  * don't push this into Parameter
    //  * too complicated
    //  *
    //  * Don't set disableHMMProbabilityUpdateFlag to true unless you 
    //  * know what you're doing. Otherwise can cause inconsistent state
    //  * between our HMM state and Jahmm's probability matrices/vectors.
    //  */
    // protected void updateHMM (Parameter parameter, boolean disableHMMProbabilityUpdateFlag) {

    // 	else {
    // 	    throw (new RuntimeException("ERROR: parameter in updateHMM(...) has unrecognized type."));
    // 	}

    // 	// need to propagate changes on to Jahmm's probability matrices/vectors
    // 	if (!disableHMMProbabilityUpdateFlag) {
    // 	    runHmmObject.updateTransitionProbabilities();
    // 	}
    // }
