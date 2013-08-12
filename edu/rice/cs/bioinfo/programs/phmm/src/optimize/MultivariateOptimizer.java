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
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.library.programming.BijectiveHashtable;
import edu.rice.cs.bioinfo.library.programming.BidirectionalMultimap;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;


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

    // only for randomization purposes
    //
    // only go to 80% of max for randomization purposes
    public static final double FRACTION_OF_MAXIMUM_FOR_RANDOMIZATION_PURPOSES = 0.8;

    // hmm... is GeneTreeProbability able to handle probability == 0 or 1?
    // yes, it handles this fine
    public static final double DEFAULT_MINIMUM_PROBABILITY = 0.0;
    public static final double DEFAULT_INITIAL_PROBABILITY = 1e-2;
    // initial probability always initialized to uniform
    public static final double DEFAULT_MAXIMUM_PROBABILITY = 1.0;

    // hack - to support set branch length constraints
    // PhyloNet's Rich Newick support has some craziness about '_' underscore letter.
    // use dash '-' instead.
    public static final String IDENTIFIER_SET_BRANCH_LENGTH_CONSTRAINTS = "-CONSTRAINT-SET-";

    public static final String PARENTAL_NODE_LABEL_DELIMITER = ",";

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
    protected Map<Network<Double>,Set<HiddenState>> parentalTreeClasses;
    protected BijectiveHashtable<String,Network<Double>> parentalTreeNameMap;
    protected List<ObservationMap> observation;

    // bijective map between parental tree node objects and their names
    // naming convention:
    // "<tree name>,<node name>"
    protected BijectiveHashtable<NetNode<Double>,String> parentalNodeLabelMap;

    // Use a single BidirectionalMultimap to capture many-to-many 
    // relationship between
    // length parameters and edge ids:
    // lid <-> eid
    protected BidirectionalMultimap<LengthParameter,Tuple<String,String>> lpEidMap;
    // reverse lookup from lid name to length-parameter references
    protected Hashtable<String,LengthParameter> lidLpMap;

    // One more layer of indirection.
    // Group length parameters into length-parameter-constraint-sets.
    // Each length-parameter-constraint-set has a total length constraint.
    // If length-parameter belongs to a length-parameter-constaint-set, then
    // enable fixed-total-length optimization for the set.
    protected BidirectionalMultimap<ParameterConstraintSet,LengthParameter> setLpMap;
    // By convention, first length-parameter in a length-parameter-constraint-set 
    // has relative weight 1.0 and isn't optimized.
    // Beware degree-of-freedoms == one less than number of length-parameters in
    // a length-parameter-constraint-set.
    //
    // Store this info in sid->lid map.
    protected Hashtable<ParameterConstraintSet,LengthParameter> setFirstLpMap;
    
    // kliu - implement this later
    protected boolean enableParentalTreeOptimizationFlag = true;
    protected boolean enableGeneGenealogyOptimizationFlag = true;
    protected boolean enableFrequencyOptimizationFlag = true;

    // heuristic for univariate optimization
    protected BrentOptimizer brentOptimizer;

    // keep list of Parameters
    protected Vector<Parameter> parameters;
    // use lpEidMap.keys() to get a list of all parental tree branch length parameters (LengthParameter objects)
    protected Vector<SingleBranchLengthParameter> singleBranchLengthParameters;
    protected Vector<FrequencyParameter> frequencyParameters;

    public MultivariateOptimizer (Hmm<ObservationMap> inHmm,
				  // HMM update routines are located in this object
				  runHmm inRunHmm,
				  List<HiddenState> inHiddenStates,
				  TransitionProbabilityParameters inTransitionProbabilityParameters,
				  Map<Network<Double>,Set<HiddenState>> inParentalTreeClasses,
				  BijectiveHashtable<String,Network<Double>> inParentalTreeNameMap,
				  List<ObservationMap> inObservation,
				  String inputLengthParameterToEdgeMapFilename,
				  String inputLengthParameterSetConstraintsFilename
				  ) {
	this.hmm = inHmm;
	this.runHmmObject = inRunHmm;
	this.hiddenStates = inHiddenStates;
	this.transitionProbabilityParameters = inTransitionProbabilityParameters;
	this.parentalTreeClasses = inParentalTreeClasses;
	this.parentalTreeNameMap = inParentalTreeNameMap;
	this.observation = inObservation;

	brentOptimizer = new BrentOptimizer(RELATIVE_ACCURACY, ABSOLUTE_ACCURACY);

	verifySearchSettings();

	createParentalNodeLabelMap();
	
	// this creates LengthParameter objects representing parental tree branch length parameters
	parseInputLengthParameterToEdgeMapFile(inputLengthParameterToEdgeMapFilename);
	parseParameterConstraintSetsFile(inputLengthParameterSetConstraintsFilename);

	createSingleBranchLengthParameters();
	createFrequencyParameters();

	createListOfAllParameters();
    }

    protected void createFrequencyParameters () {
	frequencyParameters = new Vector<FrequencyParameter>();
	for (TransitionProbabilityParameters.ParameterChoice parameterChoice : TransitionProbabilityParameters.ParameterChoice.values()) {
	    FrequencyParameter fp = new FrequencyParameter(parameterChoice.toString(),
							   transitionProbabilityParameters.get(parameterChoice),
							   transitionProbabilityParameters,
							   parameterChoice,
							   runHmmObject,
							   false // no need to update
							   );
	    frequencyParameters.add(fp);
	}
    }

    /**
     * In current model, hidden states share parental trees but not
     * gene trees. Each hidden state has a gene tree that is parameterized
     * separately from any other gene tree in the HMM.
     */
    protected void createSingleBranchLengthParameters () {
	singleBranchLengthParameters = new Vector<SingleBranchLengthParameter>();
	for (HiddenState hiddenState : hiddenStates) {
	    for (TNode node : hiddenState.getGeneGenealogy().postTraverse()) {
		if (node.isRoot() || node.getParent() == null) {
		    continue;
		}
		
		SingleBranchLengthParameter sblp = new SingleBranchLengthParameter(hiddenState.getName() + HiddenState.HIDDEN_STATE_NAME_DELIMITER + node.getName(),
										   node.getParentDistance(),
										   node,
										   // no need for update
										   false);
		singleBranchLengthParameters.add(sblp);
	    }
	}
    }

    protected void createListOfAllParameters () {
	parameters = new Vector<Parameter>();

	// parental tree branch length parameters
	for (LengthParameter lp : lpEidMap.keys()) {
	    parameters.add(lp);
	}
	
	for (SingleBranchLengthParameter sblp : singleBranchLengthParameters) {
	    parameters.add(sblp);
	}

	for (FrequencyParameter fp : frequencyParameters) {
	    parameters.add(fp);
	}
    }

    /**
     * Token format: <branch parental tree name>,<branch child id>,<branch parent id>
     */
    protected Tuple<String,String> parseEid (String filename, String s) {
	StringTokenizer st = new StringTokenizer(s, ",");
	if (st.countTokens() != 3) {
	    throw (new RuntimeException("ERROR: incorrectly formatted edge id in file " + filename + ": " + s));
	}

	String tid = st.nextToken();
	String cid = st.nextToken();
	String pid = st.nextToken();
	return (new Tuple<String,String>(tid + PARENTAL_NODE_LABEL_DELIMITER + cid, tid + PARENTAL_NODE_LABEL_DELIMITER + pid));
    }

    /**
     * Input file format:
     * <length-parameter-constraint-set ID> <total length> <length parameter 1 unique ID> <length parameter 2 unique ID> ...
     *
     * A length parameter ID can appear AT MOST ONCE in this input file!
     * Avoid complex constraints - sets *MUST* be disjoint!
     */
    protected void parseParameterConstraintSetsFile (String filename) {
	if ((filename == null) ||
	    (filename.trim().length() <= 0)) {
	    throw (new RuntimeException ("ERROR: invalid filename in parseParameterConstraintSetsFile()."));
	}

	if (!(new File(filename)).isFile()) {
	    throw (new RuntimeException ("ERROR: filename " + filename + " does not exist."));
	}

	setLpMap = new BidirectionalMultimap<ParameterConstraintSet,LengthParameter>();
	setFirstLpMap = new Hashtable<ParameterConstraintSet,LengthParameter>();

	try {
	    BufferedReader br = new BufferedReader(new FileReader(filename));
	    String line;
	    while ((line = br.readLine()) != null) {
		// skip empty lines
		line = line.trim();
		if (line.length() <= 0) {
		    continue;
		}

		StringTokenizer st = new StringTokenizer(line);
		if (st.countTokens() < 3) {
		    throw (new RuntimeException ("ERROR: incorrectly formatted line in filename " + filename + ": " + line));
		}

		String sid = st.nextToken();
		double value = Double.parseDouble(st.nextToken());
		ParameterConstraintSet slp = new LengthParameter(sid, value);
		
		// Disallow duplicate length-parameter-ids in input file.
		if (setLpMap.containsKey(slp)) {
		    throw (new RuntimeException ("ERROR: duplicate length-parameter-constraint-set id in filename " + filename + ": " + slp.toString()));
		}

		while (st.hasMoreTokens()) {
		    String lid = st.nextToken();
		    if (!lidLpMap.containsKey(lid)) {
			// barf
			throw (new RuntimeException("ERROR: invalid lid " + lid + " in file " + filename + "."));
		    }
		    LengthParameter lp = lidLpMap.get(lid);
		    // duplicate lids fine
		    // all goes into a HashSet anyways
		    setLpMap.put(slp, lp);

		    // Keep track of first lid in set.
		    if (!setFirstLpMap.containsKey(slp)) {
			setFirstLpMap.put(slp, lp);
			// Weight of first lid is always 1.0.
			lp.setValue(1.0);
		    }
		}
	    }
	}
	catch (IOException ioe) {
	    System.err.println (ioe);
	    System.exit(1);
	}

	// One last constraint - each length-parameter can belong to *ONLY* 
	// one length-parameter-constraint-set.
	// No complex constraints.
	if (!setLpMap.checkInjective()) {
	    throw (new RuntimeException ("ERROR: injective property violated in input file + " + filename + ". Make sure that each length-parameter belongs to at most one length-parameter-constraint-set."));
	}

	if (Constants.WARNLEVEL > 4) { System.out.println ("setLpMap at end of parseParameterConstraintSetsFile(): \n" + setLpMap.toString() + "\n"); }
    }

    /**
     * Input file format:
     * <length parameter unique ID> <initial value> <branch 1 parental tree name>,<branch 1 child id>,<branch 1 parent id> <branch 2 parental tree name>,<branch 2 child id>,<branch 2 parent id> ...
     */
    protected void parseInputLengthParameterToEdgeMapFile (String filename) {
	if ((filename == null) ||
	    (filename.trim().length() <= 0)) {
	    throw (new RuntimeException ("ERROR: invalid filename in parseInputLengthParameterToEdgeMapFile()."));
	}

	if (!(new File(filename)).isFile()) {
	    throw (new RuntimeException ("ERROR: filename " + filename + " does not exist."));
	}

	lpEidMap = new BidirectionalMultimap<LengthParameter,Tuple<String,String>>();
	lidLpMap = new Hashtable<String,LengthParameter>();

	try {
	    BufferedReader br = new BufferedReader(new FileReader(filename));
	    String line;
	    while ((line = br.readLine()) != null) {
		// skip empty lines
		line = line.trim();
		if (line.length() <= 0) {
		    continue;
		}

		StringTokenizer st = new StringTokenizer(line);
		if (st.countTokens() < 3) {
		    throw (new RuntimeException ("ERROR: incorrectly formatted line in filename " + filename + ": " + line));
		}
		String lid = st.nextToken();
		double value = Double.parseDouble(st.nextToken());
		LengthParameter lp = new LengthParameter(lid, value);
		
		//if (Constants.WARNLEVEL > 4) { System.out.println ("Parsing " + filename + " lid: |" + lid + "| value: |" + value); }

		lidLpMap.put(lp.getName(), lp);

		// Disallow duplicate length-parameter-ids in input file.
		if (lpEidMap.containsKey(lp)) {
		    throw (new RuntimeException ("ERROR: duplicate length-parameter-id in filename " + filename + ": " + lp.toString()));
		}
		
		while (st.hasMoreTokens()) {
		    Tuple<String,String> eid = parseEid (filename, st.nextToken());

		    //if (Constants.WARNLEVEL > 4) { System.out.println ("Parsing " + filename + " lid: |" + lid + "| value: |" + value + "| eid: |" + eid + "|"); }

		    // Repeated eid per length-parameter ok.
		    // All goes into a set anyways (keyed on length-parameter).
		    lpEidMap.put(lp, eid);
		}
	    }
	}
	catch (IOException ioe) {
	    System.err.println (ioe);
	    System.exit(1);
	}

	// force all eids to exist
	// furthermore, force all edges to be listed in lpEidMap
	if (!checkLpEidMap()) {
	    throw (new RuntimeException("ERROR: invalid lpEidMap in parseInputLengthParameterToEdgeMapFile."));
	}

	if (Constants.WARNLEVEL > 4) { System.out.println ("lpEidMap at end of parseInputLengthParameterToEdgeMapFile: \n" + lpEidMap.toString() + "\n"); }
    }

    protected boolean checkLpEidMap () {
	// also make sure that every edge in the tree is contained in lpEidMap
	// 
	// actually - not necessary
	// excluded edges -> not optimized
	//
	// Just warn.
	for (Network<Double> parentalTree : parentalTreeClasses.keySet()) {
	    for(NetNode<Double> node : parentalTree.dfs()) {
		for(NetNode<Double> parent : node.getParents()) {
		    Tuple<String,String> eid = new Tuple<String,String>(parentalNodeLabelMap.get(node), parentalNodeLabelMap.get(parent));
		    
		    // testing
		    System.out.println ("Edge: |" + eid + "|");
		    
		    if (!lpEidMap.containsValue(eid)) {
			if (Constants.WARNLEVEL > 1) { System.err.println ("WARNING: no map entry for edge with id " + eid + ". No length optimization will be performed for this edge."); }
			//return (false);
		    }
		}
	    }
	}

	for (Tuple<String,String> eid : lpEidMap.values()) {
	    if (!checkEid(eid)) {
		if (Constants.WARNLEVEL > 1) { System.err.println ("ERROR: invalid map entry with edge id " + eid + "."); }
		return (false);
	    }
	}

	return (true);
    }

    protected boolean checkEid (Tuple<String,String> eid) {
	if (!(parentalNodeLabelMap.containsValue(eid.Item1) && parentalNodeLabelMap.containsValue(eid.Item2))) {
	    return (false);
	}

	// guaranteed bijective by guard in constructor
	NetNode<Double> c = parentalNodeLabelMap.rget(eid.Item1);
	NetNode<Double> p = parentalNodeLabelMap.rget(eid.Item2);

	// argh
	for (NetNode<Double> cp : c.getParents()) {
	    if (checkNodesEqual(p, cp)) {
		return (true);
	    }
	}
	
	return (false);
    }

    /**
     * annoying
     * kludged equals() function.
     * Warning - throws RuntimeException if node<->label map lookup fails.
     */
    protected boolean checkNodesEqual (NetNode<Double> x, NetNode<Double> y) {
	if (parentalNodeLabelMap == null) {
	    throw (new NodeEqualityTestException("ERROR: node<->label map not initialized in checkNodesEqual()."));
	}

	if (!parentalNodeLabelMap.containsKey(x)) {
	    throw (new NodeEqualityTestException("ERROR: node<->label map doesn't contain node " + x.getName() + "."));
	}

	if (!parentalNodeLabelMap.containsKey(y)) {
	    throw (new NodeEqualityTestException("ERROR: node<->label map doesn't contain node " + y.getName() + "."));
	}

	return (parentalNodeLabelMap.get(x).equals(parentalNodeLabelMap.get(y)));
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

    protected void createParentalNodeLabelMap () {
	parentalNodeLabelMap = new BijectiveHashtable<NetNode<Double>,String>();
	for (String parentalTreeName : parentalTreeNameMap.keys()) {
	    Network<Double> parentalTree = parentalTreeNameMap.get(parentalTreeName);
	    for (NetNode<Double> node : parentalTree.dfs()) {
		String label = parentalTreeName + PARENTAL_NODE_LABEL_DELIMITER + node.getName();
		// strict!
		if (parentalNodeLabelMap.containsKey(label)) {
		    throw (new RuntimeException("ERROR: duplicate label " + label + " in parental node <-> label map. Check inputs for duplicate parental tree names and/or parental node names."));
		}
		parentalNodeLabelMap.put(node, label);
	    }
	}
    }

    /**
     * Handles all parameters.
     */
    protected void updateHMM () {
	// duplicated work here for each member of a length-parameter-constraint-set
	// oh well
	// ok since each call below just rebalances weight to satisfy constraint on a constraint-set
	// repeated re-balancings don't change anything, so long as individual
	// parameter values are unchanged
	for (Parameter p : parameters) {
	    // manually update Jahmm's probability matrices/vectors *once* at the end
	    updateHMM(p, true);
	}

	// see above comment
	runHmmObject.updateTransitionProbabilities();
    }

    protected void updateHMM (Parameter parameter) {
	updateHMM(parameter, false);
    }

    /**
     * Update associated model values associated with a single parameter
     * e.g., for LengthParameter objects that belong to a ParameterConstraintSet
     * or multiple branches share a LengthParameter
     * or multiple parental trees have branches that share a LengthParameter
     * etc.
     * 
     * don't push this into Parameter
     * too complicated
     *
     * Don't set disableHMMProbabilityUpdateFlag to true unless you 
     * know what you're doing. Otherwise can cause inconsistent state
     * between our HMM state and Jahmm's probability matrices/vectors.
     */
    protected void updateHMM (Parameter parameter, boolean disableHMMProbabilityUpdateFlag) {
	if (parameter instanceof SingleBranchLengthParameter) {
	    SingleBranchLengthParameter sblp = (SingleBranchLengthParameter) parameter;
	    sblp.updateModelState();
	}
	else if (parameter instanceof FrequencyParameter) {
	    FrequencyParameter fp = (FrequencyParameter) parameter;
	    fp.updateModelState();
	}
	else if (parameter instanceof LengthParameter) {
	    LengthParameter inputLengthParameter = (LengthParameter) parameter;
	    // do appropriate multiplexing here
	    if (setLpMap.containsValue(inputLengthParameter)) {
		ParameterConstraintSet slp = getParameterConstraintSet(inputLengthParameter);
		for (LengthParameter member : setLpMap.get(slp)) {
		    updateLengthParameterBranchLengthsHelper(member);
		}
	    }
	    else { 
		// not a member of a constraint-set
		updateLengthParameterBranchLengthsHelper(inputLengthParameter);
	    }
	}
	else {
	    throw (new RuntimeException("ERROR: parameter in updateHMM(...) has unrecognized type."));
	}

	// need to propagate changes on to Jahmm's probability matrices/vectors
	if (!disableHMMProbabilityUpdateFlag) {
	    runHmmObject.updateTransitionProbabilities();
	}
    }

    /**
     * Compute the current value for a member of a 
     * length-parameter-constraint-set.
     * Would be more elegant to do some sort of caching.
     */
    protected double computeValueForMemberOfParameterConstraintSet (LengthParameter lp) {
	// paranoid
	if (!setLpMap.containsValue(lp)) {
	    throw (new RuntimeException ("ERROR: computeValueForMemberOfParameterConstraintSet() called with length-parameter that is not a member of a length-parameter-constraint-set."));
	}

	// sum up weights from members of length-parameter-constraint-set
	ParameterConstraintSet slp = getParameterConstraintSet(lp);
	double totalWeight = 0.0;
	for (LengthParameter members : setLpMap.get(slp)) {
	    // no need to any recursion here due to
	    // injective property of setLpMap
	    //
	    // Each length-parameter belongs to *at most* one length-parameter-constraint-set.
	    // All members uniquely belong to this length-parameter-constraint-set.
	    totalWeight += members.getValue();
	}
	// solve for x
	double frac = slp.getValue() / totalWeight;

	// substitute
	return (lp.getValue() * frac);
    }

    /**
     * Update all network branch lengths associated with a length-parameter.
     * If length-parameter belongs to a length-parameter-constraint-set,
     * need to maintain total length in set.
     * 
     * WARNING: only updates inputLengthParameter!
     */
    protected void updateLengthParameterBranchLengthsHelper (LengthParameter inputLengthParameter) {
	// if length-parameter is a member of a length-parameter-constraint-set
	// need to "atomically" update all members in the set

	//if (Constants.WARNLEVEL > 4) { System.out.println ("updateNetworkBranchLengthsHelper() for length-parameter with id " + inputLengthParameter.getName() + "."); }

	// Get relevant branches.
	// Weird? No branches returned from map?
	for (Tuple<String,String> eid : lpEidMap.get(inputLengthParameter)) {

	    //if (Constants.WARNLEVEL > 4) { System.out.println ("Processing eid: " + eid.toString()); }

	    // Update relevant branch.
	    if (!checkEid(eid)) {
		throw (new RuntimeException("ERROR: unknown node labels in eid " + eid + " in updateNetworkBranchLengths()."));
	    } 

	    NetNode<Double> child = parentalNodeLabelMap.rget(eid.Item1);
	    NetNode<Double> parent = parentalNodeLabelMap.rget(eid.Item2);
	    double length = 0.0;
	    // simple addition, nothing fancy
	    for (LengthParameter lp : lpEidMap.rget(eid)) {
		// if length-parameter belongs to length-parameter-constraint-set,
		// no storage - just re-compute value based on weighted ratios
		// each time
		if (setLpMap.containsValue(lp)) {
		    length += computeValueForMemberOfParameterConstraintSet(lp);
		}
		else {
		    // straightforward for length-parameters
		    // that aren't subject to length-parameter-constraint-sets
		    length += lp.getValue();
		}
	    }
	    child.setParentDistance(parent, length);

	    //System.out.println ("setParentDistance: " + child.getName() + " " + parent.getName() + " " + length);
	}

	//if (Constants.WARNLEVEL > 4) { System.out.println ("updateNetworkBranchLengthsHelper() for length-parameter with id " + inputLengthParameter.getName() + " DONE."); }
    }

    /**
     * Get length-parameter-constraint-set for a length-parameter.
     * Enforces injective property of setLpMap via RuntimeExceptions.
     * Returns null if length-parameter doesn't belong to a length-parameter-constraint-set.
     */
    protected ParameterConstraintSet getParameterConstraintSet (LengthParameter lp) {
	if (!setLpMap.containsValue(lp)) {
	    return (null);
	}

	// paranoid
	// despite injective constraint on relation
	if (setLpMap.rget(lp).size() != 1) {
	    throw (new RuntimeException("ERROR: injective property violated in setLpMap for length-parameter with id " + lp.getName() + "."));
	}
	ParameterConstraintSet slp = setLpMap.rget(lp).iterator().next();
	
	return (slp);
    }

    /**
     * Convenience function.
     * Check to see if a length-parameter is the first in a
     * length-parameter-constraint-set.
     */
    protected boolean checkFirstInParameterConstraintSet (LengthParameter lp) {
	if (setLpMap.containsValue(lp)) {
	    ParameterConstraintSet slp = getParameterConstraintSet(lp);
	    // paranoid
	    if (!setFirstLpMap.containsKey(slp)) {
		throw (new RuntimeException ("ERROR: setFirstLpMap doesn't contain sid " + slp.getName() + "."));
	    }
	    
	    if (setFirstLpMap.get(slp).equals(lp)) {
		return (true);
	    }
	}

	return (false);
    }

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
	// tweak these and x -> output
	double l;
	double u;

	if (Constants.WARNLEVEL > 4) { System.out.println (debugMessage + " search interval min and max: |" + min + " " + max + "|"); }

	// due to constructor checks, guaranteed
	// that min/max interval width is large enough for search-by-halves to work

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
	     "optimizeSingleLengthParameter: ");
	UnivariatePointValuePair upvp = brentOptimizer.optimize(BRENT_METHOD_SINGLE_ROUND_MAXIMUM_ITERATIONS,
								f,
								GoalType.MAXIMIZE,
								searchInterval.Item1.doubleValue(), // l
								searchInterval.Item3.doubleValue(), // u
							        searchInterval.Item2.doubleValue()); // x

	if (Constants.WARNLEVEL > 4) { System.out.println ("Brent optimization point: |" + upvp.getPoint() + "| likelihood: |" + upvp.getValue() + "|"); }

	double brentOptimizedLogLikelihood = upvp.getValue();
	// see function comments
	if (brentOptimizedLogLikelihood > logLikelihood) {
	    if (Constants.WARNLEVEL > 4) { System.out.println ("INFO: Brent's method resulted in strict improvement in likelihood. Updating."); }

	    // update 
	    p.setValue(upvp.getPoint());
	    updateHMM(p);
	    logLikelihood = brentOptimizedLogLikelihood;
	}
	else {
	    // no update - info instead
	    if (Constants.WARNLEVEL > 4) { System.out.println ("INFO: Round " + round + " optimized point and log likelihood for length parameter " + p.getName() + " resulted in log likelihood " + brentOptimizedLogLikelihood + " which isn't better than current log likelihood " + logLikelihood + ". Not updating branch length nor best round log likelihood."); }
	}

	return (logLikelihood);
    }

    /**
     * For diagnostic purposes.
     */
    // protected void debugParameters () {
    // 	for (Parameter p : parameters) {
    // 	    System.out.println (p.toString());
    // 	}
    // }

    /**
     * For diagnostic purposes.
     */
    protected void debugModel () {
	System.out.println ("============================================");
	System.out.println ("Current PhyloNet-HMM state: ");
	System.out.println ();
	System.out.println ("Hidden states: ");
	for (int i = 0; i < hiddenStates.size(); i++) {
	    System.out.println("Hidden state " + i + ":");
	    System.out.println(hiddenStates.get(i).toString());
	}
	System.out.println ();
	System.out.println ("Parameters: ");
	for (Parameter parameter : parameters) {
	    System.out.println (parameter.toString());
	}
	System.out.println ();
	System.out.println ("HMM transition and emission probabilities: ");
	System.out.println (hmm.toString());
	System.out.println ("============================================");
    }

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
	    debugModel();
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
	    debugModel();
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
	double roundLogLikelihood = initialLogLikelihood;
	while (!checkConvergence(round, prevLogLikelihood, currLogLikelihood)) {
	    if (Constants.WARNLEVEL > 1) { System.out.println ("Processing round " + round + "."); }

	    // iterate through parameters
	    // may want to make order random??
	    int pCount = 0;
	    for (Parameter p : parameters) {
		// by convention, skip the first lengthParameter
		// in a length-parameter-constraint-set
		if ((p instanceof LengthParameter) && 
		    (!(p instanceof SingleBranchLengthParameter)) && 
		    checkFirstInParameterConstraintSet((LengthParameter) p)) {
		    continue;
		}

		// skip parental tree branch length optimization if disabled
		if (!enableParentalTreeOptimizationFlag && 
		    (p instanceof LengthParameter) &&
		    (!(p instanceof SingleBranchLengthParameter))) {
		    continue;
		}

		// skip gene genealogy branch length optimization if disabled
		if (!enableGeneGenealogyOptimizationFlag &&
		    (p instanceof SingleBranchLengthParameter)) {
		    continue;
		}

		// skip frequency optimization if disabled
		if (!enableFrequencyOptimizationFlag &&
		    (p instanceof FrequencyParameter)) {
		    continue;
		}

		if (Constants.WARNLEVEL > 1) { System.out.println ("Processing parameter " + p.getName() + " count " + pCount + "."); }

		// single branch length optimization
		// does update appropriately
		    roundLogLikelihood = optimizeSingleParameter(p, roundLogLikelihood, round);

		if (Constants.WARNLEVEL > 1) { System.out.println ("Processing length parameter " + p.getName() + " count " + pCount + " DONE."); }
		pCount++;
	    }

	    // paranoid
	    if (Double.isNaN(roundLogLikelihood)) {
		System.err.println ("ERROR: unable to evaluate round " + round + " log likelihood. Aborting and returning NaN.");
		return (Double.NaN);
	    }

	    // end of round
	    // update likelihoods
	    // if we did the swap at the start of each iteration, wouldn't need the third variable - meh
	    prevLogLikelihood = currLogLikelihood;
	    currLogLikelihood = roundLogLikelihood;

	    if (Constants.WARNLEVEL > 1) { System.out.println ("Processing round " + round + " DONE."); }

	    round++;
	}

	// recompute for debugging output
	double finalLikelihood = computeHMMLikelihood();

	if (Constants.WARNLEVEL > 1) { 
	    debugModel();
	    System.out.println ("Pass log likelihood: |" + finalLikelihood + "|");	    
	    // // Also print out length-parameters.
	    // System.out.println ("Pass parameters:");
	    // debugParameters();
	    System.out.println ("Processing pass " + pass + " with initial search setting " + initialSearchSettings.toString() + " DONE.");
	}

	return (finalLikelihood);
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
	    for (LengthParameter parentalLengthParameter : lpEidMap.keys()) {
		if (checkFirstInParameterConstraintSet(parentalLengthParameter)) {
		    // canonical!
		    parentalLengthParameter.setValue(1.0);
		}
		else {
		    parentalLengthParameter.setValue(parentalLengthParameter.getDefaultInitialValue());
		}
	    }
	}

	if (enableGeneGenealogyOptimizationFlag) {
	    for (SingleBranchLengthParameter sblp : singleBranchLengthParameters) {
		sblp.setValue(sblp.getDefaultInitialValue());
	    }
	}	    

	if (enableFrequencyOptimizationFlag) {
	    for (FrequencyParameter fp : frequencyParameters) {
		fp.setValue(fp.getDefaultInitialValue());
	    }
	}

	updateHMM();
    }

    /**
     * Can't let all branch lengths approach max.
     * Can cause model likelihood calculation to barf.
     * Cap the maximum randomized branch length.
     */
    protected void initializeRandom () {
	if (enableParentalTreeOptimizationFlag) {
	    for (LengthParameter parentalLengthParameter : lpEidMap.keys()) {
		if (checkFirstInParameterConstraintSet(parentalLengthParameter)) {
		    // canonical!
		    parentalLengthParameter.setValue(1.0);
		}
		else {
		    double parameterMinimumValue = parentalLengthParameter.getMinimumValue();
		    double parameterMaximumValue = FRACTION_OF_MAXIMUM_FOR_RANDOMIZATION_PURPOSES * parentalLengthParameter.getMaximumValue();
		    double randomDistance = Math.random() * (parameterMaximumValue - parameterMinimumValue) + parameterMinimumValue;
		    parentalLengthParameter.setValue(randomDistance);
		}
	    }
	}

	if (enableGeneGenealogyOptimizationFlag) {
	    for (SingleBranchLengthParameter sblp : singleBranchLengthParameters) {
		double parameterMinimumValue = sblp.getMinimumValue();
		double parameterMaximumValue = FRACTION_OF_MAXIMUM_FOR_RANDOMIZATION_PURPOSES * sblp.getMaximumValue();
		double randomDistance = Math.random() * (parameterMaximumValue - parameterMinimumValue) + parameterMinimumValue;
		sblp.setValue(randomDistance);
	    }
	}	    

	// warning - random number may be larger than permitted
	if (enableFrequencyOptimizationFlag) {
	    for (FrequencyParameter fp : frequencyParameters) {
		double parameterMinimumValue = fp.getMinimumValue();
		double parameterMaximumValue = FRACTION_OF_MAXIMUM_FOR_RANDOMIZATION_PURPOSES * fp.getMaximumValue();
		double randomProbability = Math.random() * (parameterMaximumValue - parameterMinimumValue) + parameterMinimumValue;
		fp.setValue(randomProbability);
	    }
	}

	updateHMM();
    }

    /**
     * f 
     * No clone constructors for any of the network data structures.
     * Just cache and restore all branch lengths and probabilities by hand.
     * For child_node <- parent_node edge, cache key is the string "label(child_node) label(parent_node)".
     *
     * Need to pass in two empty maps.
     */
    protected void cacheParameterValues () {
	for (Parameter parameter : parameters) {
	    parameter.cacheValue();
	}
    }

    protected void restoreParameterValues () {
	for (Parameter parameter : parameters) {
	    parameter.restoreCachedValue();
	}

	// finally, update internal state of HMM and propagate through to Jahmm's probability matrices/vectors
	updateHMM();
    }

    /**
     * WARNING - modifies HMM object and associated objects.
     * Returns likelihood of optimized model given observations.
     */
    public double optimize () {
	// initial search settings matter quite a bit
	// try fewer rounds, but a couple of random restarts
	InitialSearchSettings[] initialSearchSettingsForEachPass = new InitialSearchSettings[] 
	    { InitialSearchSettings.CURRENT, // think about how to get parameter settings in line with inputs???
	      InitialSearchSettings.DEFAULT,
	      InitialSearchSettings.RANDOM,
	      InitialSearchSettings.RANDOM,
	      InitialSearchSettings.RANDOM
	    };

	// testing
	// InitialSearchSettings[] initialSearchSettingsForEachPass = new InitialSearchSettings[] 
	//     { InitialSearchSettings.RANDOM,
	//       InitialSearchSettings.RANDOM
	//     };

	double finalLikelihood = 1.0; // impossible log likelihood - why won't the compiler allow this to be uninitialized?

	for (int pass = 0; pass < initialSearchSettingsForEachPass.length; pass++) {
	    // cache all continuous model parameters
	    // could save a cache if we just did a restore - meh
	    //
	    // Maintain the invariant that branch-length-constraint-set constraints are
	    // satisfied prior to each cache operation.
	    // Thus, no need for cache to worry about branch-length-constraint-set weights.
	    cacheParameterValues();

	    double passLikelihood = singlePassOptimization(pass, initialSearchSettingsForEachPass[pass]);

	    // update
	    if (pass == 0) {
		finalLikelihood = passLikelihood;

		if (Constants.WARNLEVEL > 1) { 
		    System.out.println ("Updating with likelihood " + passLikelihood + ".");
		}
	    }
	    else {
		if (passLikelihood > finalLikelihood) {
		    finalLikelihood = passLikelihood;
		    // no need for restore

		    if (Constants.WARNLEVEL > 1) { 
			System.out.println ("Updating with likelihood " + passLikelihood + ".");
		    }
		}
		else {
		    // need to restore
		    restoreParameterValues();

		    if (Constants.WARNLEVEL > 1) { 
			System.out.println ("Not updating.");
		    }
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
	    double originalSetting = parameter.getValue();
	   
	    // update
	    parameter.setValue(x);
	    // update associated model values associated with a single parameter
	    // e.g., for LengthParameter objects that belong to a ParameterConstraintSet
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

    public class NodeEqualityTestException extends RuntimeException {
	public NodeEqualityTestException (String message) {
	    super(message);
	}
    }

}




// Shoot - need to name everything.
// Since might need to fix gene genealogy branch lengths.
// Also neater this way.
// This means that all trees must be named.
