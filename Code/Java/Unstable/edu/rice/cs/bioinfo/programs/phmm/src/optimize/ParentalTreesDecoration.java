/**
 * Container class for decoration/annotation of
 * parental tree object (represented using Network<Double> objects)
 * and their branch length parameters.
 *
 * For convenience.
 */

package optimize;

import java.util.Hashtable;
import java.util.Vector;
import java.util.StringTokenizer;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.File;
import util.Constants;
import runHmm.runHmm;
import be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.library.programming.BijectiveHashtable;
import edu.rice.cs.bioinfo.library.programming.BidirectionalMultimap;
import edu.rice.cs.bioinfo.library.programming.Tuple;


public class ParentalTreesDecoration {
    public static final String PARENTAL_NODE_LABEL_DELIMITER = ",";

    // indirection between parental tree name and parental tree objects themselves
    protected BijectiveHashtable<String,Network<Double>> parentalTreeNameMap;

    // for convenience, map from a parental tree node back to the parental tree that contains it
    // no similar method in PhyloNet's network structure
    // although one exists in PhyloNet's tree structure 
    protected Hashtable<NetNode<Double>,Network<Double>> parentalNodeTreeMap;

    // bijective map between parental tree node objects and their names
    // naming convention:
    // "<tree name>,<node name>"
    protected BijectiveHashtable<NetNode<Double>,String> parentalNodeLabelMap;

    // Use a single BidirectionalMultimap to capture many-to-many 
    // relationship between
    // length parameters and edge ids:
    // lid <-> eid
    protected BidirectionalMultimap<ParentalBranchLengthParameter,Tuple<String,String>> lpEidMap;
    // reverse lookup from lid name to length-parameter references
    protected Hashtable<String,ParentalBranchLengthParameter> lidLpMap;

    // One more layer of indirection.
    // Group length parameters into length-parameter-constraint-sets.
    // Each length-parameter-constraint-set has a total length constraint.
    // If length-parameter belongs to a length-parameter-constaint-set, then
    // enable fixed-total-length optimization for the set.
    protected BidirectionalMultimap<ParameterConstraintSet,ParentalBranchLengthParameter> setLpMap;
    // By convention, first length-parameter in a length-parameter-constraint-set 
    // has relative weight 1.0 and isn't optimized.
    // Beware degree-of-freedoms == one less than number of length-parameters in
    // a length-parameter-constraint-set.
    //
    // Store this info in sid->lid map.
    protected Hashtable<ParameterConstraintSet,ParentalBranchLengthParameter> setFirstLpMap;

    // shoot - since ParentalBranchLengthParameter object construction occurs in here
    protected runHmm runHmmObject;

    protected CalculationCache calculationCache;

    public ParentalTreesDecoration(BijectiveHashtable<String,Network<Double>> inParentalTreeNameMap,
				   String inputParentalBranchLengthParameterToEdgeMapFilename,
				   String inputParentalBranchLengthParameterSetConstraintsFilename,
				   runHmm inRunHmmObject,
				   CalculationCache inCalculationCache
				   ) {
	this.parentalTreeNameMap = inParentalTreeNameMap;
	this.runHmmObject = inRunHmmObject;
	this.calculationCache = inCalculationCache;

	createParentalNodeLabelMap();
	createParentalNodeTreeMap();
	
	// no model update
	// this creates LengthParameter objects representing parental tree branch length parameters
	parseInputParentalBranchLengthParameterToEdgeMapFile(inputParentalBranchLengthParameterToEdgeMapFilename);
	parseParentalBranchLengthParameterConstraintSetsFile(inputParentalBranchLengthParameterSetConstraintsFilename);
    }

    /**
     * Get parental branch length parameters created during class construction.
     * WARNING - allocated new Vector for each call.
     */
    public Vector<ParentalBranchLengthParameter> getParentalBranchLengthParameters () {
	return (new Vector<ParentalBranchLengthParameter>(lpEidMap.keys()));
    }

    protected void createParentalNodeTreeMap () {
	parentalNodeTreeMap = new Hashtable<NetNode<Double>,Network<Double>>();
	for (Network<Double> parentalTree : parentalTreeNameMap.values()) {
	    for (NetNode<Double> node : parentalTree.dfs()) {
		parentalNodeTreeMap.put(node, parentalTree);
	    }
	}
    }

    /**
     * Create map between parental nodes and their labels,
     * Also create a map between parental node and the tree that contains it.
     */
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
    protected void parseParentalBranchLengthParameterConstraintSetsFile (String filename) {
	if ((filename == null) ||
	    (filename.trim().length() <= 0)) {
	    throw (new RuntimeException ("ERROR: invalid filename in parseParameterConstraintSetsFile()."));
	}

	if (!(new File(filename)).isFile()) {
	    throw (new RuntimeException ("ERROR: filename " + filename + " does not exist."));
	}

	setLpMap = new BidirectionalMultimap<ParameterConstraintSet,ParentalBranchLengthParameter>();
	setFirstLpMap = new Hashtable<ParameterConstraintSet,ParentalBranchLengthParameter>();

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
		ParameterConstraintSet slp = new ParameterConstraintSet(sid, value);
		
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
		    ParentalBranchLengthParameter lp = lidLpMap.get(lid);
		    // duplicate lids fine
		    // all goes into a HashSet anyways
		    setLpMap.put(slp, lp);

		    // Keep track of first lid in set.
		    if (!setFirstLpMap.containsKey(slp)) {
			setFirstLpMap.put(slp, lp);
			// Weight of first lid is always 1.0.
			// no update during parse - delay until just right after everything parsed
			lp.setValue(1.0, true, true, false);
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
    protected void parseInputParentalBranchLengthParameterToEdgeMapFile (String filename) {
	if ((filename == null) ||
	    (filename.trim().length() <= 0)) {
	    throw (new RuntimeException ("ERROR: invalid filename in parseInputParentalBranchLengthParameterToEdgeMapFile()."));
	}

	if (!(new File(filename)).isFile()) {
	    throw (new RuntimeException ("ERROR: filename " + filename + " does not exist."));
	}

	lpEidMap = new BidirectionalMultimap<ParentalBranchLengthParameter,Tuple<String,String>>();
	lidLpMap = new Hashtable<String,ParentalBranchLengthParameter>();

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
		ParentalBranchLengthParameter lp = new ParentalBranchLengthParameter(ParentalBranchLengthParameter.class.getName() + HiddenState.HIDDEN_STATE_NAME_DELIMITER + lid, 
										     value,
										     runHmmObject,
										     // meh
										     this,
										     calculationCache,
										     true,
										     true,
										     // no update during parse
										     false
										     );
		
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
	    throw (new RuntimeException("ERROR: invalid lpEidMap in parseInputParentalBranchLengthParameterToEdgeMapFile."));
	}

	if (Constants.WARNLEVEL > 4) { System.out.println ("lpEidMap at end of parseInputParentalBranchLengthParameterToEdgeMapFile: \n" + lpEidMap.toString() + "\n"); }
    }

    // parentalTreeClasses.keySet()

    protected boolean checkLpEidMap () {
	// also make sure that every edge in the tree is contained in lpEidMap
	// 
	// actually - not necessary
	// excluded edges -> not optimized
	//
	// Just warn.
	for (Network<Double> parentalTree : parentalTreeNameMap.values()) {
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

    /**
     * Argh - expose this utility method.
     */
    public boolean checkEid (Tuple<String,String> eid) {
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

    public void updateBranchesBasedOnParentalBranchLengthParameters (ParentalBranchLengthParameter inputParentalBranchLengthParameter) {
	if (setLpMap.containsValue(inputParentalBranchLengthParameter)) {
	    // bleh - should really package up these indirection layers into a single class
	    ParameterConstraintSet slp = getParameterConstraintSet(inputParentalBranchLengthParameter);
	    for (ParentalBranchLengthParameter member : setLpMap.get(slp)) {
		updateBranchesBasedOnParentalBranchLengthParametersHelper(member);
	    }
	}
	else { 
	    // not a member of a constraint-set
	    updateBranchesBasedOnParentalBranchLengthParametersHelper(inputParentalBranchLengthParameter);
	}

    }

    /**
     * Update all network branch lengths associated with a length-parameter.
     * If length-parameter belongs to a length-parameter-constraint-set,
     * need to maintain total length in set.
     * 
     * WARNING: only updates inputLengthParameter!
     */
    protected void updateBranchesBasedOnParentalBranchLengthParametersHelper (ParentalBranchLengthParameter inputParentalBranchLengthParameter) {
	// if length-parameter is a member of a length-parameter-constraint-set
	// need to "atomically" update all members in the set

	//if (Constants.WARNLEVEL > 4) { System.out.println ("updateBranchesBasedOnParentalBranchLengthParametersHelper() for length-parameter with id " + inputParentalBranchLengthParameter.getName() + "."); }

	// Get relevant branches.
	// Weird? No branches returned from map?
	for (Tuple<String,String> eid : lpEidMap.get(inputParentalBranchLengthParameter)) {

	    //if (Constants.WARNLEVEL > 4) { System.out.println ("Processing eid: " + eid.toString()); }

	    // Update relevant branch.
	    if (!checkEid(eid)) {
		throw (new RuntimeException("ERROR: unknown node labels in eid " + eid + " in updateNetworkBranchLengths()."));
	    } 

	    NetNode<Double> child = parentalNodeLabelMap.rget(eid.Item1);
	    NetNode<Double> parent = parentalNodeLabelMap.rget(eid.Item2);
	    double length = 0.0;
	    // simple addition, nothing fancy
	    for (ParentalBranchLengthParameter lp : lpEidMap.rget(eid)) {
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

	    // invalidate cache for the associated tree
	    // paranoid
	    calculationCache.cacheProbabilityOfGeneGenealogyInParentalTree.clear(parentalNodeTreeMap.get(child));
	    calculationCache.cacheProbabilityOfGeneGenealogyInParentalTree.clear(parentalNodeTreeMap.get(parent));

	    //System.out.println ("setParentDistance: " + child.getName() + " " + parent.getName() + " " + length);
	}

	//if (Constants.WARNLEVEL > 4) { System.out.println ("updateBranchesBasedOnParentalBranchLengthParametersHelper() for length-parameter with id " + inputParentalBranchLengthParameter.getName() + " DONE."); }
    }

    /**
     * Compute the current value for a member of a 
     * length-parameter-constraint-set.
     * Would be more elegant to do some sort of caching.
     */
    protected double computeValueForMemberOfParameterConstraintSet (ParentalBranchLengthParameter lp) {
	// paranoid
	if (!setLpMap.containsValue(lp)) {
	    throw (new RuntimeException ("ERROR: computeValueForMemberOfParameterConstraintSet() called with length-parameter that is not a member of a length-parameter-constraint-set."));
	}

	// sum up weights from members of length-parameter-constraint-set
	// bleh - should really package up these indirection layers into a single class
	ParameterConstraintSet slp = getParameterConstraintSet(lp);
	double totalWeight = 0.0;
	for (ParentalBranchLengthParameter members : setLpMap.get(slp)) {
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
     * Convenience function.
     * Check to see if a length-parameter is the first in a
     * length-parameter-constraint-set.
     */
    public boolean checkFirstInParameterConstraintSet (ParentalBranchLengthParameter lp) {
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

    /**
     * Get length-parameter-constraint-set for a length-parameter.
     * Enforces injective property of setLpMap via RuntimeExceptions.
     * Returns null if length-parameter doesn't belong to a length-parameter-constraint-set.
     *
     * Argh - expose this utility function.
     */
    public ParameterConstraintSet getParameterConstraintSet (ParentalBranchLengthParameter lp) {
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

    public class NodeEqualityTestException extends RuntimeException {
	public NodeEqualityTestException (String message) {
	    super(message);
	}
    }


}


										     // parentalNodeLabelMap,
										     // lpEidMap,
										     // setLpMap,
										     // setFirstLpMap,
