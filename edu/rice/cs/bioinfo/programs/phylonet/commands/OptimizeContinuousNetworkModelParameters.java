/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * kl23
 * Heuristically run Brent's optimization method for unimodal continuous
 * optimization of continuous network model parameters.
 * Iteratively optimizes each parameter until convergence.
 *
 * Usage: OptimizeContinuousNetworkModelParameters <species network name> (<comma-separated list of gene tree topology names>) <optional arguments> <optional output filename>
 * where <optional arguments> is any of:
 * -a <allele-taxon map>
 * -p to enable detailed verbose debugging output
 * -d to disable branch length optimization
 * -e to disable reticulate edge probability optimization
 *
 * Totally change input format for branch length optimization.
 * Input file (-l option) specifies mapping from
 * a length parameter to the list of edges that include that 
 * length parameter.
 * For example, branch b_1 = l_1 + l_2,
 * branch b_2 = l_2, and branch b_3 = l_2.
 * Allows multiple branches to share length parameters.
 * Also, if a branch doesn't have an associated length parameter,
 * then its input extended-newick branch length is never changed.
 * 
 * Requires unique ID for each branch. Force input network
 * to have unique names for all nodes.
 *
 * Also allow fixed totals between sets of length parameters.
 *
 * WARNING - network branch lengths matter too!
 *
 * Completely convert to GeneTreeProbabilityYF, which uses
 * Yufeng Wu's ancestral configuration approach to speed
 * up calculation vs. MUL-tree approach of Yu et al.
 */

package edu.rice.cs.bioinfo.programs.phylonet.commands;

//import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.NetworkNonEmpty;
//import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
//import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.*;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF.CoalescePattern;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
//import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;


// kliu - additional imports
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;
import org.apache.commons.math3.optimization.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optimization.GoalType;

import java.io.*;
import java.util.*;

// command name for Nexus file
@CommandName("optimizecontinuousnetworkmodelparameters")
public class OptimizeContinuousNetworkModelParameters extends CommandBaseFileOut {
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
    public static final double DEFAULT_MAXIMUM_BRANCH_LENGTH = 1e1;

    // default weights for members of branch-length-constraint-sets
    public static final double DEFAULT_INITIAL_BRANCH_LENGTH_CONSTRAINT_SET_MEMBER_WEIGHT = 1.0;

    // only for randomization purposes
    public static final double MAXIMUM_BRANCH_LENGTH_FOR_RANDOMIZATION_PURPOSES = 2.0;

    // hmm... is GeneTreeProbabilityYF able to handle probability == 0 or 1?
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

    protected BrentOptimizer brentOptimizer;
    protected GeneTreeProbabilityYF gtpyf;
    protected HashMap<String,List<String>> _taxonMap = null;
    protected boolean  _printDetail = false;
    // Lookup nodes based on unique name.
    // nodes parameterized by unique identifer
    // need to store labels external to network
    // 
    // WARNING - need to maintain uniqueness of keys into this map by construction.
    // Since NetNode<T> doesn't have meaningful equals() or hashCode() functions!
    protected Hashtable<NetNode<CoalescePattern[]>,String> nodeLabelMap;
    // reverse above map
    protected Hashtable<String,NetNode<CoalescePattern[]>> labelNodeMap;
    protected Network<CoalescePattern[]> _speciesNetwork;
    //protected NetworkNonEmpty originalSpeciesNetwork;
    protected List<Tree> _geneTrees;
    // appears to be support for fractional weighting of each tree
    // not just unit counting
    protected List<Double> _geneTreeCounts;
    //protected List<Integer> _geneTreeCounts;

    // // length parameter -> edge ids
    // protected Hashtable<LengthParameter,Set<Tuple<String,String>>> lpEidMap;
    // // edge id -> length parameters
    // protected Hashtable<Tuple<String,String>,Set<LengthParameter>> eidLidMap;

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
    protected BidirectionalMultimap<LengthParameterConstraintSet,LengthParameter> setLpMap;
    // By convention, first length-parameter in a length-parameter-constraint-set 
    // has relative weight 1.0 and isn't optimized.
    // Beware degree-of-freedoms == one less than number of length-parameters in
    // a length-parameter-constraint-set.
    //
    // Store this info in sid->lid map.
    protected Hashtable<LengthParameterConstraintSet,LengthParameter> setFirstLpMap;

    // later, similarly between length parameters
    // and length-parameter-constraint-sets (with constraints on total length)

    // // Use integer labels in place of references to NetNode<T> objects.
    // // Since equals() and hashCode() functions aren't defined in NetNode<T>!
    // //
    // // for set branch length constraints
    // // Need three static hashes:
    // // 1. edge (node id,parent id) -> constraint set id i
    // protected Hashtable<Tuple<String,String>,String> edgeSidMap;
    // // 2. Constraint set id i -> vector of edges (node id,parent id) 
    // // (reverse of hash #1)
    // // by convention, always set first node's weight to 1.0.
    // protected Hashtable<String,Vector<Tuple<String,String>>> sidEdgesMap;
    // // 3. Constraint set id i -> total length b_l
    // protected Hashtable<String,Double> sidWeightMap;
    // // Need one dynamic hash. Just keep set-member weights separate from species network data structure.
    // // Use branch length information consistently.
    // protected Hashtable<Tuple<String,String>,Double> edgeWeightMap;

    protected boolean enableBranchLengthOptimizationFlag = true;
    protected boolean enableBranchProbabilityOptimizationFlag = true;

    // cache it
    protected RnNewickPrinter<CoalescePattern[]> rnNewickPrinter;

    protected double _bootstrap = 100;

    // public OptimizeContinuousNetworkModelParameters (SyntaxCommand motivatingCommand, 
    // 						     ArrayList<Parameter> params,
    // 						     Map<String,NetworkNonEmpty> sourceIdentToNetwork, 
    // 						     Proc3<String, Integer, Integer> errorDetected){
    //  super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);

    public OptimizeContinuousNetworkModelParameters (SyntaxCommand motivatingCommand, 
						     ArrayList<Parameter> params,
						     Map<String,NetworkNonEmpty>  sourceIdentToNetwork, 
						     Proc3<String, Integer, Integer> errorDetected,
						     RichNewickReader<Networks> rnReader){
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);

	// paranoid
	verifySearchSettings();

	gtpyf = new GeneTreeProbabilityYF();
	brentOptimizer = new BrentOptimizer(RELATIVE_ACCURACY, ABSOLUTE_ACCURACY);

	rnNewickPrinter = new RnNewickPrinter<CoalescePattern[]>();
    }

    @Override
    protected int getMinNumParams(){
        return 2;
    }

    @Override
    protected int getMaxNumParams(){
	// kliu was 5 before change in branch length parameterization input format.
        return 9;
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
     * crap
     * no available storage in BniNetNode
     * just use external Hashtable
     *
     * Only useNodeNameFlag if node names are guaranteed to be unique and non-empty!!!
     */
    protected void assignUniqueNodeLabels (boolean useNodeNameFlag) {
	// Strict! Force all nodes to have unique, non-empty names.
	// Since we need to be able to specify nodes/edges in input files.
	if (useNodeNameFlag && !checkUniqueNonemptyNodeNames()) {
	    throw (new RuntimeException ("ERROR: if using node names to uniquely label nodes in assignUniqueNodeLabels(boolean), all nodes in input species network must have unique, non-empty names."));
	}

	// initialize external maps
	nodeLabelMap = new Hashtable<NetNode<CoalescePattern[]>,String>();
	labelNodeMap = new Hashtable<String,NetNode<CoalescePattern[]>>();

	int id = 0;
	for(NetNode<CoalescePattern[]> node : _speciesNetwork.dfs()) {
	    if (useNodeNameFlag) {
		nodeLabelMap.put(node, node.getName());
		labelNodeMap.put(node.getName(), node);
	    }
	    else {
		nodeLabelMap.put(node, (new Integer(id)).toString());
		labelNodeMap.put((new Integer(id)).toString(), node);
	    }

	    id++;
	}
    }

    /**
     * annoying
     * only used by f() UnivariateFunction implementation classes below
     * kludged equals() function.
     * Warning - throws RuntimeException if node<->label map lookup fails.
     */
    protected boolean checkNodesEqual (NetNode<CoalescePattern[]> x, NetNode<CoalescePattern[]> y) {
	if (nodeLabelMap == null) {
	    throw (new NodeEqualityTestException("ERROR: node<->label map not initialized in checkNodesEqual()."));
	}

	if (!nodeLabelMap.containsKey(x)) {
	    throw (new NodeEqualityTestException("ERROR: node<->label map doesn't contain node " + x.getName() + "."));
	}

	if (!nodeLabelMap.containsKey(y)) {
	    throw (new NodeEqualityTestException("ERROR: node<->label map doesn't contain node " + y.getName() + "."));
	}

	return (nodeLabelMap.get(x).equals(nodeLabelMap.get(y)));
    }

    protected Network<CoalescePattern[]> convertSpeciesNetwork (NetworkNonEmpty inSpeciesNetwork) {
	NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
	Network<CoalescePattern[]> result = transformer.makeNetwork(inSpeciesNetwork);
	return (result);
    }

    protected Tuple<List<Tree>, List<Double>> convertGeneTrees (List<NetworkNonEmpty> inGeneTrees) {
        List<Tree> nbGeneTrees = new ArrayList<Tree>();
        List<Double> nbCounter = new ArrayList<Double>();
        for(NetworkNonEmpty geneTree : inGeneTrees){
            double prob = geneTree.TreeProbability.execute(new TreeProbabilityAlgo<Double, RuntimeException>() {
                @Override
                public Double forEmpty(TreeProbabilityEmpty empty) {
                    return 1.0;
                }

                @Override
                public Double forNonEmpty(TreeProbabilityNonEmpty nonEmpty) {
                    return Double.parseDouble(nonEmpty.ProbString);
                }
            });

            String phylonetGeneTree = NetworkTransformer.toENewickTree(geneTree);
            NewickReader nr = new NewickReader(new StringReader(phylonetGeneTree));
            STITree<Double> newtr = new STITree<Double>(true);
            if(_bootstrap<100){
                if(Trees.handleBootStrapInTree(newtr, _bootstrap)==-1){
                    throw new IllegalArgumentException("Input gene tree " + newtr + " have nodes that don't have bootstrap value");
                }

            }
            try
            {
                nr.readTree(newtr);
                    }
            catch(Exception e)
            {
                errorDetected.execute(e.getMessage(),
                        this._motivatingCommand.getLine(), this._motivatingCommand.getColumn());
            }
            boolean found = false;
            int index = 0;
            for(Tree tr: nbGeneTrees){
                if(Trees.haveSameRootedTopology(tr, newtr)){
                    found = true;
                    break;
                }
                index++;
            }
            if(found){
                nbCounter.set(index, nbCounter.get(index)+prob);
            }
            else{
                nbGeneTrees.add(newtr);
                nbCounter.add(prob);
            }
        }
	
	return (new Tuple<List<Tree>, List<Double>>(nbGeneTrees, nbCounter));
    }

    /**
     * Not a tuple - just a pair. Oh well.
     */
    // protected Tuple<List<Tree>, List<Integer>> convertGeneTrees (List<NetworkNonEmpty> inGeneTrees) {
    // 	// kliu - count duplicate gene trees
    //     List<Tree> geneTrees = new ArrayList<Tree>();
    //     List<Integer> counter = new ArrayList<Integer>();
    //     for(NetworkNonEmpty geneTree : inGeneTrees){
    // 	    // kliu - use this to convert between NetworkNonEmpty and STITree.
    //         String phylonetGeneTree = NetworkTransformer.toENewickTree(geneTree);
    // 	    // actually, simplest to just use this to construct 
    // 	    // STITree from a newick string, use this to pass into 
    // 	    // GeneTreeProbabilityYF.calculateGTDistribution()
    // 	    //
    // 	    // and use edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.ExNewickReader to 
    // 	    // read in an equivalent Network<T> object for a newick (or extended-newick) string
    //         NewickReader nr = new NewickReader(new StringReader(phylonetGeneTree));
    //         STITree<Integer> newtr = new STITree<Integer>(true);
    //         try
    //         {
    //             nr.readTree(newtr);
    //         }
    //         catch(Exception e)
    //         {
    //             errorDetected.execute(e.getMessage(),
    //                     this._motivatingCommand.getLine(), this._motivatingCommand.getColumn());
    //         }
    //         boolean found = false;
    //         int index = 0;
    //         for(Tree tr: geneTrees){
    //             if(Trees.haveSameRootedTopology(tr, newtr)){
    //                 found = true;
    //                 break;
    //             }
    //             index++;
    //         }
    //         if(found){
    //             counter.set(index, counter.get(index)+1);
    //         }
    //         else{
    //             geneTrees.add(newtr);
    //             counter.add(1);
    //         }
    //     }
	
    // 	return (new Tuple<List<Tree>, List<Integer>>(geneTrees, counter));
    // }

    @Override
    protected boolean checkParamsForCommand(){
        boolean noError = true;

        NetworkNonEmpty inSpeciesNetwork = this.assertAndGetNetwork(0);
        noError = noError && inSpeciesNetwork != null;

	//originalSpeciesNetwork = inSpeciesNetwork;

	// Just convert it up front. Prefer Cuong's data structure anyways.
	//
	// create a copy of the network represented in Cuong's data structure
	// using an input copy of the network represented in Matt's data structure
	// makes sense because Yun's code likely totally reliant on Cuong's data structure
        _speciesNetwork = convertSpeciesNetwork(inSpeciesNetwork);
	
	// follow Yun's lead - per discussion with Matt
	//
	// assign unique identifier to each node for equality test
	// ugly kludge - cuong's data structure needs to be updated with equals() and hashValue() method at a minimum
	// also compareTo()
	assignUniqueNodeLabels(true);

        ParameterIdentList geneTreeParam = this.assertParameterIdentList(1);
        noError = noError && geneTreeParam != null;
        List<NetworkNonEmpty> inGeneTrees = new LinkedList<NetworkNonEmpty>();
        for(String ident : geneTreeParam.Elements)
        {
            noError = noError && this.assertNetworkExists(ident, geneTreeParam.getLine(), geneTreeParam.getColumn());
            if(noError)
            {
                inGeneTrees.add(this.sourceIdentToNetwork.get(ident));
            }
        }

	// ditto with the gene trees
	Tuple<List<Tree>, List<Double>> result = convertGeneTrees(inGeneTrees);
	_geneTrees = result.Item1;
	_geneTreeCounts = result.Item2;
	
        //ParamExtractorAllelMap aParam = new ParamExtractorAllelMap("a", this.params, this.errorDetected);
        ParamExtractorAllelListMap aParam = new ParamExtractorAllelListMap("a", this.params, this.errorDetected);
        if(aParam.ContainsSwitch){
            noError = noError && aParam.IsValidMap;
            if(aParam.IsValidMap){
                _taxonMap = aParam.ValueMap;
            }
        }

        ParamExtractor pParam = new ParamExtractor("p", this.params, this.errorDetected);
        if(pParam.ContainsSwitch)
        {
            _printDetail = true;
        }

	gtpyf.setPrintDetails(_printDetail);

        ParamExtractor bParam = new ParamExtractor("b", this.params, this.errorDetected);
        if(bParam.ContainsSwitch)
        {
            if(bParam.PostSwitchParam != null)
            {
                try
                {
                    _bootstrap = Double.parseDouble(bParam.PostSwitchValue);
                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unrecognized bootstrap value " + bParam.PostSwitchValue, bParam.PostSwitchParam.getLine(), bParam.PostSwitchParam.getColumn());
                }
            }
            else
            {
                errorDetected.execute("Expected value after switch -b.", bParam.SwitchParam.getLine(), bParam.SwitchParam.getColumn());
            }
        }

	// disable branch length optimization
	ParamExtractor dParam = new ParamExtractor("d", this.params, this.errorDetected);
        if(dParam.ContainsSwitch)
        {
            enableBranchLengthOptimizationFlag = false;
        }

	// disable branch probability optimization
	ParamExtractor eParam = new ParamExtractor("e", this.params, this.errorDetected);
        if(eParam.ContainsSwitch)
        {
            enableBranchProbabilityOptimizationFlag = false;
        }

	// Input mapping between a length parameter and branches to which the length parameter contributes to.
	ParamExtractor lParam = new ParamExtractor("l", this.params, this.errorDetected);
        if(lParam.ContainsSwitch)
        {
	    // kliu - nice, this works just fine
	    // testing
	    //System.err.println("hoo! |" + lParam.PostSwitchValue + "|");

	    String inputLengthParameterToEdgeMapFile = lParam.PostSwitchValue;
	    // Process input map from edge id to list of constituent length parameters
	    parseInputLengthParameterToEdgeMapFile(inputLengthParameterToEdgeMapFile);
        }
	else {
	    noError = false;
	    // eh, stricter still
	    throw (new RuntimeException("ERROR: must specify l parameter to OptimizeContinuousNetworkModelParameters command."));
	}

	// WARNING - order matters! Below length-parameter-set-constraint file
	// must reference based on length-parameter ids in above file!
	//
	// Input file specifying constraints on total length of sets of length parameters.
	ParamExtractor cParam = new ParamExtractor("c", this.params, this.errorDetected);
        if(cParam.ContainsSwitch)
        {
	    String inputLengthParameterSetConstraintsFile = cParam.PostSwitchValue;
	    parseLengthParameterConstraintSetsFile(inputLengthParameterSetConstraintsFile);
        }
	else {
	    noError = false;
	    // eh, stricter still
	    throw (new RuntimeException("ERROR: must specify c parameter to OptimizeContinuousNetworkModelParameters command."));
	}

        noError = noError && checkForUnknownSwitches("p", "a", "d", "e", "l", "c", "b");
        checkAndSetOutFile(aParam, pParam, bParam, dParam, eParam, lParam, cParam);

        return  noError;
    }

    /**
     * Token format: <branch child id>,<branch parent id>
     */
    protected Tuple<String,String> parseEid (String filename, String s) {
	StringTokenizer st = new StringTokenizer(s, ",");
	if (st.countTokens() != 2) {
	    throw (new RuntimeException("ERROR: incorrectly formatted edge id in file " + filename + ": " + s));
	}

	String cid = st.nextToken();
	String pid = st.nextToken();
	return (new Tuple<String,String>(cid, pid));
    }

    /**
     * Input file format:
     * <length-parameter-constraint-set ID> <total length> <length parameter 1 unique ID> <length parameter 2 unique ID> ...
     *
     * A length parameter ID can appear AT MOST ONCE in this input file!
     * Avoid complex constraints - sets *MUST* be disjoint!
     */
    protected void parseLengthParameterConstraintSetsFile (String filename) {
	if ((filename == null) ||
	    (filename.trim().length() <= 0)) {
	    throw (new RuntimeException ("ERROR: invalid filename in parseLengthParameterConstraintSetsFile()."));
	}

	if (!(new File(filename)).isFile()) {
	    throw (new RuntimeException ("ERROR: filename " + filename + " does not exist."));
	}

	setLpMap = new BidirectionalMultimap<LengthParameterConstraintSet,LengthParameter>();
	setFirstLpMap = new Hashtable<LengthParameterConstraintSet,LengthParameter>();

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
		LengthParameterConstraintSet slp = new LengthParameter(sid, value);
		
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

	if (Constants.WARNLEVEL > 4) { System.out.println ("setLpMap at end of parseLengthParameterConstraintSetsFile(): \n" + setLpMap.toString() + "\n"); }
    }

    /**
     * Input file format:
     * <length parameter unique ID> <initial value> <branch 1 child id>,<branch 1 parent id> <branch 2 child id>,<branch 2 parent id> ...
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
	for(NetNode<CoalescePattern[]> node : _speciesNetwork.dfs()) {
	    for(NetNode<CoalescePattern[]> parent : node.getParents()) {
		Tuple<String,String> eid = new Tuple<String,String>(nodeLabelMap.get(node), nodeLabelMap.get(parent));

		// testing
		System.out.println ("Edge: |" + eid + "|");

		if (!lpEidMap.containsValue(eid)) {
		    if (Constants.WARNLEVEL > 1) { System.err.println ("WARNING: no map entry for edge with id " + eid + ". No length optimization will be performed for this edge."); }
		    //return (false);
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

    protected boolean checkUniqueNonemptyNodeNames () {
	HashSet<String> hs = new HashSet<String>();
	for(NetNode<CoalescePattern[]> node : _speciesNetwork.dfs()) {
	    if ((node.getName() == null) || (node.getName().trim().length() <= 0)) {
		return (false);
	    }

	    if (!hs.contains(node.getName())) {
		hs.add(node.getName());
	    }
	    else {
		return (false);
	    }
	}

	return (true);
    }


    /**
     * Convenience function.
     * Update all branch lengths based on current state of set of length parameters.
     */
    protected void updateNetworkBranchLengths () {
	// duplicated work here for each member of a length-parameter-constraint-set
	// oh well
	for (LengthParameter lp : lpEidMap.keys()) {
	    updateNetworkBranchLengths(lp);
	}
    }

    /**
     * Compute the current value for a member of a 
     * length-parameter-constraint-set.
     * Would be more elegant to do some sort of caching.
     */
    protected double computeValueForMemberOfLengthParameterConstraintSet (LengthParameter lp) {
	// paranoid
	if (!setLpMap.containsValue(lp)) {
	    throw (new RuntimeException ("ERROR: computeValueForMemberOfLengthParameterConstraintSet() called with length-parameter that is not a member of a length-parameter-constraint-set."));
	}

	// sum up weights from members of length-parameter-constraint-set
	LengthParameterConstraintSet slp = getLengthParameterConstraintSet(lp);
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
     * If inputLengthParameter belongs to a length-parameter-constraint-set,
     * updates all members of its length-parameter-constraint-set.
     */
    protected void updateNetworkBranchLengths (LengthParameter inputLengthParameter) {
	if (setLpMap.containsValue(inputLengthParameter)) {
	    LengthParameterConstraintSet slp = getLengthParameterConstraintSet(inputLengthParameter);
	    for (LengthParameter member : setLpMap.get(slp)) {
		updateNetworkBranchLengthsHelper(member);
	    }
	}
	else {
	    updateNetworkBranchLengthsHelper(inputLengthParameter);
	}
    }

    /**
     * Update all network branch lengths associated with a length-parameter.
     * If length-parameter belongs to a length-parameter-constraint-set,
     * need to maintain total length in set.
     * 
     * WARNING: only updates inputLengthParameter!
     */
    protected void updateNetworkBranchLengthsHelper (LengthParameter inputLengthParameter) {
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

	    NetNode<CoalescePattern[]> child = labelNodeMap.get(eid.Item1);
	    NetNode<CoalescePattern[]> parent = labelNodeMap.get(eid.Item2);
	    double length = 0.0;
	    // simple addition, nothing fancy
	    for (LengthParameter lp : lpEidMap.rget(eid)) {
		// if length-parameter belongs to length-parameter-constraint-set,
		// no storage - just re-compute value based on weighted ratios
		// each time
		if (setLpMap.containsValue(lp)) {
		    length += computeValueForMemberOfLengthParameterConstraintSet(lp);
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

    protected boolean checkEid (Tuple<String,String> eid) {
	if (!(labelNodeMap.containsKey(eid.Item1) && labelNodeMap.containsKey(eid.Item2))) {
	    return (false);
	}

	NetNode<CoalescePattern[]> c = labelNodeMap.get(eid.Item1);
	NetNode<CoalescePattern[]> p = labelNodeMap.get(eid.Item2);

	// argh
	for (NetNode<CoalescePattern[]> cp : c.getParents()) {
	    if (checkNodesEqual(p, cp)) {
		return (true);
	    }
	}
	
	return (false);
    }

    /**
     * Get length-parameter-constraint-set for a length-parameter.
     * Enforces injective property of setLpMap via RuntimeExceptions.
     * Returns null if length-parameter doesn't belong to a length-parameter-constraint-set.
     */
    protected LengthParameterConstraintSet getLengthParameterConstraintSet (LengthParameter lp) {
	if (!setLpMap.containsValue(lp)) {
	    return (null);
	}

	// paranoid
	// despite injective constraint on relation
	if (setLpMap.rget(lp).size() != 1) {
	    throw (new RuntimeException("ERROR: injective property violated in setLpMap for length-parameter with id " + lp.getName() + "."));
	}
	LengthParameterConstraintSet slp = setLpMap.rget(lp).iterator().next();
	
	return (slp);
    }

    /**
     * Convenience function.
     * Check to see if a length-parameter is the first in a
     * length-parameter-constraint-set.
     */
    protected boolean checkFirstInLengthParameterConstraintSet (LengthParameter lp) {
	if (setLpMap.containsValue(lp)) {
	    LengthParameterConstraintSet slp = getLengthParameterConstraintSet(lp);
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
     * Helper function.
     * Only update branch length if likelihood is strictly better.
     */
    protected double optimizeSingleLengthParameter (LengthParameter lp, double logLikelihood, int round) {
	DistanceParameterUnivariateFunction f = new DistanceParameterUnivariateFunction(lp);
	// l, x, u
	Tuple3<Double,Double,Double> searchInterval = getSearchInterval
	    (f, 
	     lp.getValue(), 
	     DEFAULT_MINIMUM_BRANCH_LENGTH, 
	     DEFAULT_MAXIMUM_BRANCH_LENGTH,
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
	    lp.setValue(upvp.getPoint());
	    updateNetworkBranchLengths(lp);
	    logLikelihood = brentOptimizedLogLikelihood;
	}
	else {
	    // no update - info instead
	    if (Constants.WARNLEVEL > 4) { System.out.println ("INFO: Round " + round + " optimized point and log likelihood for length parameter " + lp.getName() + " resulted in log likelihood " + brentOptimizedLogLikelihood + " which isn't better than current log likelihood " + logLikelihood + ". Not updating branch length nor best round log likelihood."); }
	}

	return (logLikelihood);
    }

    /**
     * Helper function.
     * Only update branch probability if likelihood is strictly better.
     */
    protected double optimizeSingleBranchProbability (NetNode<CoalescePattern[]> node, NetNode<CoalescePattern[]> parent, double logLikelihood, int round) {
	ProbabilityParameterUnivariateFunction f = new ProbabilityParameterUnivariateFunction(node, parent);
	// l, x, u
	Tuple3<Double,Double,Double> searchInterval = getSearchInterval
	    (f, 
	     node.getParentProbability(parent), 
	     DEFAULT_MINIMUM_PROBABILITY, 
	     DEFAULT_MAXIMUM_PROBABILITY,
	     "optimizeSingleBranchProbability: ");
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
	    // annoying - see comments in ProbabilityParameterUnivariateFunction.set()
	    // f already stored references to node->parent edge
	    f.set(upvp.getPoint());
	    logLikelihood = brentOptimizedLogLikelihood;
	}
	else {
	    // no update - info instead
	    if (Constants.WARNLEVEL > 4) { System.out.println ("INFO: Round " + round + " optimized point and log likelihood for probability on edge " + parent.getName() + "->" + node.getName() + " resulted in log likelihood " + brentOptimizedLogLikelihood + " which isn't better than current log likelihood " + logLikelihood + ". Not updating probability nor best round log likelihood."); }
	}

	return (logLikelihood);
    }

    /**
     * Convenience function.
     */
    protected String getNetworkString () {
	// man, how do I empty the buffer/reset the writer?
	// don't think there's support for this
	StringWriter sw = new StringWriter();
	rnNewickPrinter.print(_speciesNetwork, sw);
	String result = sw.toString();
	return (result);
    }

    protected void initializeDefaultDistanceAndUniformProbability () {
	// first set all length parameters
	if (enableBranchLengthOptimizationFlag) {
	    for (LengthParameter lp : lpEidMap.keys()) {
		if (checkFirstInLengthParameterConstraintSet(lp)) {
		    // canonical!
		    lp.setValue(1.0);
		}
		else {
		    lp.setValue(DEFAULT_INITIAL_BRANCH_LENGTH);
		}
	    }
	    updateNetworkBranchLengths();
	}

	// kliu - set initial parameter settings
	// per Yun, leaf branches don't factor into probability calculation
	// since they cannot have deep coalescences and hybridization can't occur between leaves
	//
	// Above assumes that only a single individual sampled per population.
	// 
	// adapted from Matt's code, since I'm too lazy to 
	// study Cuong's classes at length
	//
	// DFS from the root
	for(NetNode<CoalescePattern[]> node : _speciesNetwork.dfs()) {
	    // kliu - if multiple individuals sampled from a population,
	    // then leaf branches need optimization too.
	    for(NetNode<CoalescePattern[]> parent : node.getParents()) {		    
		// getParentProbability()/setParentProbability() for hybridization probabilities
		// default to equal hybridization probabilities
		// 
		// assume tree nodes initialized with parent probability 1
		if (enableBranchProbabilityOptimizationFlag && node.isNetworkNode()) {
		    // root has in-degree 0 - need above guard
		    node.setParentProbability(parent, 1.0 / node.getIndeg());
		}
		
		// testing
		// System.out.println ("parent support in initializeDefaultDistanceAndUniformProbability: " + node.getParentSupport(parent));
	    }
	}

    }

    /**
     * Can't let all branch lengths approach max.
     * Can cause model likelihood calculation to barf.
     * Cap the maximum randomized branch length.
     */
    protected void initializeRandom () {
	// first set all length parameters
	if (enableBranchLengthOptimizationFlag) {
	    for (LengthParameter lp : lpEidMap.keys()) {
		if (checkFirstInLengthParameterConstraintSet(lp)) {
		    // canonical!
		    lp.setValue(1.0);
		}
		else {
		    double randomDistance = Math.random() * (MAXIMUM_BRANCH_LENGTH_FOR_RANDOMIZATION_PURPOSES - DEFAULT_MINIMUM_BRANCH_LENGTH) + DEFAULT_MINIMUM_BRANCH_LENGTH;
		    lp.setValue(randomDistance);
		}
	    }
	    updateNetworkBranchLengths();
	}

	// kliu - set initial parameter settings
	// per Yun, leaf branches don't factor into probability calculation
	// since they cannot have deep coalescences and hybridization can't occur between leaves
	for(NetNode<CoalescePattern[]> node : _speciesNetwork.dfs()) {
	    // kliu - if multiple individuals sampled from a population,
	    // then leaf branches need optimization too.
	    double norm = 0.0;
	    for(NetNode<CoalescePattern[]> parent : node.getParents()) {
		// getParentProbability()/setParentProbability() for hybridization probabilities
		// default to equal hybridization probabilities
		// 
		// assume tree nodes initialized with parent probability 1
		if (enableBranchProbabilityOptimizationFlag && node.isNetworkNode()) {
		    // root has in-degree 0 - need above guard
		    double randomDouble = Math.random();
		    node.setParentProbability(parent, randomDouble);
		    norm += randomDouble;
		}
	    }

	    if (enableBranchProbabilityOptimizationFlag && node.isNetworkNode()) {
		// need to re-normalize all probabilities to produce a proper probability distribution.
		for(NetNode<CoalescePattern[]> parent : node.getParents()) {
		    node.setParentProbability(parent, node.getParentProbability(parent) / norm);
		}
	    }
	}	

    }

    /**
     * For diagnostic purposes.
     */
    protected void debugLengthParameters () {
	for (LengthParameter lp : lpEidMap.keys()) {
	    System.out.println (lp.toString());
	}
    }

    /**
     * Perform one optimization pass from a single starting point.
     * Either starting point uses current settings, random settings, or default initial rate/uniform probabilities.
     */
    protected double singlePassOptimization (int pass, InitialSearchSettings initialSearchSettings) {
	double inputLogLikelihood = computeGTProb(true);

	// testing
	//System.err.println ("Input log likelihood: |" + inputLogLikelihood + "|");
	//System.exit(1);

	if (Constants.WARNLEVEL > 1) { 
	    System.out.println ("Processing pass " + pass + " with initial search setting " + initialSearchSettings.toString() + ".");
	}

	if (Constants.WARNLEVEL > 1) { 
	    System.out.println ("Input network: |" + getNetworkString() + "|");
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
	    initializeDefaultDistanceAndUniformProbability();
	    break;
	}

	double initialLogLikelihood = computeGTProb();

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
	    System.out.println ("Initial network: |" + getNetworkString() + "|");
	    System.out.println ("Initial log likelihood: |" + initialLogLikelihood + "|");
	    // Also print out length-parameters.
	    System.out.println ("Initial length-parameters:");
	    debugLengthParameters();
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

	    // iterate through length parameters first
	    int lpCount = 0;
	    for (LengthParameter lp : lpEidMap.keys()) {
		// by convention, skip the first lengthParameter
		// in a length-parameter-constraint-set
		if (checkFirstInLengthParameterConstraintSet(lp)) {
		    continue;
		}

		if (Constants.WARNLEVEL > 1) { System.out.println ("Processing length parameter " + lp.getName() + " count " + lpCount + "."); }

		// single branch length optimization
		// does update appropriately
		if (enableBranchLengthOptimizationFlag) {
		    roundLogLikelihood = optimizeSingleLengthParameter(lp, roundLogLikelihood, round);
		}

		if (Constants.WARNLEVEL > 1) { System.out.println ("Processing length parameter " + lp.getName() + " count " + lpCount + " DONE."); }
		lpCount++;
	    }

	    // then, dfs to iterate through network node probabilities
	    double nodeCount = 0;
	    for(NetNode<CoalescePattern[]> node : _speciesNetwork.dfs()) {
	    // kliu - if multiple individuals sampled from a population,
	    // then leaf branches need optimization too.
		if (Constants.WARNLEVEL > 1) { System.out.println ("Processing node " + node.getName() + " node count " + nodeCount + "."); }

		Iterator<NetNode<CoalescePattern[]>> iter = node.getParents().iterator();
		while (iter.hasNext()) {
		    NetNode<CoalescePattern[]> parent = iter.next();
		    Tuple<String,String> eid = new Tuple<String,String>(nodeLabelMap.get(node), nodeLabelMap.get(parent));
		    
		    // if no optimization to be done, roundLogLikelihood
		    // always remains at initialLogLikelihood
		    
		    // Fix u searching in getSearchInterval!
		    // Too difficult to somehow only do halving on smaller incoming hybridization probability.
		    
		    // skip one of the parent probabilities
		    // free probability parameters == num parents - 1
		    //
		    // kliu - do it for all parents - since need incoming hybridization probabilities
		    // to sum to one
		    //
		    // getParentProbability()/setParentProbability() for hybridization probabilities
		    // default to equal hybridization probabilities
		    // 
		    // assume tree nodes initialized with parent probability 1
		    //
		    // skip last parent for this reason
		    if (enableBranchProbabilityOptimizationFlag && node.isNetworkNode() && iter.hasNext()) {
			roundLogLikelihood = optimizeSingleBranchProbability(node, parent, roundLogLikelihood, round);
		    }
		}

		if (Constants.WARNLEVEL > 1) { System.out.println ("Processing node " + node.getName() + " node count " + nodeCount + " DONE."); }
		nodeCount++;
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
	double finalLikelihood = computeGTProb(true);

	if (Constants.WARNLEVEL > 1) { 
	    System.out.println ("Pass network: |" + getNetworkString() + "|");
	    System.out.println ("Pass log likelihood: |" + finalLikelihood + "|");	    
	    // Also print out length-parameters.
	    System.out.println ("Pass length-parameters:");
	    debugLengthParameters();
	    System.out.println ("Processing pass " + pass + " with initial search setting " + initialSearchSettings.toString() + " DONE.");
	}

	return (finalLikelihood);
    }

    /**
     * f 
     * No clone constructors for any of the network data structures.
     * Just cache and restore all branch lengths and probabilities by hand.
     * For child_node <- parent_node edge, cache key is the string "label(child_node) label(parent_node)".
     *
     * Need to pass in two empty maps.
     */
    protected void cacheLengthParametersAndBranchProbabilities (Hashtable<LengthParameter,Double> lengthParameterCache, Hashtable<Tuple<String,String>,Double> probabilityCache) {
	// clear output caches just in case
	lengthParameterCache.clear();
	probabilityCache.clear();

	// paranoid
	for(NetNode<CoalescePattern[]> node : _speciesNetwork.dfs()) {
	    if (!nodeLabelMap.containsKey(node)) {
		throw (new RuntimeException("ERROR: unable to lookup label for node " + node.getName() + " in cacheBranchDistancesAndProbabilities()."));
	    }
	}

	// just do it for all nodes
	for(NetNode<CoalescePattern[]> node : _speciesNetwork.dfs()) {
	    for(NetNode<CoalescePattern[]> parent : node.getParents()) {
		Tuple<String,String> eid = new Tuple<String,String>(nodeLabelMap.get(node), nodeLabelMap.get(parent));
		if (node.isNetworkNode()) {
		    probabilityCache.put(eid, new Double(node.getParentProbability(parent)));
		}
	    }
	}
	
	// now for all length parameters
	for (LengthParameter lp : lpEidMap.keys()) {
	    lengthParameterCache.put(lp, new Double(lp.getValue()));
	}
    }

    protected void restoreLengthParametersAndBranchProbabilities (Hashtable<LengthParameter,Double> lengthParameterCache, Hashtable<Tuple<String,String>,Double> probabilityCache) {
	// paranoid
	for(NetNode<CoalescePattern[]> node : _speciesNetwork.dfs()) {
	    if (!nodeLabelMap.containsKey(node)) {
		throw (new RuntimeException("ERROR: unable to lookup label for node " + node.getName() + " in restoreBranchDistancesAndProbabilities()."));
	    }
	}

	// restore branch probabilities
	for(NetNode<CoalescePattern[]> node : _speciesNetwork.dfs()) {
	    for(NetNode<CoalescePattern[]> parent : node.getParents()) {
		Tuple<String,String> eid = new Tuple<String,String>(nodeLabelMap.get(node), nodeLabelMap.get(parent));
		if (node.isNetworkNode()) {
		    if (!probabilityCache.containsKey(eid)) {
			throw (new RuntimeException("ERROR: no cached probability for key " + eid.Item1 + " " + eid.Item2 + "."));
		    }
		    node.setParentProbability(parent, probabilityCache.get(eid).doubleValue());
		}
	    }
	}

	// restore all length parameters
	for (LengthParameter lp : lengthParameterCache.keySet()) {
	    lp.setValue(lengthParameterCache.get(lp).doubleValue());
	}

	// finally, update all branch lengths based on restored length parameter values
	updateNetworkBranchLengths();
    }

    @Override
    protected String produceResult () {
	Hashtable<LengthParameter,Double> lengthParameterCache = new Hashtable<LengthParameter,Double>();
	Hashtable<Tuple<String,String>,Double> probabilityCache = new Hashtable<Tuple<String,String>,Double>();

	// initial search settings matter quite a bit
	// try fewer rounds, but a couple of random restarts
	InitialSearchSettings[] initialSearchSettingsForEachPass = new InitialSearchSettings[] 
	    { InitialSearchSettings.CURRENT, 
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
	    cacheLengthParametersAndBranchProbabilities(lengthParameterCache, probabilityCache);

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
		    restoreLengthParametersAndBranchProbabilities(lengthParameterCache, probabilityCache);

		    if (Constants.WARNLEVEL > 1) { 
			System.out.println ("Not updating.");
		    }
		}
	    }
	}

	String resultString = 
	    "Final network: |" + getNetworkString() + "|\n" +
	    "Final log likelihood: |" + finalLikelihood + "|\n";

	// return output string later
	return (resultString);
    }

    /**
     * Need to satisfy four constraints:
     * 0. min, max make sense - see verifySearchSettings()
     * 1. x \in [l, u]
     * 2. x \in [min, max] and similarly for l and u
     * 3. f(x) is better than both f(l) and f(u)
     *
     * No log use here. Consistent with GeneTreeProbabilityYF library.
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

	Tuple3<Double,Double,Double> result = new Tuple3<Double,Double,Double>
	    (new Double(l), new Double(x), new Double(u));
	return (result);
    }

    protected double computeGTProb () {
	return (computeGTProb(false));
    }

    protected double computeGTProb (boolean debugFlag) {
	gtpyf.emptyState();

	// testing
	//NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
        //Network testSpeciesNetwork = transformer.makeNetwork(originalSpeciesNetwork);

	// need to use these!
	List<Tree> bGeneTrees = new ArrayList<Tree>();
        List<List<Integer>> nbTree2bTrees = new ArrayList<List<Integer>>();
        for(Tree nbgt: _geneTrees){
            List<Integer> bTrees = new ArrayList<Integer>();
            for(Tree bgt: Trees.getAllBinaryResolution(nbgt)){
                int index = 0;
                for(Tree exBgt: bGeneTrees){
                    if(Trees.haveSameRootedTopology(bgt,exBgt)){
                        break;
                    }
                    index++;
                }
                if(index==bGeneTrees.size()){
                    bGeneTrees.add(bgt);
                }
                bTrees.add(index);
            }
            nbTree2bTrees.add(bTrees);
        }

        List<Double> probList;
        //if(_multree){
        //    GeneTreeProbability gtp = new GeneTreeProbability();
        //    probList = gtp.calculateGTDistribution(speciesNetwork, bGeneTrees, _taxonMap, false);
        //}
        //else{
        //    GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        //    probList = gtp.calculateGTDistribution(speciesNetwork, bGeneTrees, _taxonMap, 0);
	//probList = gtpyf.calculateGTDistribution(_speciesNetwork, _geneTrees, _taxonMap, 0);

	// kliu - weird - child2parent doesn't look right - trying to index 5 but fails
	// should exist
	probList = gtpyf.calculateGTDistribution(_speciesNetwork, bGeneTrees, _taxonMap, 0);

	//GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
	//probList = gtp.calculateGTDistribution(_speciesNetwork, _geneTrees, _taxonMap, 0);


	    //}
        Iterator<Double> nbCounterIt = _geneTreeCounts.iterator();
        Iterator<List<Integer>> bGTIDs = nbTree2bTrees.iterator();
        double total = 0;
        for(Tree nbgt: _geneTrees){
            for(TNode node: nbgt.getNodes()){
                node.setParentDistance(TNode.NO_DISTANCE);
            }
            double maxProb = 0;
            for(int id: bGTIDs.next()){
                maxProb = Math.max(maxProb, probList.get(id));
            }
            double weight = nbCounterIt.next();
            total += Math.log(maxProb)*weight;
            if (debugFlag) {
		System.out.println("\n[x" + weight + "] " + nbgt.toString() + " : " + maxProb);
	    }
        }
	
	return (total);
        //result.append("\n" + "Total log probability: " + total);
    }

    /**
     * Use code from ComputeGTProb.
     * WARNING - returns log likelihood!
     */
    // protected double computeGTProb (boolean debugFlag) {
    // 	gtpyf.emptyState();

    // 	// calculation under model from Yu et al. 2012
    // 	// crap - this method requires Network<CoalescePattern[]>
    // 	// Yun uses Double to store hybridization probabilities during calculation
    //     Iterator<Double> probList = gtpyf.calculateGTDistribution(_speciesNetwork, _geneTrees, _taxonMap, 0).iterator();
    //     Iterator<Integer> countIter = _geneTreeCounts.iterator();
    //     double total = 0.0;
    //     for (Tree gt: _geneTrees){
    //         for (TNode node: gt.getNodes()){
    //             node.setParentDistance(TNode.NO_DISTANCE);
    //         }

    // 	    if (probList.hasNext() && countIter.hasNext()) {
    // 		double prob = probList.next();
    // 		int count = countIter.next();

    // 		// useful debugging info
    // 		if (debugFlag) {
    // 		    System.out.println("[x" + count + "] " + gt.toString() + " : " + prob);
    // 		}

    // 		total += Math.log(prob)*count;
    // 	    }
    // 	    else {
    // 		System.err.println ("ERROR: expecting another likelihood/count in computeGTProb(). Aborting and returning NaN.");
    // 		return (Double.NaN);
    // 	    }	    
    //     }

    //     return (total);
    // }

    public class DistanceParameterUnivariateFunction implements UnivariateFunction {
	protected LengthParameter lp;

	public DistanceParameterUnivariateFunction
	    (LengthParameter inlp) {
	    set(inlp);
	}

	public void set (LengthParameter inlp) {
	    lp = inlp;
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
	    double originalSetting = lp.getValue();
	   
	    // update
	    lp.setValue(x);
	    updateNetworkBranchLengths(lp);

	    // evaluate f(x)
	    double result = computeGTProb();

	    // restore original setting
	    lp.setValue(originalSetting);
	    updateNetworkBranchLengths(lp);

	    // above restore original setting op
	    // not really necessary if always accept new branch length even if no likelihood improvement
	    // since Brent's method is guaranteed to never make the likelihood worse
	    // meh, no update if no likelihood improvement
	    // -> above is necessary

	    return (result);
	}	
    }

    /**
     * Network edge probabilities need to be treated a little differently since they are probabilities \in [0,1], not rates.
     * For node n with parent p,
     * permit setting current p->n branch's probability to x \in [0,1], then for all other parents p'
     * allot 1-x probability weighted by current probability distribution on p'.
     * Best approach I can think of for setting probabilities.
     * Works like you'd expect for node with in-degree two, since there's effectively only a single
     * probability parameter. If in-degree > 2, then it gets a little
     * funkier since there's more than one probability parameter - can't change one without
     * re-scaling others to produce a valid probability distribution at node n.
     */
    public class ProbabilityParameterUnivariateFunction implements UnivariateFunction {
	// change a single continuous parameter on the edge parent->node
	protected NetNode<CoalescePattern[]> node;
	protected NetNode<CoalescePattern[]> parent;

	// cache parent node order
	protected Vector<NetNode<CoalescePattern[]>> parentOrder;

	public ProbabilityParameterUnivariateFunction
	    (NetNode<CoalescePattern[]> inNode, NetNode<CoalescePattern[]> inParent) {
	    set(inNode, inParent);
	}

	public void set (NetNode<CoalescePattern[]> inNode, NetNode<CoalescePattern[]> inParent) {
	    node = inNode;
	    parent = inParent;

	    cacheParentOrder();
	}

	protected void cacheParentOrder () {
	    parentOrder = new Vector<NetNode<CoalescePattern[]>>();
	    for (NetNode<CoalescePattern[]> p : node.getParents()) {
		parentOrder.add(p);
	    }
	}

	/**
	 * Cache parent probabilities.
	 * Output vector order corresponds to parentOrder.
	 */
	protected Vector<Double> cacheParentProbabilities () {
	    Vector<Double> probabilities = new Vector<Double>();
	    for (NetNode<CoalescePattern[]> p : parentOrder) {
		probabilities.add(new Double(node.getParentProbability(p)));
	    }
	    return (probabilities);
	}

	protected double getParentProbabilitiesExcludingCurrentParent () {
	    double total = 0.0;
	    for (NetNode<CoalescePattern[]> p : parentOrder) {
		// see kludge above
		// internally assign unique ids to all nodes
		// kludged BniNetNode.equals() method
		if (!checkNodesEqual(p, parent)) {
		    total += node.getParentProbability(p);
		}
	    }

	    return (total);
	}
	
	/**
	 * Somewhat annoying - see class comments.
	 * Another annoyance - caller needs access to this to repeat the operation if
	 * update rule satisfied. Expose method to caller. Oh well.
	 */
	public void set (double x) {
	    double remainingProbability = getParentProbabilitiesExcludingCurrentParent();

	    // re-scale by (1 - x) / remainingProbability 
	    // so that it modified distribution sums to 1
	    // while retaining relative proportion among probabilities other than x_old -> x
	    double factor = (1 - x) / remainingProbability;

	    for (NetNode<CoalescePattern[]> p : parentOrder) {
		// see kludge above
		// internally assign unique ids to all nodes
		// kludged BniNetNode.equals() method
		if (!checkNodesEqual(p, parent)) {
		    node.setParentProbability(p, node.getParentProbability(p) * factor);
		}
	    }

	    // now just set x_old -> x for current edge (node, parent)
	    node.setParentProbability(parent, x);
	}

	/**
	 * No safety check on input distribution.
	 * Entries in originalProbabilities correspond to order in parentOrder.
	 */
	protected void set (Vector<Double> originalProbabilities) {
	    // paranoid
	    // map is syntactically cleaner, with a little extra overhead
	    if (originalProbabilities.size() != parentOrder.size()) {
		throw (new UnivariateFunctionEvaluationException("ERROR: cached parent probabilities and cached parent order have different number of entries."));
	    }

	    for (int i = 0; i < parentOrder.size(); i++) {
		NetNode<CoalescePattern[]> p = parentOrder.get(i);
		node.setParentProbability(p, originalProbabilities.get(i).doubleValue());
	    }
	}

	/**
	 * Evaluate likelihood function f(x) for a particular x'.
	 * WARNING: f(x) is always a log likelihood, since computeGTProb() always returns a log likelihood!
	 */
	public double value (double x) {
	    //if (Constants.WARNLEVEL > 6) { System.out.println ("ProbabilityParameterUnivariateFunction.value() 1: |" + getNetworkString() + "|"); }

	    // cache original setting
	    Vector<Double> originalProbabilities = cacheParentProbabilities();

	    // set new settings
	    set(x);

	    //if (Constants.WARNLEVEL > 6) { System.out.println ("ProbabilityParameterUnivariateFunction.value() 2: |" + getNetworkString() + "|"); }

	    // evaluate f(x)
	    double result = computeGTProb();

	    //if (Constants.WARNLEVEL > 6) { System.out.println ("ProbabilityParameterUnivariateFunction.value() result: |" + result + "|"); }

	    // restore original setting
	    set(originalProbabilities);

	    //if (Constants.WARNLEVEL > 6) { System.out.println ("ProbabilityParameterUnivariateFunction.value() 3: |" + getNetworkString() + "|"); }

	    // above restore original setting op
	    // not really necessary if always accept new branch length even if no likelihood improvement
	    // since Brent's method is guaranteed to never make the likelihood worse
	    // meh, no update if no likelihood improvement
	    // -> above is necessary

	    return (result);
	}	
    }    

    /**
     * Just need some storage for constraint.
     * Immutable after creation.
     */
    public static class LengthParameterConstraintSet {
	protected String name;
	protected double value;

	public LengthParameterConstraintSet (String inName, double inValue) {
	    setName(inName);
	    setValue(inValue);
	}

	public String getName () {
	    return (name);
	}

	public double getValue () {
	    return (value);
	}

	// immutable once created
	protected void setName (String inName) {
	    name = inName;
	}

	protected void setValue (double inValue) {
	    value = inValue;
	}

	// LengthParameterConstraintSet objects are hashCode() and equals() equivalent to
	// their name Strings
	
	public String toString () {
	    return (name.toString() + ": " + Double.toString(value));
	}

	public int hashCode () {
	    return (name.hashCode());
	}

	// derp
	public boolean equals (Object obj) {
	    if (!(obj instanceof LengthParameterConstraintSet)) {
		return (false);
	    }
	    LengthParameterConstraintSet lpcs = (LengthParameterConstraintSet) obj;
	    return (this.getName().equals(lpcs.getName()));
	}
    }
    
    /**
     * Simple state container for length parameter.
     * Reduce indirection between IDs and objects a bit.
     *
     * Don't clone.
     */
    public static class LengthParameter extends LengthParameterConstraintSet {
	public LengthParameter (String inName, double inValue) {
	    super(inName, inValue);
	}

	public void setName (String inName) {
	    super.setName(inName);
	}

	public void setValue (double inValue) {
	    super.setValue(inValue);
	}

	// derp
	public boolean equals (Object obj) {
	    if (!(obj instanceof LengthParameter)) {
		return (false);
	    }
	    LengthParameter lp = (LengthParameter) obj;
	    return (this.getName().equals(lp.getName()));
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












    // /**
    //  * Convenience function.
    //  * Update all branch lengths that are subject to set branch length constraints.
    //  */
    // protected void updateAllSetBranchLengths () {
    // 	Enumeration<String> sids = sidEdgesMap.keys();
    // 	while (sids.hasMoreElements()) {
    // 	    String sid = sids.nextElement();
    // 	    updateSetBranchLengths(sid);
    // 	}
    // }

    // protected void updateSetBranchLengths (String sid) {
    // 	if ((sid == null) || (sid.equals(""))) {
    // 	    throw (new RuntimeException("ERROR: empty sid in updateSetBranchLengths."));
    // 	}

    // 	if (!sidEdgesMap.containsKey(sid) || !sidWeightMap.containsKey(sid)) {
    // 	    throw (new RuntimeException("ERROR: unable to locate hash entries for sid " + sid + " in updateSetBranchLengths."));
    // 	}

    // 	// Get total weight.
    // 	double totalWeight = 0.0;
    // 	for (Tuple<String,String> eid : sidEdgesMap.get(sid)) {
    // 	    totalWeight += edgeWeightMap.get(eid).doubleValue();
    // 	}
    // 	// solve
    // 	double frac = sidWeightMap.get(sid).doubleValue() / totalWeight;

    // 	// now update all members
    // 	for (Tuple<String,String> eid : sidEdgesMap.get(sid)) {
    // 	    NetNode<CoalescePattern[]> node = labelNodeMap.get(eid.Item1);
    // 	    NetNode<CoalescePattern[]> parent = labelNodeMap.get(eid.Item2);
    // 	    node.setParentDistance(parent, edgeWeightMap.get(eid).doubleValue() * frac);
    // 	}
    // }


    // /**
    //  * Need three static hashes:
    //  * 1. Node n -> constraint set id i
    //  * 2. Constraint set id i -> vector of nodes n (reverse of hash #1)
    //  * 3. Constraint set id i -> total length b_l
    //  */
    // protected void preprocessSetBranchLengthConstraints () {
    // 	// initialize hashes
    // 	edgeSidMap = new Hashtable<Tuple<String,String>,String>();
    // 	sidEdgesMap = new Hashtable<String,Vector<Tuple<String,String>>>();
    // 	sidWeightMap = new Hashtable<String,Double>();
    // 	edgeWeightMap = new Hashtable<Tuple<String,String>,Double>();

    // 	for(NetNode<CoalescePattern[]> node : _speciesNetwork.dfs()) {
    // 	    // only enable set branch length constraint for non-network nodes with a single parent
    // 	    // hack, since extended newick doesn't support arbitrary node labels
    // 	    // node name contains string "-CONSTRAINT-SET-<idstring>"
    // 	    if (!node.isNetworkNode() && (node.getParentNumber() == 1) && (node.getName().contains(IDENTIFIER_SET_BRANCH_LENGTH_CONSTRAINTS))) {
    // 		// Meh, can just retain IDENTIFIER_SET_BRANCH_LENGTH_CONSTRAINTS as prefix of set id.
    // 		String sid = node.getName().substring(node.getName().lastIndexOf(IDENTIFIER_SET_BRANCH_LENGTH_CONSTRAINTS));

    // 		// paranoid, despite above guard
    // 		for (NetNode<CoalescePattern[]> parent : node.getParents()) {
    // 		    Tuple<String,String> eid = new Tuple<String,String>(nodeLabelMap.get(node), nodeLabelMap.get(parent));
    // 		    double inputWeight = node.getParentDistance(parent);
		    
    // 		    // now create hash entries
    // 		    // dynamic hash first
    // 		    //
    // 		    // Use input weights.
    // 		    edgeWeightMap.put(eid, new Double(inputWeight));
		    
    // 		    edgeSidMap.put(eid, sid);
		    
    // 		    if (!sidEdgesMap.containsKey(sid)) {
    // 			sidEdgesMap.put(sid, new Vector<Tuple<String,String>>());
    // 		    }
    // 		    sidEdgesMap.get(sid).add(eid);

    // 		    // Add up all lengths in branch-length-constraint-set ->
    // 		    // maintain that total thereafter.
    // 		    if (!sidWeightMap.containsKey(sid)) {
    // 			sidWeightMap.put(sid, new Double(0.0));
    // 		    }
    // 		    sidWeightMap.put(sid, new Double(sidWeightMap.get(sid).doubleValue() + inputWeight));
    // 		}
    // 	    }
    // 	}

    // 	// normalize all weights in a branch-length-constraint-set so that first weight is 1.0 by default
    // 	Enumeration<String> sids = sidEdgesMap.keys();
    // 	while (sids.hasMoreElements()) {
    // 	    String sid = sids.nextElement();
    // 	    // paranoid
    // 	    if (sidEdgesMap.get(sid).size() > 0) {
    // 		int i = 0;
    // 		double norm = edgeWeightMap.get(sidEdgesMap.get(sid).firstElement()).doubleValue();
    // 		for (Tuple<String,String> eid : sidEdgesMap.get(sid)) {
    // 		    edgeWeightMap.put(eid, new Double(edgeWeightMap.get(eid).doubleValue() / norm));
    // 		    i++;
    // 		}
    // 	    }
    // 	}

    // 	// no need to do this
    // 	//updateAllSetBranchLengths();
    // }

	// Create hashes to enable set constrain optimization.
	// Must happen after assignUniqueNodeLabels() in above convertSpeciesNetwork() call.
	// Since we need unique node labels to do hash map lookups.
	// Unfortunately, NetNode<T> doesn't define reasonable equals() or hashCode() functions.
	//preprocessSetBranchLengthConstraints();

    // /**
    //  * Diagnostic function.
    //  */
    // protected void printBranchLengthConstraintSetWeights () {
    // 	System.out.println ("Branch-length-constraint-set weights: ");
    // 	Enumeration<Tuple<String,String>> eids = edgeWeightMap.keys();
    // 	int i = 0;
    // 	while (eids.hasMoreElements()) {
    // 	    Tuple<String,String> eid = eids.nextElement();
    // 	    System.out.println ("(" + labelNodeMap.get(eid.Item1).getName() + "," + labelNodeMap.get(eid.Item2).getName() + ") | (" + eid.Item1.toString() + "," + eid.Item2.toString() + ") | " + edgeWeightMap.get(eid).toString());
    // 	    i++;
    // 	}
    // 	System.out.println ("Total number of branch-length-constraint-set weights: " + i);
    // }


     // *
     // * For all nodes in same constraint set s that has total length b_l, 
     // * branch lengths become relative weights.
     // * (w_1 * x) + ... + (w_{|s|} * x) = b_l
     // * Need to solve for x and substitute during each f likelihood evaluation.
     // * Meh, don't bother keeping one of the weights to 1.
     // *
     // * Degree-of-freedoms is one less than the number of weights that we optimize,
     // * but that's ok. Saves a hash, somewhat increases the range of possible values
     // * for branches in s.
     // *
     // * Need to just create a new f function that does solve for x, update all branch lengths in s, and
     // * evaluate likelihood.
     // *
     // * Need three static hashes:
     // * 1. Node n -> constraint set id i
     // * 2. Constraint set id i -> vector of nodes n (reverse of hash #1)
     // * 3. Constraint set id i -> total length b_l
     // *
     // * Create all read-only during pre-processing step.

	// Also set all initial weights to default setting.
	// Enumeration<Tuple<String,String>> eids = edgeWeightMap.keys();
	// while (eids.hasMoreElements()) {
	//     Tuple<String,String> eid = eids.nextElement();
	//     // all equal weights - doesn't really matter
	//     // first member of branch-length-constraint-set gets weight 1.0 by convention
	//     if (sidEdgesMap.get(edgeSidMap.get(eid)).firstElement().equals(eid)) {
	// 	edgeWeightMap.put(eid, new Double(1.0));
	//     }
	//     else {
	// 	edgeWeightMap.put(eid, new Double(DEFAULT_INITIAL_BRANCH_LENGTH_CONSTRAINT_SET_MEMBER_WEIGHT));
	//     }
	// }

	// // now set all constraint-set branch lengths based on weights (all weights initialized to 1.0)
	// updateAllSetBranchLengths();

	// // Also set all initial weights to random setting.
	// Enumeration<Tuple<String,String>> eids = edgeWeightMap.keys();
	// while (eids.hasMoreElements()) {
	//     Tuple<String,String> eid = eids.nextElement();
	//     // all equal weights - doesn't really matter
	//     // first member of branch-length-constraint-set gets weight 1 by convention
	//     if (sidEdgesMap.get(edgeSidMap.get(eid)).firstElement().equals(eid)) {
	// 	edgeWeightMap.put(eid, new Double(1.0));
	//     }
	//     else {
	// 	double randomDistance = Math.random() * (MAXIMUM_BRANCH_LENGTH_FOR_RANDOMIZATION_PURPOSES - DEFAULT_MINIMUM_BRANCH_LENGTH) + DEFAULT_MINIMUM_BRANCH_LENGTH;
	// 	// weights and branch lengths are both in the interval (0,\infty)
	// 	edgeWeightMap.put(eid, new Double(randomDistance));
	//     }
	// }

	// // now set all constraint-set branch lengths based on weights (all weights initialized to 1.0)
	// updateAllSetBranchLengths();

		    //&&
		    // undefined support or support < 1.0 -> optimize branch length
		    //(Double.isNaN(node.getParentSupport(parent)) || (node.getParentSupport(parent) < 1.0)) 
		    //&& 
		    // not subject to branch-length-constraint-set constraints
		    // *or* not the first member of vector of branch-length-constraint-set
		    // 
		    // latter is due to degrees-of-freedom == sizeof(branch-length-constraint-set) - 1
		    //(!edgeSidMap.containsKey(eid) || !(sidEdgesMap.get(edgeSidMap.get(eid)).firstElement().equals(eid)))
