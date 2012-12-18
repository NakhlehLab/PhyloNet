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
 */

package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.Constants;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

// kliu - additional imports
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;
import org.apache.commons.math3.optimization.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optimization.GoalType;

import java.io.*;
import java.util.*;

@CommandName("optimizecontinuousnetworkmodelparameters")
public class OptimizeContinuousNetworkModelParameters extends CommandBaseFileOut {
    protected static final double RELATIVE_ACCURACY = 1e-12;
    protected static final double ABSOLUTE_ACCURACY = 1e-8;
    protected static final int BRENT_METHOD_SINGLE_ROUND_MAXIMUM_ITERATIONS = 100;
    protected static final int MAXIMUM_NUM_ROUNDS = 100;
    protected static final double MINIMUM_LOG_LIKELIHOOD_DELTA_FOR_CONVERGENCE = ABSOLUTE_ACCURACY;
    protected static final double MINIMUM_FACTOR_WIDTH_FOR_MINIMUM_MAXIMUM_INTERVAL = 4.0;

    // search defaults
    // don't allow zero branch lengths
    public static final double DEFAULT_MINIMUM_BRANCH_LENGTH = 1e-20;
    public static final double DEFAULT_INITIAL_BRANCH_LENGTH = 1e-4;
    public static final double DEFAULT_MAXIMUM_BRANCH_LENGTH = 1e20;

    // hmm... is GeneTreeProbability able to handle probability == 0 or 1?
    // yes, it handles this fine
    public static final double DEFAULT_MINIMUM_PROBABILITY = 0.0;
    // initial probability always initialized to uniform
    public static final double DEFAULT_MAXIMUM_PROBABILITY = 1.0;

    // to support fixed network probabilities
    // over-parameterization issue with optimization of model of Yu et al. 2012
    public static final boolean ENABLE_BRANCH_PROBABILITY_OPTIMIZATION_FLAG = true;

    protected BrentOptimizer brentOptimizer;
    protected GeneTreeProbability gtp;
    protected HashMap<String,String> _taxonMap = null;
    protected boolean  _printDetail = false;
    // nodes parameterized by unique identifer
    // need to store labels external to network
    protected Hashtable<NetNode<Double>,Integer> nodeLabelMap;
    protected Network<Double> _speciesNetwork;
    protected List<Tree> _geneTrees;
    protected List<Integer> _geneTreeCounts;

    // cache it
    protected RnNewickPrinter<Double> rnNewickPrinter;

    public OptimizeContinuousNetworkModelParameters (SyntaxCommand motivatingCommand, 
						     ArrayList<Parameter> params,
						     Map<String,NetworkNonEmpty> sourceIdentToNetwork, 
						     Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader){
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);

	// paranoid
	verifySearchSettings();

	gtp = new GeneTreeProbability();
	brentOptimizer = new BrentOptimizer(RELATIVE_ACCURACY, ABSOLUTE_ACCURACY);

	rnNewickPrinter = new RnNewickPrinter<Double>();
    }

    @Override
    protected int getMinNumParams(){
        return 2;
    }

    @Override
    protected int getMaxNumParams(){
        return 5;
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
     */
    protected void assignUniqueNodeLabels (Network<Double> n) {
	// initialize external map
	nodeLabelMap = new Hashtable<NetNode<Double>,Integer>();

	int id = 0;
	for(NetNode<Double> node : n.dfs()) {
	    nodeLabelMap.put(node, new Integer(id));
	    id++;
	}
    }

    /**
     * annoying
     * only used by f() UnivariateFunction implementation classes below
     * kludged equals() function.
     * Warning - throws RuntimeException if node<->label map lookup fails.
     */
    protected boolean checkNodesEqual (NetNode<Double> x, NetNode<Double> y) {
	if (nodeLabelMap == null) {
	    throw (new NodeEqualityTestException("ERROR: node<->label map not initialized in checkNodesEqual()."));
	}

	if (!nodeLabelMap.containsKey(x)) {
	    throw (new NodeEqualityTestException("ERROR: node<->label map doesn't contain node " + x.getName() + "."));
	}

	if (!nodeLabelMap.containsKey(y)) {
	    throw (new NodeEqualityTestException("ERROR: node<->label map doesn't contain node " + y.getName() + "."));
	}

	return (nodeLabelMap.get(x).intValue() == nodeLabelMap.get(y).intValue());
    }

    protected Network<Double> convertSpeciesNetwork (NetworkNonEmpty inSpeciesNetwork) {
	NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
	Network<Double> result = transformer.makeNetwork(inSpeciesNetwork);
	// follow Yun's lead - per discussion with Matt
	//
	// assign unique identifier to each node for equality test
	// ugly kludge - cuong's data structure needs to be updated with equals() and hashValue() method at a minimum
	// also compareTo()
	assignUniqueNodeLabels(result);
	return (result);
    }

    /**
     * Not a tuple - just a pair. Oh well.
     */
    protected Tuple<List<Tree>, List<Integer>> convertGeneTrees (List<NetworkNonEmpty> inGeneTrees) {
	// kliu - count duplicate gene trees
        List<Tree> geneTrees = new ArrayList<Tree>();
        List<Integer> counter = new ArrayList<Integer>();
        for(NetworkNonEmpty geneTree : inGeneTrees){
            String phylonetGeneTree = NetworkTransformer.toENewickTree(geneTree);
            NewickReader nr = new NewickReader(new StringReader(phylonetGeneTree));
            STITree<Integer> newtr = new STITree<Integer>(true);
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
            for(Tree tr: geneTrees){
                if(Trees.haveSameRootedTopology(tr, newtr)){
                    found = true;
                    break;
                }
                index++;
            }
            if(found){
                counter.set(index, counter.get(index)+1);
            }
            else{
                geneTrees.add(newtr);
                counter.add(1);
            }
        }
	
	return (new Tuple<List<Tree>, List<Integer>>(geneTrees, counter));
    }

    @Override
    protected boolean checkParamsForCommand(){
        boolean noError = true;

        NetworkNonEmpty inSpeciesNetwork = this.assertAndGetNetwork(0);
        noError = noError && inSpeciesNetwork != null;

	// Just convert it up front. Prefer Cuong's data structure anyways.
	//
	// create a copy of the network represented in Cuong's data structure
	// using an input copy of the network represented in Matt's data structure
	// makes sense because Yun's code likely totally reliant on Cuong's data structure
        _speciesNetwork = convertSpeciesNetwork(inSpeciesNetwork);

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
	Tuple<List<Tree>, List<Integer>> result = convertGeneTrees(inGeneTrees);
	_geneTrees = result.Item1;
	_geneTreeCounts = result.Item2;
	
        ParamExtractorAllelMap aParam = new ParamExtractorAllelMap("a", this.params, this.errorDetected);
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

        noError = noError && checkForUnknownSwitches("p", "a");
        checkAndSetOutFile(aParam, pParam);

        return  noError;
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
    protected double optimizeSingleBranchDistance (NetNode<Double> node, NetNode<Double> parent, double logLikelihood, int round) {
	DistanceParameterUnivariateFunction f = new DistanceParameterUnivariateFunction(node, parent);
	// l, x, u
	Tuple3<Double,Double,Double> searchInterval = getSearchInterval
	    (f, 
	     node.getParentDistance(parent), 
	     DEFAULT_MINIMUM_BRANCH_LENGTH, 
	     DEFAULT_MAXIMUM_BRANCH_LENGTH,
	     "optimizeSingleBranchDistance: ");
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
	    node.setParentDistance(parent, upvp.getPoint());
	    logLikelihood = brentOptimizedLogLikelihood;
	}
	else {
	    // no update - info instead
	    if (Constants.WARNLEVEL > 4) { System.out.println ("INFO: Round " + round + " optimized point and log likelihood for branch length on edge " + parent.getName() + "->" + node.getName() + " resulted in log likelihood " + brentOptimizedLogLikelihood + " which isn't better than current log likelihood " + logLikelihood + ". Not updating branch length nor best round log likelihood."); }
	}

	return (logLikelihood);
    }

    /**
     * Helper function.
     * Only update branch probability if likelihood is strictly better.
     */
    protected double optimizeSingleBranchProbability (NetNode<Double> node, NetNode<Double> parent, double logLikelihood, int round) {
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

    @Override
    protected String produceResult () {
	if (Constants.WARNLEVEL > 4) { 
	    System.out.println ("Input network: |" + getNetworkString() + "|");
	    System.out.println ("Input log likelihood: |" + computeGTProb() + "|");
	}
	
	// Cuong's network representation:
	// Nodes contain name, one datum, and links to children and parents (also "hidden"
	// members with distances, support, and probabilities - unnecessary).
	// Distances always associated with child->parent reference.
	// Man, children and parent lists are LinkedLists. Fine if node degrees are small.
	// Wrapper class just adds in root designation for a particular node
	// and utility functions.
	// Nice and simple.

	// kliu - set initial parameter settings
	// per Yun, leaf branches don't factor into probability calculation
	// since they cannot have deep coalescences and hybridization can't occur between leaves
	// 
	// adapted from Matt's code, since I'm too lazy to 
	// study Cuong's classes at length
	//
	// DFS from the root
	for(NetNode<Double> node : _speciesNetwork.dfs()) {
	    // kliu - do it for all parents - since need incoming hybridization probabilities
	    // to sum to one
	    if (!node.isLeaf()) {
		for(NetNode<Double> parent : node.getParents()) {
		    node.setParentDistance(parent, DEFAULT_INITIAL_BRANCH_LENGTH);
		    // getParentProbability()/setParentProbability() for hybridization probabilities
		    // default to equal hybridization probabilities
		    // 
		    // assume tree nodes initialized with parent probability 1
		    if (ENABLE_BRANCH_PROBABILITY_OPTIMIZATION_FLAG && node.isNetworkNode()) {
			// root has in-degree 0 - need above guard
			node.setParentProbability(parent, 1.0 / node.getIndeg());
		    }
		}
	    }
	}

	double initialLogLikelihood = computeGTProb();

	// looks good
	//
	// child.inDeg >= 2 -> only then parent probabilities printed
	// defaults to probability zero
	if (Constants.WARNLEVEL > 4) { 
	    System.out.println ("Initial network: |" + getNetworkString() + "|");
	    System.out.println ("Initial log likelihood: |" + initialLogLikelihood + "|");
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
	    for(NetNode<Double> node : _speciesNetwork.dfs()) {
		// only internal branches contribute to model likelihood
		if (!node.isLeaf()) {
		    Iterator<NetNode<Double>> iter = node.getParents().iterator();
		    while (iter.hasNext()) {
			NetNode<Double> parent = iter.next();
			// single branch length optimization
			// does update appropriately
			roundLogLikelihood = optimizeSingleBranchDistance(node, parent, roundLogLikelihood, round);

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
			if (ENABLE_BRANCH_PROBABILITY_OPTIMIZATION_FLAG && node.isNetworkNode() && iter.hasNext()) {
			    roundLogLikelihood = optimizeSingleBranchProbability(node, parent, roundLogLikelihood, round);
			}
		    }
		}
	    }

	    // end of round
	    // update likelihoods
	    // if we did the swap at the start of each iteration, wouldn't need the third variable - meh
	    prevLogLikelihood = currLogLikelihood;
	    currLogLikelihood = roundLogLikelihood;

	    round++;
	}

	// recompute for debugging output
	double finalLikelihood = computeGTProb(true);

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
     * No log use here. Consistent with GeneTreeProbability library.
     *
     * Output is (l', x', u') search interval that satisfies above four constraints.
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

	// constraint #3:
	// expand search interval endpoints by half/double until satisfied
	while (!(f.value(x) > f.value(l))) {
	    l /= 2.0;

	    if (min > l) {
		l = min;
		break;
	    }

	    if (Constants.WARNLEVEL > 4) { System.out.println (debugMessage + " l x u fx fl fu log likelihoods: |" + l + " " + x + " " + u + " " + f.value(l) + " " + f.value(x) + " " + f.value(u) + "|"); }
	}

	while (!(f.value(x) > f.value(u))) {
	    u *= 2.0;

	    if (u > max) {
		u = max;
		break;
	    }

	    if (Constants.WARNLEVEL > 4) { System.out.println (debugMessage + " l x u fx fl fu log likelihoods: |" + l + " " + x + " " + u + " " + f.value(l) + " " + f.value(x) + " " + f.value(u) + "|"); }
	}

	if (Constants.WARNLEVEL > 4) { System.out.println (debugMessage + " l x u fx fl fu log likelihoods after satisfying constraint #3 in getSearchInterval() : |" + l + " " + x + " " + u + " " + f.value(l) + " " + f.value(x) + " " + f.value(u) + "|"); }

	Tuple3<Double,Double,Double> result = new Tuple3<Double,Double,Double>
	    (new Double(l), new Double(x), new Double(u));
	return (result);
    }

    protected double computeGTProb () {
	return (computeGTProb(false));
    }

    /**
     * Use code from ComputeGTProb.
     * WARNING - returns log likelihood!
     */
    protected double computeGTProb (boolean debugFlag) {
	gtp.emptyState();

	// calculation under model from Yu et al. 2012
	// crap - this method requires Network<Double>
	// Yun uses Double to store hybridization probabilities during calculation
        Iterator<Double> probList = gtp.calculateGTDistribution(_speciesNetwork, _geneTrees, _taxonMap, _printDetail).iterator();
        Iterator<Integer> countIter = _geneTreeCounts.iterator();
        double total = 0.0;
        for (Tree gt: _geneTrees){
            for (TNode node: gt.getNodes()){
                node.setParentDistance(TNode.NO_DISTANCE);
            }

	    if (probList.hasNext() && countIter.hasNext()) {
		double prob = probList.next();
		int count = countIter.next();

		// useful debugging info
		if (debugFlag) {
		    System.out.println("[x" + count + "] " + gt.toString() + " : " + prob);
		}

		total += Math.log(prob)*count;
	    }
	    else {
		System.err.println ("ERROR: expecting another likelihood/count in computeGTProb(). Aborting and returning NaN.");
		return (Double.NaN);
	    }	    
        }

        return (total);
    }

    public class DistanceParameterUnivariateFunction implements UnivariateFunction {
	// change distance on the edge parent->node
	protected NetNode<Double> node;
	protected NetNode<Double> parent;

	public DistanceParameterUnivariateFunction
	    (NetNode<Double> inNode, NetNode<Double> inParent) {
	    set(inNode, inParent);
	}

	public void set (NetNode<Double> inNode, NetNode<Double> inParent) {
	    node = inNode;
	    parent = inParent;
	}

	/**
	 * Evaluate likelihood function f(x) for a particular x'.
	 * WARNING: f(x) is always a log likelihood, since computeGTProb() always returns a log likelihood!
	 */
	public double value (double x) {
	    // cache original setting, set new setting
	    double originalSetting = node.getParentDistance(parent);
	   
	    // update
	    node.setParentDistance(parent, x);

	    // evaluate f(x)
	    double result = computeGTProb();

	    // restore original setting
	    node.setParentDistance(parent, originalSetting);

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
	protected NetNode<Double> node;
	protected NetNode<Double> parent;

	// cache parent node order
	protected Vector<NetNode<Double>> parentOrder;

	public ProbabilityParameterUnivariateFunction
	    (NetNode<Double> inNode, NetNode<Double> inParent) {
	    set(inNode, inParent);
	}

	public void set (NetNode<Double> inNode, NetNode<Double> inParent) {
	    node = inNode;
	    parent = inParent;

	    cacheParentOrder();
	}

	protected void cacheParentOrder () {
	    parentOrder = new Vector<NetNode<Double>>();
	    for (NetNode<Double> p : node.getParents()) {
		parentOrder.add(p);
	    }
	}

	/**
	 * Cache parent probabilities.
	 * Output vector order corresponds to parentOrder.
	 */
	protected Vector<Double> cacheParentProbabilities () {
	    Vector<Double> probabilities = new Vector<Double>();
	    for (NetNode<Double> p : parentOrder) {
		probabilities.add(new Double(node.getParentProbability(p)));
	    }
	    return (probabilities);
	}

	protected double getParentProbabilitiesExcludingCurrentParent () {
	    double total = 0.0;
	    for (NetNode<Double> p : parentOrder) {
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

	    for (NetNode<Double> p : parentOrder) {
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
		NetNode<Double> p = parentOrder.get(i);
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
