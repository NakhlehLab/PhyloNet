package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core.TwoStageState;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.*;

/**
 * Created by dw20 on 5/11/17.
 */
public abstract class NetworkFromGTTTwoStageState extends NetworkFromGTT implements TwoStageState {

    // gene tree correspondences
    protected List<Tuple<MutableTuple<Tree, Double>, Set<Integer>>> _treeCorrespondences;
    // gene tree info for pseudo likelihood
    protected List<String> _allTriplets;
    protected List<List<Tuple<double[][],Double>>> _tripletFrequencies;
    protected GTTLikelihood _pseudoCalculation;

    public NetworkFromGTTTwoStageState(List<List<MutableTuple<Tree, Double>>> inputTrees,
                                       Map<String, List<String>> taxonMap, long seed, int parallel) {
        super(inputTrees, taxonMap, seed, parallel);
    }

    /**
     * Calculate the log likelihood(this).
     * @return
     */
    public double calculateLikelihood() {
        _calculation.setParallel(_numThreads);
        double logL = _calculation.computeProbability(_speciesNet, _geneTrees, _taxonMap, _treeCorrespondences);
        return logL;
    }

    /**
     * Calculate the log pseudo likelihood(this).
     * @return
     */
    public double calculatePseudoLikelihood() {
        _pseudoCalculation.setParallel(_numThreads);
        double logL = _pseudoCalculation.computeProbability(_speciesNet, _allTriplets, _taxonMap, _tripletFrequencies);
        return logL;
    }

    /**
     * summarize gene trees
     */
    protected void processGeneTrees() {
        // get gene trees and parameters
        _geneTrees = new ArrayList<>();
        _gtsForStartingNet = new ArrayList<>();
        _treeCorrespondences = new ArrayList<>();
        _allTriplets  = new ArrayList<>();
        _tripletFrequencies = new ArrayList<>();
        // -- summarize gene trees
        _pseudoCalculation.summarizeData(
                _inputWeightedGTs, _allele2SpeciesMap, _gtsForStartingNet, _allTriplets, _tripletFrequencies);
        _gtsForStartingNet.clear();
        _calculation.summarizeData(
                _inputWeightedGTs, _allele2SpeciesMap, _gtsForStartingNet, _geneTrees, _treeCorrespondences);
        // print trees & weights
        if(_printDetails) {
            reportGeneTrees();
        }
    }

}
