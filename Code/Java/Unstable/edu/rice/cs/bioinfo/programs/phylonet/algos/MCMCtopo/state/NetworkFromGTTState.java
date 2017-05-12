package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.*;

/**
 * Created by dw20 on 5/11/17.
 */
public abstract class NetworkFromGTTState extends NetworkFromGTT {

    // summarized gene trees
    protected List<Tree> _geneTrees;
    // original input gene trees
    protected List<List<MutableTuple<Tree,Double>>> originalGeneTrees;
    // input gene trees for processing (a copy of original input gene trees)
    protected List<List<MutableTuple<Tree,Double>>> inputWeightedGTs;
    // gene trees used for starting network
    protected List<MutableTuple<Tree, Double>> gtsForStartingNet;
    // gene tree correspondences
    protected List<Tuple<MutableTuple<Tree, Double>, Set<Integer>>> treeCorrespondences;

    public NetworkFromGTTState(List<List<MutableTuple<Tree, Double>>> inputTrees,
                               Map<String, List<String>> taxonMap, long seed, int parallel) {
        super(taxonMap, seed, parallel);
        // gene trees
        this.originalGeneTrees = inputTrees;
        try{
            this.inputWeightedGTs = new ArrayList<>();
            for(List<MutableTuple<Tree,Double>> treeList : inputTrees) {
                List<MutableTuple<Tree,Double>> wgtList = new ArrayList<MutableTuple<Tree,Double>>();
                for(MutableTuple<Tree,Double> mt : treeList) {
                    Tree newt = new STITree(mt.Item1.toNewick());
                    wgtList.add(new MutableTuple<Tree, Double>(newt, new Double(mt.Item2)));
                }
                this.inputWeightedGTs.add(wgtList);
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Calculate the log likelihood(this).
     * @return
     */
    public double calculateLikelihood() {
        _calculation.setParallel(_numThreads);
        double logL = _calculation.computeProbability(_speciesNet, _geneTrees, _taxonMap, treeCorrespondences);
        return logL;
    }

    /**
     * summarize gene trees
     */
    protected void processGeneTrees() {
        // get gene trees and parameters
        _geneTrees = new ArrayList<>();
        gtsForStartingNet = new ArrayList<>();
        treeCorrespondences = new ArrayList<>();
        // -- summarize gene trees
        _calculation.summarizeData(inputWeightedGTs, null, gtsForStartingNet, _geneTrees, treeCorrespondences);
        // print trees & weights
        if(_printDetails) {
            for(MutableTuple<Tree, Double> tuple : gtsForStartingNet){
                System.out.println(tuple.Item1.toNewick() + "  " + tuple.Item2.toString());
            }
        }
    }

    /**
     * gets starting network as MDC tree
     */
    protected Network getStartingNetwork(Network net) {
        if(net != null) {
            Utils.setBranchLengths(net);
            return net;
        }
        String start = Utils.getStartNetwork(gtsForStartingNet,
                _taxonMap, new HashSet<String>(), _speciesNet);
        if(_printDetails) {
            System.out.println("Get starting network as MDC tree: " + start);
        }
        return Networks.readNetwork(start);
    }

}
