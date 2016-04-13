package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core.NetworkPrior;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core.State;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.*;

/**
 * Created by wendingqiao on 11/3/14.
 */
public abstract class NetworkFromGTTMultiPerLocus<T> implements State {

    // species network
    protected Network _speciesNet;
    // summarized gene trees
    protected List<Tree> _geneTrees;
    // original input gene trees
    protected List<List<MutableTuple<Tree,Double>>> originalGeneTrees;
    // parameters for summarizing gene trees
    protected boolean printDetails = false;
    // input gene trees for processing (a copy of original input gene trees)
    protected List<List<MutableTuple<Tree,Double>>> inputWeightedGTs;
    // gene trees used for starting network
    protected List<MutableTuple<Tree, Double>> gtsForStartingNet;
    // gene tree correspondences
    protected List<Tuple<MutableTuple<Tree, Double>, Set<Integer>>> treeCorrespondences;
    // for likelihood calculation
    protected GTTLikelihood_MultiPerLocus calculation;
    // random
    protected long _seed;
    protected Random _random;
    // for prior
    private List<NetNode> startLeaves;
    // parallel
    private int numThreads;
    // taxon map
    private Map<String, List<String>> _taxonMap;

    /**
     * Constructor
     * @param start
     * @param inputTrees
     * @param seed
     */
    public NetworkFromGTTMultiPerLocus(Network start,
                                       List<List<MutableTuple<Tree, Double>>> inputTrees,
                                       Map<String, List<String>> taxonMap,
                                       long seed,
                                       int parallel) {
        this._taxonMap = taxonMap;
        this.numThreads = parallel;
        // gene trees
        this.originalGeneTrees = inputTrees;
        this.calculation = new GTTLikelihood_MultiPerLocus();
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
            processGeneTrees();
            // set starting network
            _speciesNet = getStartingNetwork(start);
            // random
            this._seed = seed;
            this._random = new Random(seed);
            this.startLeaves = IterableHelp.toList(_speciesNet.getLeaves());

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }


    /**
     * Calculate the log likelihood(this).
     * @return
     */
    public double calculateLikelihood() {
        calculation.setParallel(numThreads);
        double logL = calculation.computeProbability(_speciesNet, _geneTrees, _taxonMap, treeCorrespondences);
        return logL;
    }


    /**
     * Calculates prior value
     * @param k
     * @return
     */
    public double calculatePrior(double k) {
        return (new NetworkPrior(_speciesNet, k, startLeaves).calculate());
    }

    /**
     * Store state as a sample of GTT
     * @return
     */
    public List<Tuple<String,Double>> storeState(double posterior) {
        List<Tuple<String,Double>> list = new ArrayList<>();
        list.add(new Tuple(_speciesNet.toString(), posterior));
        return list;
    }

    /**
     * Presents network as a string
     * @return
     */
    public String toString() {
        return _speciesNet.toString();
    }

    /**
     * gets number of reticulation nodes
     * @return
     */
    public int numOfReticulation() {
        return _speciesNet.getReticulationCount();
    }

    public void setNetwork(Network net) {
        this._speciesNet = net;
    }

    public Network getNetwork() {
        return this._speciesNet;
    }

    /**
     * summarize gene trees
     */
    private void processGeneTrees() {
        // get gene trees and parameters
        _geneTrees = new ArrayList<>();
        gtsForStartingNet = new ArrayList<>();
        treeCorrespondences = new ArrayList<>();
        // -- summarize gene trees
        calculation = new GTTLikelihood_MultiPerLocus();
        calculation.summarizeData(inputWeightedGTs, null, gtsForStartingNet, _geneTrees, treeCorrespondences);
        // print trees & weights
        if(printDetails) {
            for(MutableTuple<Tree, Double> tuple : gtsForStartingNet){
                System.out.println(tuple.Item1.toNewick() + "  " + tuple.Item2.toString());
            }
        }
    }

    /**
     * gets starting network as MDC tree
     */
    private Network getStartingNetwork(Network net) {
        if(net != null) {
            setBranchLengths(net);
            return net;
        }
        String start = calculation.getStartNetwork(gtsForStartingNet,
                _taxonMap, new HashSet<String>(), _speciesNet);
        if(printDetails) System.out.println("Get starting network as MDC tree: " + start);
        return Networks.readNetwork(start);
    }

    private void setBranchLengths(Network net) {
        for(Object n1 : net.bfs()) {
            NetNode node1 = (NetNode)n1;
            for(Object n2 : node1.getParents()) {
                NetNode node2 = (NetNode)n2;
                if( Double.isNaN( node1.getParentDistance(node2) ) ) {
                    node1.setParentDistance(node2, 1.0);
                }
                if(node1.isNetworkNode() && Double.isNaN( node1.getParentProbability(node2) ) ) {
                    node1.setParentProbability(node2, 0.5);
                }
            }
        }
    }


}
