package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core.NetworkPrior;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core.State;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.proposal.NetworkProposal;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.*;

/**
 * Created by dw20 on 5/11/17.
 */
public abstract class NetworkFromGTT<T> implements State {

    // species network
    protected Network _speciesNet;
    // for likelihood calculation
    protected GTTLikelihood _calculation;

    // parameters for summarizing gene trees
    protected boolean _printDetails = false;
    // random
    protected long _seed;
    protected Random _random;
    // for prior
    protected List<NetNode> _startLeaves;
    // parallel
    protected int _numThreads;
    // taxon map
    protected Map<String, List<String>> _taxonMap;
    // allele species map
    protected Map<String, String> _allele2SpeciesMap;
    // operations
    protected NetworkProposal _operation;

    // summarized gene trees
    protected List<Tree> _geneTrees;
    // original input gene trees
    protected List<List<MutableTuple<Tree,Double>>> _originalGeneTrees;
    // input gene trees for processing (a copy of original input gene trees)
    protected List<List<MutableTuple<Tree,Double>>> _inputWeightedGTs;
    // gene trees used for starting network
    protected List<MutableTuple<Tree, Double>> _gtsForStartingNet;

    /**
     * Constructor
     */
    public NetworkFromGTT(List<List<MutableTuple<Tree, Double>>> inputTrees,
                          Map<String, List<String>> taxonMap,
                          long seed,
                          int parallel) {
        this._taxonMap = taxonMap;
        if(_taxonMap != null) {
            this._allele2SpeciesMap = new HashMap<String, String>();
            for(String key : taxonMap.keySet()) {
                for(String val : taxonMap.get(key)) {
                    _allele2SpeciesMap.put(val, key);
                }
            }
        }
        this._numThreads = parallel;
        this._seed = seed;
        this._random = new Random(seed);
        // gene trees
        this._originalGeneTrees = inputTrees;
        try{
            this._inputWeightedGTs = new ArrayList<>();
            for(List<MutableTuple<Tree,Double>> treeList : inputTrees) {
                List<MutableTuple<Tree,Double>> wgtList = new ArrayList<MutableTuple<Tree,Double>>();
                for(MutableTuple<Tree,Double> mt : treeList) {
                    Tree newt = new STITree(mt.Item1.toNewick());
                    wgtList.add(new MutableTuple<Tree, Double>(newt, new Double(mt.Item2)));
                }
                this._inputWeightedGTs.add(wgtList);
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Calculates prior value
     * @param k
     * @return
     */
    public double calculatePrior(double k) {
        return (new NetworkPrior(_speciesNet, k, _startLeaves).calculate());
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
     * Propose a new state based on current state.
     * @return   hastings ratio.
     */
    public double propose() {
        return _operation.propose(_speciesNet);
    }

    public void undo() {
        _operation.undo();
    }

    /**
     * Report the name of operation used
     */
    public String getOperation() {
        return _operation.getOperationName();
    }

    /**
     * gets summarized gene trees and weights
     */
    protected void reportGeneTrees() {
        for(MutableTuple<Tree, Double> tuple : _gtsForStartingNet){
            System.out.println(tuple.Item1.toNewick() + "\t" + tuple.Item2.toString());
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
        String start = Utils.getStartNetwork(_gtsForStartingNet,
                _taxonMap, new HashSet<String>(), _speciesNet);
        if(_printDetails) {
            System.out.println("Get starting network as MDC tree: " + start);
        }
        return Networks.readNetwork(start);
    }
}
