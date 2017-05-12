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

    /**
     * Constructor
     */
    public NetworkFromGTT(Map<String, List<String>> taxonMap,
                          long seed,
                          int parallel) {
        this._taxonMap = taxonMap;
        this._allele2SpeciesMap = new HashMap<String, String>();
        for(String key : taxonMap.keySet()) {
            for(String val : taxonMap.get(key)) {
                _allele2SpeciesMap.put(val, key);
            }
        }
        this._numThreads = parallel;
        this._seed = seed;
        this._random = new Random(seed);
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
}
