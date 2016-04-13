package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.proposal.NetworkProposal;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.*;

/**
 * Created by wendingqiao on 11/2/14.
 */
public class NetworkFromGTTMultiPerLocusState<T> extends NetworkFromGTTMultiPerLocus {

    // operations
    protected NetworkProposal _operation;

    public NetworkFromGTTMultiPerLocusState(Network start,
                                            List<List<MutableTuple<Tree, Double>>> inputTrees,
                                            long seed,
                                            int reti,
                                            int parallel,
                                            double[] weights,
                                            Map<String, List<String>> taxonMap
    ) {
        super(start, inputTrees, taxonMap, seed, parallel);
        _operation = new NetworkProposal(weights, reti, _seed);
    }

    /**
     * Propose a new state based on current state.
     * @return   hastings ratio.
     */
    public double propose() {
        //System.out.println(_speciesNet.toString());
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
