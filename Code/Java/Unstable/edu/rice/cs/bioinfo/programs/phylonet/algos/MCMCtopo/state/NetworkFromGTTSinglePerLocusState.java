package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.proposal.NetworkProposal;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.List;
import java.util.Map;

/**
 * Created by wendingqiao on 11/2/14.
 */
public class NetworkFromGTTSinglePerLocusState<T> extends NetworkFromGTT {

    public NetworkFromGTTSinglePerLocusState(Network start,
                                             List<List<MutableTuple<Tree, Double>>> inputTrees,
                                             long seed,
                                             int reti,
                                             int parallel,
                                             double[] weights,
                                             Map<String, List<String>> taxonMap
                                             ) {
        super(start, inputTrees, taxonMap, seed, parallel);
        this.calculation = new GTTLikelihood_SinglePerLocus();
        _operation = new NetworkProposal(weights, reti, _seed);
    }

}
