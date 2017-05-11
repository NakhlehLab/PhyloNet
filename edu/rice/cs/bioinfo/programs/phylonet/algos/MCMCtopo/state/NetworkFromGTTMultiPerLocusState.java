package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.proposal.NetworkProposal;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.*;

/**
 * Created by wendingqiao on 11/2/14.
 */
public class NetworkFromGTTMultiPerLocusState<T> extends NetworkFromGTT {

    public NetworkFromGTTMultiPerLocusState(Network start,
                                            List<List<MutableTuple<Tree, Double>>> inputTrees,
                                            long seed,
                                            int reti,
                                            int parallel,
                                            double[] weights,
                                            Map<String, List<String>> taxonMap
    ) {
        super(start, inputTrees, taxonMap, seed, parallel);
        this.calculation = new GTTLikelihood_MultiPerLocus();
        _operation = new NetworkProposal(weights, reti, _seed);
    }

}
