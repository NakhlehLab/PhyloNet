package edu.rice.cs.bioinfo.programs.phylonet.algos.network;


import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 3/25/13
 * Time: 3:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class InferNetworkMLFromGTTBLWithParametricBootstrap extends InferNetworkMLFromGTWithParametricBootstrap {
    protected Func2<Network, Integer, List> getSimulator(final Map<String, List<String>> species2alleles, final String MSPath){
        return new Func2<Network, Integer, List>() {
            public List execute(Network network, Integer numGTs) {
                SimGTInNetworkByMS simulator = new SimGTInNetworkByMS();
                List<List<MutableTuple>> gts = new ArrayList<List<MutableTuple>>();
                for(Tree tr: simulator.generateGTs(network, species2alleles, numGTs, MSPath)){
                    gts.add(Arrays.asList(new MutableTuple(tr, 1.0)));
                }
                return gts;
            }
        };
    }



    protected Func1<List, Tuple<Network,Double>> getInferenceMethod(final Map<String,List<String>> species2alleles, final int maxReticulations, final boolean optimization){
        return new Func1<List, Tuple<Network,Double>>() {
            public Tuple<Network,Double> execute(List gts) {
                InferNetworkMLFromGTTBL_SingleTreePerLocus inference = new InferNetworkMLFromGTTBL_SingleTreePerLocus();
                inference.setSearchParameter(_maxRounds,_maxTryPerBranch,_improvementThreshold,_maxBranchLength,_Brent1,_Brent2,_maxExaminations,_maxFailure,_moveDiameter, _reticulationDiameter, _numThreads, _startNetwork, _fixedHybrid, _operationWeights, _numRuns, _optimizeBL, _seed);
                LinkedList<Tuple<Network,Double>> resultList = new LinkedList<>();
                inference.inferNetwork(gts, species2alleles, maxReticulations, 1, optimization, resultList);
                Tuple<Network,Double> result = resultList.get(0);
                return result;
            }
        };
    }
}
