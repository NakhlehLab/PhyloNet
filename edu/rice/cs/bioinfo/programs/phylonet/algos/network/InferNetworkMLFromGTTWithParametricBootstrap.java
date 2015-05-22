package edu.rice.cs.bioinfo.programs.phylonet.algos.network;


import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkBootstrap;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 3/25/13
 * Time: 3:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class InferNetworkMLFromGTTWithParametricBootstrap extends InferNetworkMLFromGTWithParametricBootstrap{
    protected Func2<Network, Integer, List> getSimulator(final Map<String, List<String>> species2alleles, final String simulatorPath){
        return new Func2<Network, Integer, List>() {
            public List execute(Network network, Integer numGTs) {
                SimGTInNetwork simulator = new SimGTInNetwork();
                simulator.setSeed(_seed);
                _seed = _seed * 2;
                List<List<MutableTuple>> gts = new ArrayList<List<MutableTuple>>();
                for(Tree tr: simulator.generateGTs(network, species2alleles, numGTs)){
                    gts.add(Arrays.asList(new MutableTuple(tr, 1.0)));
                }
                return gts;
            }
        };
    }



    protected Func1<List, Tuple<Network,Double>> getInferenceMethod(final Map<String,List<String>> species2alleles, final int maxReticulations){
        return new Func1<List, Tuple<Network,Double>>() {
            public Tuple<Network,Double> execute(List gts) {
                InferNetworkMLFromGTT_SingleTreePerLocus inference = new InferNetworkMLFromGTT_SingleTreePerLocus();
                inference.setSearchParameter(_maxRounds,_maxTryPerBranch,_improvementThreshold,_maxBranchLength,_Brent1,_Brent2,_maxExaminations,_maxFailure,_moveDiameter,_reticulationDiameter, _numThreads, _startNetwork, _fixedHybrid, _operationWeight, _numRuns, _seed);
                LinkedList<Tuple<Network,Double>> resultList = new LinkedList<>();
                inference.inferNetwork(gts, species2alleles, maxReticulations, 1, resultList);
                Tuple<Network,Double> result = resultList.get(0);
                //System.out.println(network2String(network));
                return result;
            }
        };
    }
}
