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
 * Created by Yun Yu
 * Date: 3/25/13
 * Time: 3:58 PM
 *
 * This class is a subclass of InferNetworkMLFromGTWithParametricBootstrap. It uses only the topologies of the gene trees to infer species network
 */
public class InferNetworkMLFromGTTWithParametricBootstrap extends InferNetworkMLFromGTWithParametricBootstrap{

    /**
     * This function is to obtain the simulator for generating gene trees down a species network
     * In this case, an internal simulator is used which generates gene trees without branch lengths
     *
     * @param species2alleles    mapping from species to alleles sampled from it
     * @param simulatorPath      not used
     */
    protected Func2<Network, Integer, List> getSimulator(final Map<String, List<String>> species2alleles, final String simulatorPath){
        return new Func2<Network, Integer, List>() {
            public List execute(Network network, Integer numGTs) {
                SimGTInNetwork simulator = new SimGTInNetwork();
                simulator.setSeed(_seed);
                if(_seed != null)_seed = _seed * 2;
                List<List<MutableTuple>> gts = new ArrayList<List<MutableTuple>>();
                for(Tree tr: simulator.generateGTs(network, species2alleles, numGTs)){
                    gts.add(Arrays.asList(new MutableTuple(tr, 1.0)));
                }
                return gts;
            }
        };
    }



    /**
     * This function is to obtain the method used for inferring species network
     * Only single gene tree per locus is considered
     *
     * @param species2alleles    mapping from species to alleles sampled from it
     * @param maxReticulations   the maximal number of networks examined during the search; can be used as one of the termination criterion of the search
     */
    protected Func1<List, Tuple<Network,Double>> getInferenceMethod(final Map<String,List<String>> species2alleles, final int maxReticulations, final boolean postOptimize){
        return new Func1<List, Tuple<Network,Double>>() {
            public Tuple<Network,Double> execute(List gts) {
                InferNetworkMLFromGTT_SingleTreePerLocus inference = new InferNetworkMLFromGTT_SingleTreePerLocus();
                inference.setSearchParameter(_maxRounds,_maxTryPerBranch,_improvementThreshold,_maxBranchLength,_Brent1,_Brent2,_maxExaminations,_maxFailure,_moveDiameter,_reticulationDiameter, _numThreads, _startNetwork, _fixedHybrid, _operationWeights, _numRuns, _optimizeBL, _seed);
                LinkedList<Tuple<Network,Double>> resultList = new LinkedList<>();
                inference.inferNetwork(gts, species2alleles, maxReticulations, 1, postOptimize, resultList);
                Tuple<Network,Double> result = resultList.get(0);
                return result;
            }
        };
    }
}
