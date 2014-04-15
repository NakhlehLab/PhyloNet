package edu.rice.cs.bioinfo.programs.phylonet.algos.network;


import edu.rice.cs.bioinfo.library.phylogenetics.parametricbootstrap.network.NetworkParametricBootstrap;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkBootstrap;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 3/25/13
 * Time: 3:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class InferILSNetworkMLBootstrap extends InferILSNetworkProbabilisticallyParallel {

    public Tuple<Network,Double> inferNetworkWithBootstrap(List<MutableTuple<Tree,Double>> gts, Map<String,List<String>> species2alleles, int maxReticulations, int numRounds){
        Tuple<Network, Double> result = this.inferNetwork(gts, species2alleles, maxReticulations, 1).get(0);
        //System.out.println("\n"+network2String(result.Item1) + " : " + result.Item2);
        Func2<Network, Integer, List<MutableTuple<Tree,Double>>> simulator = getSimulator(species2alleles);
        Func2<Network, List<Network>, Network> estimator = getEstimator();
        Func1<List<MutableTuple<Tree,Double>>, Network> inference = getInferenceMethod(species2alleles, maxReticulations);
        NetworkParametricBootstrap<Network, List<MutableTuple<Tree,Double>>, Network, Integer> networkBootstrap = new NetworkParametricBootstrap<Network, List<MutableTuple<Tree,Double>>, Network, Integer>(simulator, inference, estimator);
        Network inferredNetwork = networkBootstrap.performBootstrapping(result.Item1, gts.size(), numRounds);
        return new Tuple<Network, Double>(inferredNetwork, result.Item2);
    }


    private Func2<Network, Integer, List<MutableTuple<Tree,Double>>> getSimulator(final Map<String, List<String>> species2alleles){
        return new Func2<Network, Integer, List<MutableTuple<Tree,Double>>>() {
            public List<MutableTuple<Tree,Double>> execute(Network network, Integer numGTs) {
                SimGTInNetwork simulator = new SimGTInNetwork();
                List<MutableTuple<Tree,Double>> gts = new ArrayList<MutableTuple<Tree, Double>>();
                for(Tree tr: simulator.generateGTs(network, species2alleles, numGTs)){
                    gts.add(new MutableTuple<Tree, Double>(tr, 1.0));
                }
                return gts;
            }
        };
    }

    private Func2<Network, Integer, List<MutableTuple<Tree,Double>>> getMSSimulator(final Map<String, List<String>> species2alleles){
        return new Func2<Network, Integer, List<MutableTuple<Tree,Double>>>() {
            public List<MutableTuple<Tree,Double>> execute(Network network, Integer numGTs) {
                SimGTInNetworkByMS simulator = new SimGTInNetworkByMS();
                List<MutableTuple<Tree,Double>> gts = new ArrayList<MutableTuple<Tree, Double>>();
                for(Tree tr: simulator.generateGTs(network, species2alleles, numGTs, "/research/tools/MS/msdir/")){
                    gts.add(new MutableTuple<Tree, Double>(tr, 1.0));
                }
                return gts;
            }
        };
    }

    private Func2<Network, List<Network>, Network> getEstimator(){
        return new Func2<Network, List<Network>, Network>() {
            public Network execute(Network targetNet, List<Network> refNets) {
                NetworkBootstrap nb = new NetworkBootstrap();
                targetNet = nb.computeNetworkBranchSupport(targetNet, refNets);
                return targetNet;
            }
        };
    }

    private Func1<List<MutableTuple<Tree,Double>>, Network> getInferenceMethod(final Map<String,List<String>> species2alleles, final int maxReticulations){
        return new Func1<List<MutableTuple<Tree,Double>>, Network>() {
            public Network execute(List<MutableTuple<Tree,Double>> gts) {
                InferILSNetworkProbabilisticallyParallel inference = new InferILSNetworkProbabilisticallyParallel();
                inference.setSearchParameter(_maxRounds,_maxTryPerBranch,_improvementThreshold,_maxBranchLength,_Brent1,_Brent2,_maxExaminations,_maxFailure,_diameterLimit, _numThreads, _startNetwork, _fixedHybrid, _operationWeight, _numRuns, _seed);
                Network network = inference.inferNetwork(gts, species2alleles, maxReticulations, 1).get(0).Item1;
                System.out.println(network2String(network));
                return network;
            }
        };
    }
}
