package edu.rice.cs.bioinfo.programs.phylonet.algos.gibbssampling;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityPseudo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF_Cached;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by yunyu on 10/1/15.
 */
public abstract class GibbsSamplingForPruningNetworksFromGTTPseudo extends GibbsSamplingForPruningNetworks {

    protected List<Tuple<NetNode,NetNode>> getAllParametersToSample(Network network, List<Tree> gts, Map<String, String> allele2species) {
        List<Tuple<NetNode, NetNode>> parametersToSample = new ArrayList<>();

        Map<NetNode, Set<String>> node2leaves = new HashMap<>();
        for (Object nodeO : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeO;
            Set<String> leaves = new HashSet<>();
            if (node.isLeaf()) {
                leaves.add(node.getName());
            }
            for (Object childO : node.getChildren()) {
                NetNode childNode = (NetNode) childO;
                Set<String> childLeaves = node2leaves.get(childNode);
                leaves.addAll(childLeaves);
                Tuple<NetNode, NetNode> edge = new Tuple<>(node, childNode);
                if (childLeaves.size() != 1) {
                    parametersToSample.add(edge);
                }
            }
            if (node.isNetworkNode()) {
                Tuple<NetNode, NetNode> edge = new Tuple<>(null, node);
                parametersToSample.add(edge);
            }
            node2leaves.put(node, leaves);
        }
        return parametersToSample;

    }


    protected double computeInitialLikelihood(Network network, List allTriplets, Map<String, List<String>> species2alleles, List tripleFrequencies){
        GeneTreeProbabilityPseudo calculator = new GeneTreeProbabilityPseudo();
        if(_numProcessors!=0){
            int batchSize = allTriplets.size()/_numProcessors;
            if(allTriplets.size()%_numProcessors != 0){
                batchSize++;
            }
            calculator.setBatchSize(batchSize);
        }
        calculator.initialize(network);
        double[][] probs = new double[allTriplets.size()][3];
        Thread[] myThreads = new Thread[_numProcessors];
        //System.out.println(speciesNetwork);
        if(_numProcessors>1) {
            calculator.setParallel(true);
            for (int i = 0; i < _numProcessors; i++) {
                myThreads[i] = new MyThread(network, calculator, allTriplets, probs);
                myThreads[i].start();
            }
            for (int i = 0; i < _numProcessors; i++) {
                try {
                    myThreads[i].join();
                } catch (InterruptedException ignore) {
                }
            }
        }else{
            try {
                calculator.computePseudoLikelihood(network, allTriplets, probs);
            }catch (Exception e){
                System.out.println(network);
                System.err.println(e.getMessage());
                e.getStackTrace();
                System.exit(-1);

            }
        }
        double totalProb = computeLikelihood(probs, tripleFrequencies);
        return totalProb;
    }


    protected Func<Double> likelihoodUpdater(Network network, List triplets, Map<String, List<String>> species2alleles,List tripletFrequencies, Tuple<NetNode, NetNode> editedEdge){
        return new Func() {
            @Override
            public Double execute() {
                double lnLikelihood = computeInitialLikelihood(network, triplets, null, tripletFrequencies);
                return lnLikelihood;
            }
        };
    }


    protected class MyThread extends Thread{
        Network _network;
        GeneTreeProbabilityPseudo _calculator;
        List<String> _allTriplets;
        double[][] _probs;


        public MyThread(Network network,GeneTreeProbabilityPseudo calculator,List<String> allTriplets, double[][] probs){
            _network = network;
            _calculator = calculator;
            _allTriplets = allTriplets;
            _probs = probs;
        }

        public void run(){
            _calculator.computePseudoLikelihood(_network, _allTriplets, _probs);
        }
    }

     protected abstract double computeLikelihood(double[][] probs, List tripletFrequencies);

}
