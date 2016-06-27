package edu.rice.cs.bioinfo.programs.phylonet.algos.gibbssampling;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF_Cached;
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
public abstract class GibbsSamplingForPruningNetworksFromGTT extends GibbsSamplingForPruningNetworks {

    protected List<Tuple<NetNode,NetNode>> getAllParametersToSample(Network network, List<Tree> gts, Map<String, String> allele2species){
        Set<String> singleAlleleSpecies = new HashSet<>();
        if(allele2species==null){
            for(Tree gt: gts){
                for(String leaf: gt.getLeaves()){
                    singleAlleleSpecies.add(leaf);
                }
            }
        }
        else{
            singleAlleleSpecies.addAll(allele2species.values());
            for(Object tuple: gts){
                Tree gt = ((MutableTuple<Tree,Double>)tuple).Item1;
                if(singleAlleleSpecies.size()==0)break;
                Set<String> speciesVisited = new HashSet<>();

                for(String leaf: gt.getLeaves()){
                    String species = allele2species.get(leaf);
                    if(singleAlleleSpecies.contains(species)){
                        if(speciesVisited.contains(species)){
                            singleAlleleSpecies.remove(species);
                        }
                        else{
                            speciesVisited.add(species);
                        }
                    }
                }

            }
        }

        List<Tuple<NetNode,NetNode>> parametersToSample = new ArrayList<>();

        Map<NetNode, Set<String>> node2leaves = new HashMap<>();
        for (Object nodeO : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeO;
            Set<String> leaves = new HashSet<>();
            if(node.isLeaf()){
                leaves.add(node.getName());
            }
            for (Object childO : node.getChildren()) {
                NetNode childNode = (NetNode) childO;
                //((NetNode) childO).setParentDistance(node, 1.0);
                Set<String> childLeaves = node2leaves.get(childNode);
                leaves.addAll(childLeaves);
                Tuple<NetNode, NetNode> edge = new Tuple<>(node, childNode);
                if (childLeaves.size() != 1 || !singleAlleleSpecies.containsAll(childLeaves)) {
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


    private double computeLikelihood(Network network, List distinctGTs, List treeCorrespondences, Tuple<NetNode, NetNode> editedEdge) {
        Set<NetNode> childNodes = new HashSet<NetNode>();
        childNodes.add(editedEdge.Item2);
        Set<NetNode> parentNodes = new HashSet<NetNode>();
        if(editedEdge.Item1==null){
            for(Object parentNode: editedEdge.Item2.getParents()){
                parentNodes.add((NetNode)parentNode);
            }
        }
        else{
            parentNodes.add(editedEdge.Item1);
        }

        GeneTreeProbabilityYF_Cached gtp = new GeneTreeProbabilityYF_Cached();
        double[] resultProbs = new double[distinctGTs.size()];
        Thread[] myThreads = new Thread[_numProcessors];
        if(_numProcessors == 1) {
            gtp.calculateGTDistribution(network, distinctGTs, childNodes, parentNodes, resultProbs);
        }
        else{
            gtp.setParallel(true);
            gtp.preProcess(network, distinctGTs, false);

            for(int i=0; i<_numProcessors; i++){
                myThreads[i] = new MyThreadFromNonScratch(gtp, network, distinctGTs, childNodes, parentNodes, resultProbs);
                myThreads[i].start();
            }

            for(int i=0; i<_numProcessors; i++){
                try {
                    myThreads[i].join();
                } catch (InterruptedException ignore) {}
            }

        }

        double probability = computeLikelihood(resultProbs, treeCorrespondences);
        return probability;
    }


    protected double computeInitialLikelihood(Network network, List distinctGTs, Map<String, List<String>> species2alleles, List treeCorrespondences){
        GeneTreeProbabilityYF_Cached gtp = new GeneTreeProbabilityYF_Cached();
        double[] resultProbs = new double[distinctGTs.size()];
        if(_numProcessors == 1){
            gtp.calculateGTDistribution(network, distinctGTs, species2alleles, resultProbs);
        }
        else{
            gtp.setParallel(true);
            gtp.preProcess(network, distinctGTs, true);
            Thread[] myThreads = new Thread[_numProcessors];
            for(int i=0; i<_numProcessors; i++){
                myThreads[i] = new MyThreadFromScratch(gtp, network, distinctGTs, species2alleles, resultProbs);
                myThreads[i].start();
            }

            for(int i=0; i<_numProcessors; i++){
                try {
                    myThreads[i].join();
                } catch (InterruptedException ignore) {}
            }
        }
        double prob = computeLikelihood(resultProbs, treeCorrespondences);
        return prob;
    }


    protected Func<Double> likelihoodUpdater(Network network, List distinctGTs, Map<String, List<String>> species2alleles,List treeCorrespondences, Tuple<NetNode, NetNode> editedEdge){
        return new Func() {
            @Override
            public Double execute() {
                double lnLikelihood = computeLikelihood(network, distinctGTs, treeCorrespondences, editedEdge);
                return lnLikelihood;
            }
        };
    }


    private class MyThreadFromScratch extends Thread{
        GeneTreeProbabilityYF_Cached _gtp;
        Network _speciesNetwork;
        List<Tree> _geneTrees;
        Map<String, List<String>> _species2alleles;
        double[] _probs;


        public MyThreadFromScratch(GeneTreeProbabilityYF_Cached gtp, Network speciesNetwork, List<Tree> geneTrees, Map<String, List<String>> species2alleles, double[] probs){
            _speciesNetwork = speciesNetwork;
            _geneTrees = geneTrees;
            _species2alleles = species2alleles;
            _probs = probs;
            _gtp = gtp;
        }


        public void run() {
            _gtp.calculateGTDistribution(_speciesNetwork, _geneTrees, _species2alleles, _probs);

        }
    }





    private class MyThreadFromNonScratch extends Thread{
        Network _speciesNetwork;
        List<Tree> _gts;
        double[] _probs;
        Set<NetNode> _childNodes;
        Set<NetNode> _parentNodes;
        GeneTreeProbabilityYF_Cached _gtp;


        public MyThreadFromNonScratch(GeneTreeProbabilityYF_Cached gtp, Network speciesNetwork, List<Tree> gts, Set<NetNode> childNodes, Set<NetNode> parentNodes, double[] probs){
            _speciesNetwork = speciesNetwork;
            _gts = gts;
            _probs = probs;
            _childNodes = childNodes;
            _parentNodes = parentNodes;
            _gtp = gtp;
        }


        public void run() {
            _gtp.calculateGTDistribution(_speciesNetwork, _gts, _childNodes, _parentNodes, _probs);

        }
    }

    protected abstract double computeLikelihood(double[] probList, List gtCorrespondences);

}
