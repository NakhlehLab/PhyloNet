package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.NetworkLikelihoodFromGTT_MultiTreesPerLocus;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by wendingqiao on 6/6/15.
 */
public class GTTLikelihood_MultiPerLocus extends NetworkLikelihoodFromGTT_MultiTreesPerLocus {

    public double computeProbability(Network<Object> speciesNetwork, List distinctTrees,
                                     Map<String,List<String>> species2alleles, List gtCorrespondences) {
        return super.computeProbability(speciesNetwork, distinctTrees, species2alleles, gtCorrespondences);
    }


    public void summarizeData(List originalGTs, Map<String,String> allele2species,
                                 List dataForStartingNetwork, List dataForInferNetwork, List treeCorrespondences){
        int treeID = 0;
        Map<String, MutableTuple<Integer,Double>> tree2Info = new HashMap<String, MutableTuple<Integer,Double>>();
        for(Object o: originalGTs) {
            List<MutableTuple<Tree,Double>> treesForOneLocus = (List<MutableTuple<Tree,Double>>)o;
            Map<String, Integer> tree2infoIndex = new HashMap<String, Integer>();
            List<MutableTuple<Integer,Double>> infoList = new ArrayList<MutableTuple<Integer, Double>>();
            for (MutableTuple<Tree, Double> gtTuple : treesForOneLocus) {
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                MutableTuple<Integer,Double> existingInfo = tree2Info.get(exp);
                if (existingInfo == null) {
                    existingInfo = new MutableTuple<Integer,Double>(treeID, gtTuple.Item2);
                    dataForStartingNetwork.add(gtTuple);
                    dataForInferNetwork.add(gtTuple.Item1);
                    tree2Info.put(exp, existingInfo);
                    tree2infoIndex.put(exp, infoList.size());
                    infoList.add(new MutableTuple(treeID, gtTuple.Item2));
                    treeID++;
                } else {
                    existingInfo.Item2 += gtTuple.Item2;
                    Integer infoID = tree2infoIndex.get(exp);
                    if (infoID == null) {
                        tree2infoIndex.put(exp, infoList.size());
                        infoList.add(new MutableTuple(existingInfo.Item1, gtTuple.Item2));
                    } else {
                        infoList.get(infoID).Item2 += gtTuple.Item2;
                    }
                }
            }
            treeCorrespondences.add(infoList);
        }
        for(MutableTuple<Integer,Double> info: tree2Info.values()){
            MutableTuple<Integer,Double> tuple = (MutableTuple<Integer,Double>)dataForStartingNetwork.get(info.Item1);
            tuple.Item2 = info.Item2;
        }
    }

    public String getStartNetwork(List gts, Map<String,List<String>> species2alleles, Set<String> hybridSpecies, Network<Object> startingNetwork){
        Map<String,String> allele2species = null;
        if(species2alleles!=null){
            allele2species = new HashMap<String, String>();
            for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                String species = entry.getKey();
                for(String allele: entry.getValue()){
                    allele2species.put(allele,species);
                }
            }
        }
        MDCInference_DP mdc = new MDCInference_DP();

        List<Solution> solutions = (allele2species==null) ?
                mdc.inferSpeciesTree(gts, false, 1, false, true, -1) :
                mdc.inferSpeciesTree(gts, allele2species, false, 1, false, true, -1);
        Solution sol = solutions.get(0);

        Tree startingTree= Trees.generateRandomBinaryResolution(sol._st);
        startingNetwork = Networks.readNetwork(startingTree.toString());

        for(String hybrid: hybridSpecies){
            createHybrid(startingNetwork, hybrid);
        }
        Networks.removeAllParameters(startingNetwork);
        for(NetNode<Object> node: startingNetwork.dfs()){
            for(NetNode<Object> parent: node.getParents()){
                node.setParentDistance(parent,1.0);
                if(node.isNetworkNode()){
                    node.setParentProbability(parent, 0.5);
                }
            }
        }
        return startingNetwork.toString();
    }

    private void createHybrid(Network<Object> network, String hybrid){
        List<Tuple<NetNode,NetNode>> edgeList = new ArrayList<Tuple<NetNode,NetNode>>();
        Tuple<NetNode,NetNode> destinationEdge = null;
        for(NetNode<Object> node: Networks.postTraversal(network)){
            for(NetNode child: node.getChildren()){
                if(child.isLeaf() && child.getName().equals(hybrid)){
                    if(node.isNetworkNode()){
                        return;
                    }
                    destinationEdge = new Tuple<NetNode, NetNode>(node, child);
                }
                else{
                    edgeList.add(new Tuple<NetNode, NetNode>(node, child));
                }
            }

        }

        int numEdges = edgeList.size();
        Tuple<NetNode,NetNode> sourceEdge = edgeList.get((int)(Math.random() * numEdges));
        NetNode insertedSourceNode = new BniNetNode();
        insertedSourceNode.adoptChild(sourceEdge.Item2, NetNode.NO_DISTANCE);
        sourceEdge.Item1.removeChild(sourceEdge.Item2);
        sourceEdge.Item1.adoptChild(insertedSourceNode, NetNode.NO_DISTANCE);
        NetNode insertedDestinationNode = new BniNetNode();
        insertedDestinationNode.adoptChild(destinationEdge.Item2, NetNode.NO_DISTANCE);
        destinationEdge.Item1.removeChild(destinationEdge.Item2);
        destinationEdge.Item1.adoptChild(insertedDestinationNode, NetNode.NO_DISTANCE);
        insertedSourceNode.adoptChild(insertedDestinationNode, NetNode.NO_DISTANCE);
    }


}
