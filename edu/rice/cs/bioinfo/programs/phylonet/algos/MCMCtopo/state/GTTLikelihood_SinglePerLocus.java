package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.NetworkLikelihoodFromGTT_SingleTreePerLocus;
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
public class GTTLikelihood_SinglePerLocus extends NetworkLikelihoodFromGTT_SingleTreePerLocus {

    public double computeProbability(Network<Object> speciesNetwork, List distinctTrees,
                                     Map<String,List<String>> species2alleles, List gtCorrespondences) {
        return super.computeProbability(speciesNetwork, distinctTrees, species2alleles, gtCorrespondences);
    }


    public void summarizeData(List originalGTs, Map<String,String> allele2species,
                                 List gtsForStartingNetwork, List gtsForInferNetwork, List treeCorrespondences){
        Map<String, MutableTuple<Tree,Double>> exp2tree = new HashMap<String, MutableTuple<Tree, Double>>();
        Map<String,Integer> exp2ID = new HashMap<String, Integer>();
        for(Object list: originalGTs) {
            for (MutableTuple<Tree, Double> gtTuple : (List<MutableTuple<Tree, Double>>)list) {
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                MutableTuple<Tree, Double> existingTuple = exp2tree.get(exp);
                if (existingTuple == null) {
                    existingTuple = gtTuple;
                    exp2tree.put(exp, existingTuple);
                    Set<Integer> binaryIDs = new HashSet<Integer>();
                    for (Tree btr : Trees.getAllBinaryResolution(gtTuple.Item1)) {
                        String btrExp = Trees.getLexicographicNewickString(btr, allele2species);
                        Integer index = exp2ID.get(btrExp);
                        if (index == null) {
                            index = gtsForInferNetwork.size();
                            gtsForInferNetwork.add(btr);
                            binaryIDs.add(index);
                            exp2ID.put(btrExp, index);
                        } else {
                            binaryIDs.add(index);
                        }
                    }
                    treeCorrespondences.add(new Tuple<MutableTuple<Tree, Double>, Set<Integer>>(gtTuple, binaryIDs));

                } else {
                    existingTuple.Item2 += gtTuple.Item2;
                }

            }
        }

        gtsForStartingNetwork.addAll(exp2tree.values());

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
