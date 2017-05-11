package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.state;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_Rooted;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by dw20 on 5/11/17.
 */
public class Utils {

    public static String getStartNetwork(List gts, Map<String,List<String>> species2alleles,
                                  Set<String> hybridSpecies, Network<Object> startingNetwork){
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
        MDCInference_Rooted mdc = new MDCInference_Rooted();

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

    private static void createHybrid(Network<Object> network, String hybrid){
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


    public static void setBranchLengths(Network net) {
        for(Object n1 : net.bfs()) {
            NetNode node1 = (NetNode)n1;
            for(Object n2 : node1.getParents()) {
                NetNode node2 = (NetNode)n2;
                if( Double.isNaN( node1.getParentDistance(node2) ) ) {
                    node1.setParentDistance(node2, 1.0);
                }
                if(node1.isNetworkNode() && Double.isNaN( node1.getParentProbability(node2) ) ) {
                    node1.setParentProbability(node2, 0.5);
                }
            }
        }
    }

}
