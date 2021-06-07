package edu.rice.cs.bioinfo.programs.phylonet.algos.gibbssampling.util;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Created by yunyu on 10/2/15.
 */
public class PruneNetwork {

    //It removes the reticulation edges that have inheritance probabilities lower than the threshold
    public static void prune(Network network, double threshold){
        List<NetNode> nodeList = new ArrayList<>();
        Set<NetNode> leaves = new HashSet<>();
        for(Object node: Networks.postTraversal(network)){
            nodeList.add((NetNode)node);
            if(((NetNode) node).isLeaf()){
               leaves.add((NetNode)node);
            }
        }
        for(NetNode node: nodeList){
            if(node.getChildCount() == 0 && !leaves.contains(node)){
                node.removeItself();
            }
            else{
                if(node.isNetworkNode()){
                    boolean removed = false;
                    for(Object parentO: node.getParents()){
                        NetNode parentNode = (NetNode)parentO;
                        if(node.getParentProbability(parentNode) <= threshold){
                            parentNode.removeChild(node);
                            removed = true;
                        }
                    }
                    if(removed){
                        node.setParentProbability((NetNode) node.getParents().iterator().next(), NetNode.NO_PROBABILITY);
                    }
                }
            }
        }
        Networks.removeBinaryNodes(network);
    }
}
