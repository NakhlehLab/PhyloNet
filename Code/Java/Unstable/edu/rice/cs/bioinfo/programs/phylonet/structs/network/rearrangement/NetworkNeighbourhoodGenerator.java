package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;



import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by Yun Yu
 *
 * This class is for generating a neighbor of a species network by mutating it
 */
public abstract class NetworkNeighbourhoodGenerator{
    protected boolean _printDetails = false;

    /**
     * This function is to undo the last rearrangement move
     */
    public abstract void undo();


    /**
     * This function is to mutate the current network in order to generate its neighbor
     */
    public abstract void mutateNetwork(Network network);


    /**
     * This function is to print an edge
     */
    protected String printEdge(Tuple<NetNode, NetNode> edge) {
        return "(" + edge.Item1.getName() + "," + edge.Item2.getName() + ")";
    }


    /**
     * This function is to check if the given species network is valid
     *
     * @param network   the given species network
     * @param leaves    the set of leaves that the given network <code>network</code> should have
     */
    protected boolean isNetworkValid(Network network, Set<String> leaves){
        int count = 0;
        for(Object leaf: network.getLeaves()){
            if(leaves.contains(((NetNode)leaf).getName())){
                count++;
            }
            else{
                return false;
            }
        }
        if(count!=leaves.size())return false;
        if(!Networks.isDisconnectedNetwork(network,null))return false;
        for(Object node: Networks.postTraversal(network)){
            double totalProb = 0;
            for (Object parent : ((NetNode) node).getParents()) {
                totalProb += ((NetNode) node).getParentProbability((NetNode) parent);
            }
            if(((NetNode)node).getChildCount()==1 && ((NetNode)node).getParentCount()<2){
                return false;
            }
            if(totalProb!=NetNode.NO_PROBABILITY && ((NetNode)node).isNetworkNode()){
                if(Math.abs(totalProb - 1) > 0.00001) {
                    throw new RuntimeException(network.toString());
                }
            }
            else if(!((NetNode)node).isRoot()){
                if(totalProb != NetNode.NO_PROBABILITY){
                    throw new RuntimeException(network.toString());
                }
            }
        }
        return true;
    }


    /**
     * This function is to check if the given species network has a set of species as hybrids
     *
     * @param network   the given species network
     * @param requiredHybridSet    the set of taxa that should be placed as hybrid species
     */
    protected boolean checkHybrid(Network network, Set<String> requiredHybridSet){
        if(requiredHybridSet==null) return true;
        Map<NetNode, Set<NetNode>> node2children = new HashMap<>();
        Set<NetNode> reticulationNodes = new HashSet<>();
        for(Object o: Networks.postTraversal(network)){
            NetNode node = (NetNode)o;
            Set<NetNode> childrenSet = new HashSet<>();
            node2children.put(node, childrenSet);
            if(node.isNetworkNode()){
                reticulationNodes.add(node);
            }
            for(Object child: node.getChildren()){
                childrenSet.add((NetNode)child);
                childrenSet.addAll(node2children.get(child));
            }
        }
        Set<String> hybridSet = new HashSet<>();
        for(NetNode reticulation: reticulationNodes){
            for(NetNode child: node2children.get(reticulation)){
                if(child.isLeaf()){
                    if(!requiredHybridSet.contains(child.getName())){
                        return false;
                    }
                    else{
                        hybridSet.add(child.getName());
                    }
                }
            }
        }
        return hybridSet.size() == requiredHybridSet.size();
    }



    /**
     * This function is to reset the list of edges in start tree
     */
    public abstract void resetList();

    public abstract int getOperationID();

}