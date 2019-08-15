package edu.rice.cs.bioinfo.programs.phylonet.algos.summarize;

/*
 *@ClassName: tripartition2
 *@Description: This class is to get the tripartitions and their support for a set of networks
 *              tripartition for $node$: X:Y|Z,
 *              X is the set of leaves under $node$,
 *              Y, Z are the sets of leaves under two siblings of $node$. If only one sibling, just Y
 *@Author: Zhen Cao
 *@Date:  2019-07-26 20:26
 *@Version: 1.0
 */

import com.google.common.collect.Lists;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

public class tripartition2 {

    private static boolean _debug = false;
    
    
    
    
    /**
     * @Description:    This function is to get the tripartitions and their support for a set of networks
     * @Param: networkMap
     * @Author: Zhen Cao
     * @Date: 2019-08-15 
     */

    public static Map<netNodeTuple, Double> summarize(Map<Network, Double> networkMap){
        Map<netNodeTuple, Double> networkNodeTupleMap = new HashMap<>();
        for(Network network: networkMap.keySet()){
            Networks.autoLabelNodes(network);

            if (_debug) System.out.println(network.toString());
            double count = networkMap.get(network);

            Set<netNodeTuple> retiSet = networkTripartition(network);
            for(netNodeTuple nodetuple: retiSet){
                boolean inMap = false;

                for (netNodeTuple nodetuple1: networkNodeTupleMap.keySet()){
                    if(nodetuple.equals(nodetuple1)){
                        double num = networkNodeTupleMap.get(nodetuple1);
                        networkNodeTupleMap.remove(nodetuple1);
                        networkNodeTupleMap.put(nodetuple1, num+count);
                        inMap = true;
                        break;
                    }
                }
                if(!inMap){
                    networkNodeTupleMap.put(nodetuple, count);
                }
            }
        }
        if (_debug) {
            for (netNodeTuple netnodetuple : networkNodeTupleMap.keySet()) {
                System.out.println(netnodetuple.toString());
                System.out.println(networkNodeTupleMap.get(netnodetuple));
            }
        }
        return networkNodeTupleMap;

    }


    
    /**
     * @Description:    This function is to get all the tripartitions for a network
     * @Param: network
     * @Author: Zhen Cao
     * @Date: 2019-08-15 
     */


    public static Set<netNodeTuple> networkTripartition(Network<Integer> network){
        Set<netNodeTuple> retiSet = new HashSet<>();
        for(NetNode<Integer> netnode: network.getNetworkNodes()){
            netNodeTuple tuple = netNodeTripartition(netnode);
            retiSet.add(tuple);

        }
        return retiSet;
    }
    
    
    /**
     * @Description:    This function is to get the tripartition for a reticulation node $node$
     * @Param: node
     * @Author: Zhen Cao
     * @Date: 2019-08-15 
     */


    public static netNodeTuple netNodeTripartition(NetNode<Integer> node){
        Set<String> retiLeaves = getLeavesUnderNode(node);
        Set<Set<String>> siblingSet = new HashSet<>();

        for(NetNode<Integer> parent: node.getParents()){
            for(NetNode<Integer> child: parent.getChildren()){
                if(!child.equals(node)){
                    siblingSet.add(getLeavesUnderNode(child));
                }
            }
        }

        return new netNodeTuple(node, siblingSet, retiLeaves);
    }
    
    
    
    
    /**
     * @Description:  This function is to get the set of leaves that can be reached from a $node$
     * @Param:  node
     * @Author: Zhen Cao
     * @Date: 2019-08-15 
     */


    public static Set<String> getLeavesUnderNode(NetNode<Integer> node){
        Set<String> leaves = new HashSet<>();
        BniNetNode<Integer> newnode = new BniNetNode();
        if(node.isLeaf()){
            leaves.add(node.getName());
        }
        else {
            for (NetNode<Integer> child:node.getChildren()){
                newnode.adoptChild(child, NetNode.NO_DISTANCE);
            }

            BniNetwork<Integer> newNetwork = new BniNetwork(newnode);
            for(NetNode<Integer> leaf:newNetwork.getLeaves()){
                leaves.add(leaf.getName());
            }

            List<NetNode<Integer>> newnodechildren = Lists.newArrayList(newnode.getChildren());
            for (NetNode<Integer> child:newnodechildren){
                newnode.removeChild(child);
            }
        }


        return leaves;
    }

    public static void main(String[] args) {

        String n = "(((A,D),(B)X#H1),(X#H1,C));";

        Map<Network, Double> networkMap = new HashMap<>();
        networkMap.put(Networks.readNetwork(n), 0.5);
        summarize(networkMap);



    }

}
