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
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BfsSearch;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.uci.ics.jung.algorithms.shortestpath.BFSDistanceLabeler;

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
            if (_debug) System.out.println("netnode:"+netnode.getName());
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
        try{
            for(NetNode<Integer> parent: node.getParents()){
                for(NetNode<Integer> child: parent.getChildren()){
                    if(!child.equals(node)){
                        siblingSet.add(getLeavesUnderNode(child));
                    }
                }
            }
        }catch (Exception e){

            e.printStackTrace();
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
        if(node.isLeaf()){
            leaves.add(node.getName());
        }
        else {

            for (NetNode<Integer> nn : new BfsSearch<Integer>(node)){
                if (nn.isLeaf()){
                    leaves.add(nn.getName());
                }
            }
        }


        return leaves;
    }
//    public static Set<String> getLeavesUnderNode(NetNode<Integer> node){
//        Set<String> leaves = new HashSet<>();
//        BniNetNode<Integer> newnode = new BniNetNode();
//        if(node.isLeaf()){
//            leaves.add(node.getName());
//        }
//        else {
//            for (NetNode<Integer> child:node.getChildren()){
//                newnode.adoptChild(child, NetNode.NO_DISTANCE);
//            }
//
//            BniNetwork<Integer> newNetwork = new BniNetwork(newnode);
//            for(NetNode<Integer> leaf:newNetwork.getLeaves()){
//                leaves.add(leaf.getName());
//            }
//
//            List<NetNode<Integer>> newnodechildren = Lists.newArrayList(newnode.getChildren());
//            for (NetNode<Integer> child:newnodechildren){
//                newnode.removeChild(child);
//            }
//        }
//
//
//        return leaves;
//    }

    public static void main(String[] args) {

        String n1 = "((L:1.0)#H1:1.0::0.6775044779684056,((F:1.0,((K:1.0,(P:1.0)#H3:1.0::0.6753614616352684):1.1244258148662702)#H2:1.5700202905515581::0.6135769094686397):5.938615025704547,(((#H1:1.0::0.3224955220315944)#H4:1.0::0.7475384042289634,(#H2:5.938905009947807::0.38642309053136026,(C:1.0,#H4:1.0::0.2524615957710365):5.913620320644154):5.912120605817229):5.94018448223261,(#H3:1.0::0.32463853836473155,O:1.0):2.8047977192873814):5.903825230963656):5.93953170330228);";
        String n2 = "((L:1.0)#H1:1.0::0.697819639034568,((F:1.0,((K:1.0,(P:1.0)#H3:1.0::0.6745110428348149):1.1100661168406565)#H2:1.7082666741217403::0.6142669809548641):5.931845445367439,((#H3:1.0::0.3254889571651851,O:1.0):2.4833570014727955,(#H2:5.929217003496626::0.3857330190451359,((#H1:1.0::0.302180360965432)#H4:1.0::0.017142685450424855,(C:1.0,#H4:1.0::0.9828573145495751):5.91361935198872):5.9135981328244815):5.925654204620761):5.9260404685243415):5.937215530824061);";
        String n3 = "((L:1.0)#H1:1.0::0.697176914245199,(((((K:1.0,(P:1.0)#H3:1.0::0.6779314464962011):1.0884181724604034)#H2:5.934653324530301::0.38618142938381267,(C:1.0,#H1:1.0::0.302823085754801):5.9084237144406755):5.9397822589449465,((#H3:1.0::0.3220685535037989)#H4:1.0::0.017142685450424855,(#H4:1.0::0.9828573145495751,O:1.0):0.00808293691880342):2.7714903993551276):5.904173024132873,(#H2:2.012640927894127::0.6138185706161874,F:1.0):5.904657745927051):5.937801402829414);";
        String n4 = "((L:1.0)#H1:1.0::0.6712728892686329,(((#H1:1.0::0.32872711073136707,(((K:1.0,(P:1.0)#H3:1.0::0.6751112732063029):1.1267080674211265)#H2:5.938355167858189::0.38140833857114587,C:1.0):5.911170667122007):5.939107582827078,(#H3:1.0::0.32488872679369707,O:1.0):2.860440724320893):5.926093346898221,(F:1.0,#H2:1.6125835293541366::0.6185916614288541):5.9381572434685745):5.939517506695277);";
        String n5 = "((L:1.0)#H1:1.0::0.6966745198922524,((F:1.0,(((P:1.0)#H3:1.0::0.683073175568802,K:1.0):1.0524332659314857)#H2:2.5623058987490537::0.6104433712226643):5.906687726566634,((#H3:1.0::0.316926824431198,O:1.0):3.4112578018682806,(#H2:5.937243226797915::0.38955662877733566,(C:1.0,#H1:1.0::0.3033254801077477):5.909482135710924):5.909111928318673):5.9042519002108):5.9371033552041155);";

        Map<Network, Double> networkMap = new HashMap<>();
        networkMap.put(Networks.readNetwork(n1), 1.0);
        networkMap.put(Networks.readNetwork(n2), 1.0);
        networkMap.put(Networks.readNetwork(n3), 1.0);
        networkMap.put(Networks.readNetwork(n4), 1.0);
        networkMap.put(Networks.readNetwork(n5), 1.0);

        summarize(networkMap);



    }

}
