package edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge;
/*
 * @ClassName:   Merge
 * @Description: Merge two subnetworks
 * @Author:      Zhen Cao
 * @Date:        7/31/23 2:03 PM
 */

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.Pipeline;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.SuperNetwork3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import sun.nio.ch.Net;

import static edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.util.Utils.getLeafSetButReti;


public class NJMergeTopology {
    private List<Network> _subnetworks = new ArrayList<>();
    private List<Network> _originalSubNetworks = new ArrayList<>();
    private int _num_taxa;
    private double[][] _matrix = new double[_num_taxa][_num_taxa];
    private List<String> _taxonList = new ArrayList<>();
    private List<Set<String>> _leaves = null;
    private String _outgroup = "";
    public static double _epsilon = 0.001;
    private int _nodeCnt = 0;
    private double _scale = 1;
    private boolean _debug = Utils._debug;


    /* Constructor */
    public NJMergeTopology(List<Network> subnetworks, double[][] matrix, List<String> taxonList) {
        for (Network net: subnetworks){
            Utils.blankSubNetInternalNodeNames(net);
            this._subnetworks.add(net);
            this._originalSubNetworks.add(net.clone());
        }
        this._matrix = matrix;
        this._taxonList = taxonList;
        _leaves = new ArrayList<>();
        _num_taxa = taxonList.size();

        for (Network net : _subnetworks) {
            Set<String> set = new HashSet<String>();
            for (Object o: net.getLeaves()){
                set.add(((NetNode) o).getName());
            }
            _leaves.add(set);
        }
    }

    public void setOutgroup(String outgroup){
        this._outgroup = outgroup;
    }

    /*Join clade A and clade B in both nets
     */
//    public List<Network> joinNodesInBothNets(Network net1, NetNode nodeAinT1, List<String> cladeAList,
//                                             Network net2, NetNode nodeBinT2, List<String> cladeBList,
//                                             boolean test) {
//
//        Set<String> cladeA = new HashSet<>(cladeAList);
//        Set<String> cladeB = new HashSet<>(cladeBList);
//        Set<String> leaves1 = Utils.getLeafSet(net1.getRoot());
//        Set<String> leaves2 = Utils.getLeafSet(net2.getRoot());
//
//        boolean cladeAisT1 = leaves1.equals(cladeA);
//        boolean cladeBisT2 = leaves2.equals(cladeB);
//
//        if (cladeAisT1 && cladeBisT2) {
//            // Join trees 1 and 2 together, i.e.,
//            // tree 1 will equal tree 2
//            if (test) {
//                return Arrays.asList(null, null);
//            }
//            //todo distances
//            NetNode root = new BniNetNode();
//            root.adoptChild(nodeAinT1, NetNode.NO_DISTANCE);
//            root.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
//            net1 = new BniNetwork((BniNetNode) root);
//            net2 = null;
//        } else if (cladeAisT1) {
//            if (test) {
//                return Arrays.asList(null, null);
//            }
//            NetNode root = new BniNetNode();
//            root.adoptChild(nodeAinT1, NetNode.NO_DISTANCE);
//            root.adoptChild(net2.getRoot(), NetNode.NO_DISTANCE);
//            net1 = new BniNetwork((BniNetNode) root);
//            net2 = null;
//        } else if (cladeBisT2) {
//            if (test) {
//                return Arrays.asList(null, null);
//            }
//            NetNode root = new BniNetNode();
//            root.adoptChild(net1.getRoot(), NetNode.NO_DISTANCE);
//            root.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
//            net1 = new BniNetwork((BniNetNode) root);
//            net2 = null;
//        } else {
//            Map<NetNode, NetNode> old2newT1 = new HashMap<>();
//            Map<NetNode, NetNode> old2newT2 = new HashMap<>();
//            NetNode nodeAinT1copy = Utils.copyNodeSubnet(nodeAinT1, old2newT1);
//            NetNode nodeBinT2copy = Utils.copyNodeSubnet(nodeBinT2, old2newT2);
//
//            NetNode parentA = (NetNode) nodeAinT1.getParents().iterator().next();
//            NetNode newParentNodeT1 = new BniNetNode();
//            parentA.removeChild(nodeAinT1);
//            parentA.adoptChild(newParentNodeT1, NetNode.NO_DISTANCE);
//            newParentNodeT1.adoptChild(nodeAinT1, NetNode.NO_DISTANCE);
//            newParentNodeT1.adoptChild(nodeBinT2copy, NetNode.NO_DISTANCE);
//
//            NetNode parentB = (NetNode) nodeBinT2.getParents().iterator().next();
//            NetNode newParentNodeT2 = new BniNetNode();
//            parentB.removeChild(nodeBinT2);
//            parentB.adoptChild(newParentNodeT2, NetNode.NO_DISTANCE);
//            newParentNodeT2.adoptChild(nodeAinT1copy, NetNode.NO_DISTANCE);
//            newParentNodeT2.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
//        }
//
//        return Arrays.asList(net1, net2);
//    }

    List<Network> joinNodesInBothNets(Network net1, NetNode nodeAinT1, Set<String> cladeA,
                                             Network net2, NetNode nodeBinT2, Set<String> cladeB
                                             ) {
        if (nodeAinT1.getName().equals("I17") && nodeBinT2.getName().equals("Z")){
            System.out.println("debug point");

        }
        Network net1copy = net1.clone();

        Set<String> leaves1 = Utils.getLeafSet(net1.getRoot());
        Set<String> leaves2 = Utils.getLeafSet(net2.getRoot());
        Set<String> cladeAset = new HashSet<>(cladeA);
        Set<String> cladeBset = new HashSet<>(cladeB);

        NetNode<NetNodeInfo> newNodeToAddinT1 = new BniNetNode();
//        newNodeToAdd.setName("I"+_nodeCnt);
        if ((cladeA.equals(leaves1) && cladeB.equals(leaves2)) || (cladeA.equals(leaves2) && cladeB.equals(leaves1))) {

            newNodeToAddinT1.adoptChild(net1.getRoot(), NetNode.NO_DISTANCE);
            newNodeToAddinT1.adoptChild(net2.getRoot(), NetNode.NO_DISTANCE);
            newNodeToAddinT1.setData(new NetNodeInfo());
            net1 = new BniNetwork((BniNetNode) newNodeToAddinT1);
            net2 = net1.clone();

        }
        else {
            // Add node B to T1
            // T1: (A,X) -> ((A,B),X)), need to check if B and X have overlapped species
            if (_debug) {
//                System.out.println("Adding node B " + nodeBinT2.getName() + " to T1");
            }
            Map<NetNode, NetNode> old2new = new HashMap<>();
            NetNode<NetNodeInfo> nodeBcopy = Utils.copyNodeSubnet(nodeBinT2, old2new);
            nodeAinT1.setName(nodeAinT1.getName());
            nodeBcopy.setName(nodeBinT2.getName());

            //check if set of species in network - A is overlapped with clade B

            if (leaves1.equals(cladeAset)) {
                newNodeToAddinT1.adoptChild(nodeAinT1, NetNode.NO_DISTANCE);
                newNodeToAddinT1.setData(new NetNodeInfo());
                newNodeToAddinT1.adoptChild(nodeBcopy, NetNode.NO_DISTANCE);
                net1 = new BniNetwork((BniNetNode) newNodeToAddinT1);
            }
            else {

                Set<String> cladeX = new HashSet<>(leaves1);
                cladeX.removeAll(cladeAset);
                if (_debug) {
                    System.out.println("check point before checkOverlap:"+nodeAinT1.getName()+","+nodeBinT2.getName());
                    System.out.println("cladeX: " + cladeX);
                }
                boolean redundantNodeManagable = checkOverlap(cladeX, cladeBset, cladeAset, net1, nodeAinT1, nodeBcopy, newNodeToAddinT1);
                if(!redundantNodeManagable){
                    return Arrays.asList(null, null);
                }


            }

            System.out.println("after adding node B to T1, net1:"+ net1.toString());
            System.out.println("after adding node B to T1, net2:"+ net2.toString());

            // Add node A to T2
            // (B,X) -> ((A,B),X)), need to check if A and X have overlapped species
            NetNode<NetNodeInfo> newNodeToAddinT2 = new BniNetNode();

            old2new = new HashMap<>();
            NetNode<NetNodeInfo> nodeAcopy = Utils.copyNodeSubnet(nodeAinT1, old2new);
            nodeAcopy.setName(nodeAinT1.getName());
            nodeBinT2.setName(nodeBinT2.getName());
            if (leaves2.equals(cladeBset)) {
                newNodeToAddinT2.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
                newNodeToAddinT2.setData(new NetNodeInfo());
                newNodeToAddinT2.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
                net2 = new BniNetwork((BniNetNode) newNodeToAddinT2);
            }
            else {//there could be several networks
                Set<String> cladeX = new HashSet<>(leaves2);
                cladeX.removeAll(cladeBset);
                if (_debug) {
                    System.out.println("check point before checkOverlap:"+nodeAinT1.getName()+","+nodeBinT2.getName());
                    System.out.println("cladeX: " + cladeX);
                }
                boolean redundantNodeManagable  = checkOverlap(cladeX, cladeAset, cladeBset, net2, nodeBinT2, nodeAcopy, newNodeToAddinT2);
                if(!redundantNodeManagable){
                    return Arrays.asList(null, null);
                }
            }
        }
        System.out.println("after adding node A to T2, net1:"+ net1.toString());
        System.out.println("after adding node A to T2, net2:"+ net2.toString());
        System.out.println("finish checkoverlap-----------------");
        return Arrays.asList(net1, net2);
    }

    // get the node that are redundant in the two networks
    //add node B to T, which has A and X, (A,X), and return redundant nodes in X (in T) and B
    //return nodeB's children to node in network
//    public Tuple<NetNode, NetNode> getRedundantMap(NetNode nodeA, NetNode nodeB, Network nodeAnet, Network nodeBnet){
//        Set<String> leaves = Utils.getLeafSet(nodeAnet.getRoot());
//        Set<String> cladeA = Utils.getLeafSet(nodeAnet.findNode(nodeA.getName()));
//        Set<String> cladeB = Utils.getLeafSet(nodeBnet.findNode(nodeB.getName()));
//        Set<String> cladeX = new HashSet<>(leaves);
//        cladeX.removeAll(cladeA);
//        Set<String> intersection = new HashSet<>(cladeX);
//        intersection.retainAll(cladeB);
//        NetNode nodeUnderA = getMRCA(nodeAnet, intersection);
//        NetNode nodeUnderB = getMRCA(nodeBnet, intersection);
//
//
//        return new Tuple<NetNode, NetNode>(nodeUnderA, nodeUnderB);
//    }


    // clade B in net2, add clade A to net2 while removing overlap intersections of clade A and clade X (net2-B)
//    public static boolean checkOverlap(Set<String> cladeX, Set<String> cladeA, Set<String> cladeB, Network net2, NetNode nodeBinT2, NetNode<NetNodeInfo> nodeAcopy, NetNode newNodeToAddinT2){
//        if (!Collections.disjoint(cladeX, cladeA)) {
//            System.out.println("checkOverlap overlap");
//            // node under nodeB to node in network
//
//            Set<String> intersection = new HashSet<>(cladeX);
//            intersection.retainAll(cladeA);
//            Set<String> unionIntersectionB = new HashSet<>(cladeB);
//            unionIntersectionB.addAll(intersection);
//            System.out.println("testjoin intersectionB:"+intersection);
//            System.out.println("testjoin unionIntersectionB:"+unionIntersectionB);
//            NetNode mrcaXB = getMRCA(net2, unionIntersectionB);
//
//            Map<NetNode, NetNode> redundantNodeMap = getRedundantNode(nodeBinT2, nodeAcopy,  net2);
//
//            if(redundantNodeMap.size()>1){
//                return false;
//            }
//            System.out.println("redundantNodeMap size:"+redundantNodeMap.size());
////            Tuple<NetNode, NetNode> redundantTuple = getRedundantNode(nodeBinT2, nodeAcopy,  net2, net1);
//            for (NetNode nodeunderA : redundantNodeMap.keySet()) {
//                System.out.println(nodeunderA.getName());
//            }
//            NetNode parentAX = null;
//            for (NetNode nodeunderA : redundantNodeMap.keySet()) {
//                parentAX =  (NetNode) nodeunderA.getParents().iterator().next();
//                parentAX.removeChild(nodeunderA);
//
//            }
//
//
//            // the child of mrcaXA that has cladeX is the child of nodeBcopy
////            System.out.println(mrcaXB==null);
//            if (mrcaXB == (null)){
//                System.out.println("check point in joinNodesInBothNets");
//                System.out.println(net2.toString());
//                System.out.println("mrcaXB is null");
//            }
//            for(Object o: mrcaXB.getChildren()){
//                NetNode child = (NetNode) o;
//                if (Utils.getLeafSet(child).containsAll(intersection)){
//                    mrcaXB.removeChild(child);
//                    mrcaXB.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
//                    parentAX.adoptChild(child, NetNode.NO_DISTANCE);
//                    break;
//                }
//            }
//            if(mrcaXB.getName().isEmpty()){
//
//            }
//            else {
//                return false;
//            }
//
//        }
//        else{
//            NetNodeInfo ni = new NetNodeInfo();
//            newNodeToAddinT2.setData(ni);
//            newNodeToAddinT2.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
//            newNodeToAddinT2.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
//
//            NetNode<NetNodeInfo> parentB = (NetNode) nodeBinT2.getParents().iterator().next();
//            parentB.removeChild(nodeBinT2);
//            parentB.adoptChild(newNodeToAddinT2, NetNode.NO_DISTANCE);
//        }
//        return true;
//    }


    public static boolean checkOverlap(Set<String> cladeX, Set<String> cladeA, Set<String> cladeB, Network net2, NetNode nodeBinT2, NetNode<NetNodeInfo> nodeAcopy, NetNode newNodeToAddinT2){
        if (!Collections.disjoint(cladeX, cladeA)) {
            System.out.println("checkOverlap overlap");
            // node under nodeB to node in network

            Set<String> intersection = new HashSet<>(cladeX);
            intersection.retainAll(cladeA);
            Set<String> unionIntersectionB = new HashSet<>(cladeB);
            unionIntersectionB.addAll(intersection);
            System.out.println("testjoin intersectionB:"+intersection);
            System.out.println("testjoin unionIntersectionB:"+unionIntersectionB);
//            NetNode mrcaXB = getMRCA(net2, unionIntersectionB);

//            NetNode intersectionMRCA = getMRCA(net2, intersection);
//            if (!Utils.ReticulationInDescendants(intersectionMRCA)){
//                NetNode parentAX =  (NetNode) intersectionMRCA.getParents().iterator().next();
//                parentAX.removeChild(intersectionMRCA);
//            }
//            else if (!Utils.ReticulationInDescendants(nodeAcopy)){
//                NetNode parentAX =  (NetNode) nodeAcopy.getParents().iterator().next();
//                parentAX.removeChild(nodeAcopy);
//            }
//            else{
//                return false;
//            }
            Map<NetNode, NetNode> redundantNodeMap = getRedundantNode(nodeBinT2, nodeAcopy,  net2);

//            if(redundantNodeMap.size()>1){
//                return false;
//            }
//            System.out.println("redundantNodeMap size:"+redundantNodeMap.size());
////            Tuple<NetNode, NetNode> redundantTuple = getRedundantNode(nodeBinT2, nodeAcopy,  net2, net1);
//            for (NetNode nodeunderA : redundantNodeMap.keySet()) {
//                System.out.println(nodeunderA.getName());
//            }
            NetNode parentAX = null;
            for (NetNode nodeunderA : redundantNodeMap.keySet()) {
                parentAX =  (NetNode) nodeunderA.getParents().iterator().next();
                parentAX.removeChild(nodeunderA);

            }


            // the child of mrcaXA that has cladeX is the child of nodeBcopy
//            System.out.println(mrcaXB==null);
//            if (mrcaXB == (null)){
//                System.out.println("check point in joinNodesInBothNets");
//                System.out.println(net2.toString());
//                System.out.println("mrcaXB is null");
//            }
//            for(Object o: mrcaXB.getChildren()){
//                NetNode child = (NetNode) o;
//                if (Utils.getLeafSet(child).containsAll(intersection)){
//                    mrcaXB.removeChild(child);
//                    mrcaXB.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
//                    parentAX.adoptChild(child, NetNode.NO_DISTANCE);
//                    break;
//                }
//            }
//            if(mrcaXB.getName().isEmpty()){
//
//            }
//            else {
//                return false;
//            }

        }
//        else{
            NetNodeInfo ni = new NetNodeInfo();
            newNodeToAddinT2.setData(ni);
            newNodeToAddinT2.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
            newNodeToAddinT2.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);

            NetNode<NetNodeInfo> parentB = (NetNode) nodeBinT2.getParents().iterator().next();
            parentB.removeChild(nodeBinT2);
            parentB.adoptChild(newNodeToAddinT2, NetNode.NO_DISTANCE);
//        }
        return true;
    }




//    public List<Network> joinNodesInBothNets(Network net1, NetNode nodeAinT1, Set<String> cladeA,
//                                             Network net2, NetNode nodeBinT2, Set<String> cladeB,
//                                             boolean test) {
//
//
//        Set<String> leaves1 = Utils.getLeafSet(net1.getRoot());
//        Set<String> leaves2 = Utils.getLeafSet(net2.getRoot());
//
//        boolean cladeAisT1 = leaves1.equals(cladeA);
//        boolean cladeBisT2 = leaves2.equals(cladeB);
//
//        if (cladeAisT1 && cladeBisT2) {
//            // Join trees 1 and 2 together, i.e.,
//            // tree 1 will equal tree 2
//            if (test) {
//                return Arrays.asList(null, null);
//            }
//            //todo distances
//            NetNode root = new BniNetNode();
//            root.adoptChild(nodeAinT1, NetNode.NO_DISTANCE);
//            root.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
//            net1 = new BniNetwork((BniNetNode) root);
//            net2 = null;
//        } else if (cladeAisT1) {
//            if (test) {
//                return Arrays.asList(null, null);
//            }
//            NetNode root = new BniNetNode();
//            root.adoptChild(nodeAinT1, NetNode.NO_DISTANCE);
//            root.adoptChild(net2.getRoot(), NetNode.NO_DISTANCE);
//            net1 = new BniNetwork((BniNetNode) root);
//            net2 = null;
//        } else if (cladeBisT2) {
//            if (test) {
//                return Arrays.asList(null, null);
//            }
//            NetNode root = new BniNetNode();
//            root.adoptChild(net1.getRoot(), NetNode.NO_DISTANCE);
//            root.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
//            net1 = new BniNetwork((BniNetNode) root);
//            net2 = null;
//        } else {
//            Map<NetNode, NetNode> old2newT1 = new HashMap<>();
//            Map<NetNode, NetNode> old2newT2 = new HashMap<>();
//            NetNode nodeAinT1copy = Utils.copyNodeSubnet(nodeAinT1, old2newT1);
//            NetNode nodeBinT2copy = Utils.copyNodeSubnet(nodeBinT2, old2newT2);
//
//            NetNode parentA = (NetNode) nodeAinT1.getParents().iterator().next();
//            NetNode newParentNodeT1 = new BniNetNode();
//            parentA.removeChild(nodeAinT1);
//            parentA.adoptChild(newParentNodeT1, NetNode.NO_DISTANCE);
//            newParentNodeT1.adoptChild(nodeAinT1, NetNode.NO_DISTANCE);
//            newParentNodeT1.adoptChild(nodeBinT2copy, NetNode.NO_DISTANCE);
//
//            NetNode parentB = (NetNode) nodeBinT2.getParents().iterator().next();
//            NetNode newParentNodeT2 = new BniNetNode();
//            parentB.removeChild(nodeBinT2);
//            parentB.adoptChild(newParentNodeT2, NetNode.NO_DISTANCE);
//            newParentNodeT2.adoptChild(nodeAinT1copy, NetNode.NO_DISTANCE);
//            newParentNodeT2.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
//        }
//
//        return Arrays.asList(net1, net2);
//    }

    public static void getNodeLeafMap(NetNode node, Map<NetNode, Set<String>> nodeLeafMap, Set<String> restrictedSet){
        Set<String> leafset = new HashSet<>();
        if (node.isLeaf()){
            if (restrictedSet == null){
                leafset.add(node.getName());
            }
            else if (restrictedSet.contains(node.getName())){
                leafset.add(node.getName());
            }

        }
        else{
            for (Object o: node.getChildren()){
                NetNode child = (NetNode) o;
                if (!nodeLeafMap.containsKey(child)){
                    getNodeLeafMap(child, nodeLeafMap, restrictedSet);
                }

                leafset.addAll(nodeLeafMap.get(child));

            }
        }
        nodeLeafMap.put(node, leafset);
    }



    //todo test this function
    // net2 has more reticulations than net1, add the reticulations that are only in net1 to net2
    public static Map<NetNode, NetNode> mapSubNets(Network net1, Network net2){
        Map<NetNode, NetNode> nodeCur2Net2 = new HashMap<>();
        List<NetNode> curReticulations = new ArrayList<>();

        Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(net1, net2);
        Tuple<Network, Map<NetNode, NetNode>> tuple = Utils.cloneNetwork(net1);
        Network curnet = tuple.Item1;

        for(Object nodeObj : curnet.dfs()) {
            NetNode node = (NetNode) nodeObj;
            if(node.isNetworkNode()) {
                curReticulations.add(node);
            }
        }
        if(Utils._debug){
            System.out.println(curnet.toString());
        }

        List<Tuple<NetNode, NetNode>> removedRetiList = new ArrayList<>();

        getAllBackbonesDfs(0, curnet, curReticulations, nodeCur2Net2, closest.Item2, removedRetiList);

        //add reticulation edges to net2
        //todo debug not the same node
        Map<NetNode, Set<String>> node2leaf1 = new HashMap<>();
        Map<NetNode, Set<String>> node2leaf2 = new HashMap<>();
        getNodeLeafMap(net1.getRoot(), node2leaf1, null);
        getNodeLeafMap(net2.getRoot(), node2leaf2, null);
        if (Utils._debug){
            System.out.println("reticulations");

        }

        for (Tuple<NetNode, NetNode> edge: removedRetiList){
            NetNode head = nodeCur2Net2.get(edge.Item1);
            NetNode tail = nodeCur2Net2.get(edge.Item2);
            if (Utils._debug){
                System.out.println(head+","+ tail);
            }

        }

        return nodeCur2Net2;
    }


    // Add newVnode to edge(nodeVParent, nodeV), newRetic to edge(nodeU, child), and edge(newVnode, newRetic)
    public static void addReticulation(NetNode<NetNodeInfo> nodeV, NetNode<NetNodeInfo> nodeVParent, NetNode<NetNodeInfo> newVnode, NetNode<NetNodeInfo> child, NetNode<NetNodeInfo> nodeU, NetNode<NetNodeInfo> newRetic){
//        NetNode<NetNodeInfo> newVnode = new BniNetNode();
//        NetNode<NetNodeInfo> newRetic = new BniNetNode();
        NetNodeInfo niRetic = new NetNodeInfo();
        newRetic.setData(niRetic);
        NetNodeInfo niV = new NetNodeInfo();
        newVnode.setData(niV);
        newVnode.adoptChild(newRetic, NetNode.NO_DISTANCE);
        newVnode.adoptChild(nodeV, NetNode.NO_DISTANCE);
        nodeVParent.adoptChild(newVnode, NetNode.NO_DISTANCE);
        nodeVParent.removeChild(nodeV);
        nodeU.adoptChild(newRetic, NetNode.NO_DISTANCE);
        newRetic.adoptChild(child, NetNode.NO_DISTANCE);
        nodeU.removeChild(child);
    }

    public static void removeReticulation(NetNode<NetNodeInfo> nodeV, NetNode<NetNodeInfo> nodeVParent, NetNode<NetNodeInfo> newVnode, NetNode<NetNodeInfo> child, NetNode<NetNodeInfo> nodeU, NetNode<NetNodeInfo> newRetic){
        newVnode.removeChild(nodeV);
        newVnode.removeChild(newRetic);
        nodeVParent.removeChild(newVnode);
        nodeVParent.adoptChild(nodeV, NetNode.NO_DISTANCE);
        nodeU.removeChild(newRetic);
        nodeU.adoptChild(child, NetNode.NO_DISTANCE);
        newRetic.removeChild(child);

    }


    public static void mergeReticulations(Set<String> nodeUleaf, Set<String> reticLeafSet, NetNode<NetNodeInfo> child, NetNode<NetNodeInfo> parent, Network net2, Map<Set<String>, List<NetNode>> set2node2, List<Network> candidateNetworks){
        Set<String> nodeULeafButReti = new HashSet<>(nodeUleaf);
        nodeULeafButReti.removeAll(reticLeafSet);
        List<NetNode> nodeUChildlist = new ArrayList<>();
        for (Set<String> temp: set2node2.keySet()) {
            if (temp.equals(nodeULeafButReti)) {
                nodeUChildlist = set2node2.get(temp);
                break;
            }
        }

        if (nodeUChildlist.size() > 0){
            for (NetNode<NetNodeInfo> nodeUChild2: nodeUChildlist){

                List<NetNode> nodeUParent2List = StreamSupport.stream(nodeUChild2.getParents().spliterator(), false)
                        .collect(Collectors.toList());

                for (NetNode<NetNodeInfo> nodeU2: nodeUParent2List) {
                    if (nodeU2.equals(parent)){
                        continue;
                    }

                    NetNode<NetNodeInfo> newVnode = new BniNetNode();
                    NetNode<NetNodeInfo> newRetic = new BniNetNode();
                    addReticulation(nodeUChild2, nodeU2, newVnode, child, parent, newRetic);
                    Network newnet = net2.clone();
                    candidateNetworks.add(newnet);
                    removeReticulation(nodeUChild2, nodeU2, newVnode, child, parent, newRetic);
                }
            }
        }
    }

    //TODO: finish this function
    //match the reticulations in two networks, add reticulations on net1 to net2
    public static List<Network> matchNetworks(Network<NetNodeInfo> net1, Network<NetNodeInfo> net2){

        Networks.removeBinaryNodes(net1);
        Networks.removeBinaryNodes(net2);
        System.out.println(net1);
        System.out.println(net2);
        if (Utils._debug){
            Networks.autoLabelNodes(net1);
            Networks.autoLabelNodes(net2);
        }

//        Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(net1, net2);
//        Network net1copy = Networks.readNetwork(net1.toString());
//        Network net2copy = Networks.readNetwork(net2.toString());
//        Networks.autoLabelNodes(net1copy);
//        Networks.autoLabelNodes(net2copy);
//        if(Utils._debug){
//            System.out.println(net1copy.toString());
//            System.out.println(net2copy.toString());
//        }

        Map<NetNode, Set<String>> node2leaf1 = new HashMap<>();
        Map<NetNode, Set<String>> node2leaf2 = new HashMap<>();

        List<Network> candidateNetworks = new ArrayList<>();

        Set<String> intersection = Utils.getLeafSet(net1.getRoot());
        intersection.retainAll(Utils.getLeafSet(net2.getRoot()));

        getNodeLeafMap(net1.getRoot(), node2leaf1, intersection);
        getNodeLeafMap(net2.getRoot(), node2leaf2, intersection);
        Network subnet1 = SuperNetwork3.getSubNetwork(Networks.readNetwork(net1.toString()), new ArrayList<>(intersection), true).Item1;
        Network subnet2 = SuperNetwork3.getSubNetwork(Networks.readNetwork(net2.toString()), new ArrayList<>(intersection), true).Item1;

        Map<Set<String>, List<NetNode>> set2node2 = Utils.invertMapUsingGroupingBy(node2leaf2);

        if(net1.getReticulationCount() == 0 || Networks.hasTheSameTopology(subnet1, subnet2)){
            candidateNetworks.add(net2);
        }
        else{
            for (NetNode<NetNodeInfo> retic: net1.getNetworkNodes()){
                Set<String> reticLeafSet = node2leaf1.get(retic);
                Iterator<NetNode<NetNodeInfo>> parentIter = retic.getParents().iterator();
                NetNode<NetNodeInfo> nodeU = parentIter.next();
                NetNode<NetNodeInfo> nodeV = parentIter.next();
                Set<String> nodeUleaf = node2leaf1.get(nodeU);
                Set<String> nodeVleaf = node2leaf1.get(nodeV);
                NetNode<NetNodeInfo> nodeUParent = nodeU.getParents().iterator().next();
                NetNode<NetNodeInfo> nodeVParent = nodeV.getParents().iterator().next();
                Set<String> nodeVParentLeaf = node2leaf1.get(nodeVParent);
                Set<String> nodeUParentLeaf = node2leaf1.get(nodeUParent);

                //Todo: debug here
                if (set2node2.containsKey(reticLeafSet)){
                    List<NetNode> node2list= set2node2.get(reticLeafSet);
                    for (NetNode<NetNodeInfo> child: node2list){
                        List<NetNode> parents = StreamSupport.stream(child.getParents().spliterator(), false)
                                .collect(Collectors.toList());
                        for (Object op: parents){
                            NetNode<NetNodeInfo> parent = (NetNode<NetNodeInfo>) op;
                            if (node2leaf2.get(parent).equals(nodeUleaf)){
                                mergeReticulations(nodeVleaf, reticLeafSet, child, parent, net2, set2node2, candidateNetworks);

                            }

                            else if (node2leaf2.get(parent).equals(nodeVleaf)){
                                //todo: debug here, shouldn't be parent of parent of reticulation
                                mergeReticulations(nodeUleaf,  reticLeafSet, child,  parent, net2, set2node2,  candidateNetworks);
                            }
                        }
                    }
                }
            }
        }
        if (Utils._debug){
            System.out.println("candidate networks after matching two networks are below");
            for (Network net: candidateNetworks){
                System.out.println(net.toString());
            }
        }

        return candidateNetworks;
    }


    public static void getAllBackbonesDfs(int index, Network cur, List<NetNode> curReticulations, Map<NetNode, NetNode> nodeCur2Net2, Network net2, List<Tuple<NetNode, NetNode>> removedRetiList) {
        if(index >= curReticulations.size()) return;

        NetNode curnode = curReticulations.get(index);
        List<NetNode> parents = new ArrayList<>();
        for(Object parentObj : curnode.getParents()) {
            NetNode parent = (NetNode) parentObj;
            parents.add(parent);
        }

        for(NetNode parent : parents) {
            double distanceBackup = curnode.getParentDistance(parent);
            parent.removeChild(curnode);

            Tuple<Network, Map<NetNode, NetNode>> tuple = Utils.cloneNetwork(cur);
            Network temp = tuple.Item1;
            Networks.removeBinaryNodes(temp);
            if (Networks.hasTheSameTopology(temp, net2)){
//                Networks.removeBinaryNodes(cur);
                Map<NetNode, NetNode> temp2net2 = Networks.mapTwoNetworks(temp, net2);
                for (NetNode newnode: tuple.Item2.keySet()){
                    nodeCur2Net2.put(tuple.Item2.get(newnode), temp2net2.get(newnode));

                }
//                nodemaps.putAll();
                removedRetiList.add(new Tuple<>(parent, curnode));
                return;
            }

            getAllBackbonesDfs(index + 1, cur, curReticulations, nodeCur2Net2, net2, removedRetiList);

            parent.adoptChild(curnode, distanceBackup);
        }
        getAllBackbonesDfs(index + 1, cur, curReticulations, nodeCur2Net2, net2, removedRetiList);
        return;
    }




    // pair subnetworks together
    public Tuple<Network, Map<NetNode, NetNode>> pairNets(Network net1, Network net2){
        Set<String> intersection = Utils.getLeafSet(net1.getRoot());
        intersection.retainAll(Utils.getLeafSet(net2.getRoot()));
        Tuple<Network, Map<NetNode, NetNode>> subnet1 = SuperNetwork3.getSubNetwork(net1, new ArrayList<>(intersection), true);
        Tuple<Network, Map<NetNode, NetNode>> subnet2 = SuperNetwork3.getSubNetwork(net2, new ArrayList<>(intersection), true);


        if (intersection.size() > 1){

            if (!Networks.hasTheSameTopology(subnet1.Item1, subnet2.Item1)){
                Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(subnet1.Item1, subnet2.Item1);
                if (closest.Item3 > 0.001){
                    return null;
                }
                if (net1.getReticulationCount() > net2.getReticulationCount()){
                    //TODO add reticulation to the subnet2, subnet2 should equal closest
                    Network net1copy = Networks.readNetwork(net1.toString());
                    for(Object o: net1.getNetworkNodes()){
                        NetNode netnode = (NetNode) o;
                        List<NetNode> parents = IterableHelp.toList(netnode.getParents());

                    }
                }
                else{
                    //TODO add reticulation to the subnet1, subnet1 should equal closest
                }


            }

        }
        return null;

    }



    // redundantNodeMap is node to node in nodesInNet
//    public static void redundantChecker(NetNode node, Map<String, NetNode> nodesInNet, Map<NetNode, NetNode> redundantNodeMap){
//        if (nodesInNet.containsKey(node.getName())){
//            redundantNodeMap.put(node, nodesInNet.get(node.getName()));
//        }
//
//        for (Object o: node.getChildren()){
//            NetNode child = (NetNode) o;
//            if (!child.isNetworkNode() && (!child.getName().equals(""))){
//                redundantChecker(child, nodesInNet, redundantNodeMap);
//            }
//        }
//    }
//
//    public static void getAllNodesUnderNode(NetNode node, List<NetNode> childList){
//        if (!node.isNetworkNode() && (!node.getName().equals(""))){
//            childList.add(node);
//        }
//        for (Object o: node.getChildren()){
//            NetNode child = (NetNode) o;
//            getAllNodesUnderNode(child, childList);
//        }
//    }


//        public static void getNodeUnderRoot(NetNode root, Set<String> clade){
//        if (nodesInNet.containsKey(node.getName())){
//            redundantNodeMap.put(node, nodesInNet.get(node.getName()));
//        }
//
//        for (Object o: node.getChildren()){
//            NetNode child = (NetNode) o;
//            if (!child.isNetworkNode() && (!child.getName().equals(""))){
//                redundantChecker(child, nodesInNet, redundantNodeMap);
//            }
//        }
//    }
    public static void getAllNodesUnderNode(NetNode node, List<NetNode> childList){
//        if (!node.isNetworkNode() && (!node.getName().equals(""))){
//            childList.add(node);
//        }
//        if (!node.isNetworkNode() && (!node.getName().isEmpty())){
            childList.add(node);
//        }
        for (Object o: node.getChildren()){
            NetNode child = (NetNode) o;
            getAllNodesUnderNode(child, childList);
        }
    }


    // redundantNodeMap is node to node in nodesInNet
    public static void redundantChecker(NetNode node, Map<String, NetNode> nodesInNet, Map<NetNode, NetNode> redundantNodeMap){
//        boolean isAdded = false;
        if (nodesInNet.containsKey(node.getName())){
            redundantNodeMap.put(node, nodesInNet.get(node.getName()));
//            isAdded = true;
            return;
        }

        for (Object o: node.getChildren()){
            NetNode child = (NetNode) o;
            redundantChecker(child, nodesInNet, redundantNodeMap);

        }
    }


    //add node B to T, which has A and X, (A,X), and return redundant nodes in X (in T) and B
    //return nodeB's children to node in network
    public static Map<NetNode, NetNode> getRedundantNode(NetNode nodeA, NetNode nodeB, Network network){
        List<NetNode> nodeInNetwork = new ArrayList<>();
        List<NetNode> nodeUnderNodeA = new ArrayList<>();
        Map<NetNode, NetNode> redundantMap = new HashMap<>();

        getAllNodesUnderNode(network.getRoot(), nodeInNetwork);
        getAllNodesUnderNode(nodeA, nodeUnderNodeA);


        Set<NetNode> nodeSetInNetwork = new HashSet<>(nodeInNetwork);
        nodeSetInNetwork.removeAll(nodeUnderNodeA);

        Map<String, NetNode> names2NodeInNetwork = new HashMap<>();
        for (NetNode n: nodeSetInNetwork){
            if (!n.getName().isEmpty()){
                names2NodeInNetwork.put(n.getName(), n);
            }

        }

        //nodeB's children to node in network
        redundantChecker(nodeB, names2NodeInNetwork, redundantMap);

        return redundantMap;

    }

    public static void makeDuplicateLeavesReticulation(Network network){
        Set<String> leaves = Utils.getLeafSet(network.getRoot());
        Map<String, Integer> cnt = new HashMap<>();
        for (Object o: network.getLeaves()){
            NetNode leaf = (NetNode) o;
            if (cnt.containsKey(leaf.getName())){
                cnt.put(leaf.getName(), cnt.get(leaf.getName())+1);
            }
            else{
                cnt.put(leaf.getName(), 1);
            }
        }
    }


    //Join cladeA and cladeB in one or both trees
//    public Tuple<List<Boolean>, NetNode<NetNodeInfo>> joinNodes(NetNode<NetNodeInfo> nodeA, NetNode<NetNodeInfo> nodeB) {
////        Set<String> cladeA = Utils.getLeafSetButReti(nodeA);
////        Set<String> cladeB = Utils.getLeafSetButReti(nodeB);
//        Set<String> cladeA = Utils.getLeafSet(nodeA);
//        Set<String> cladeB = Utils.getLeafSet(nodeB);
//        Set<String> cladeAB = new HashSet<>(cladeA);
//        cladeAB.addAll(cladeB);
//
////        cladeA = getLeafSetButReti(nodeAinTs.get(0));
////        cladeB = getLeafSetButReti(nodeBinTs.get(0));
////        if (!Collections.disjoint(cladeA, cladeB)) {
////            throw new Exception("Nodes are not disjoint on their leaf sets! "+nodeA.getName()+":"+cladeA+";"+nodeB.getName()+":"+cladeB);
////        }
//
//        List<NetNode> nodeAinTs = new ArrayList<>();
//        List<NetNode> nodeBinTs = new ArrayList<>();
//
//        for (int i = 0; i < _subnetworks.size(); i++) {
//            Set<String> leaf = this._leaves.get(i);
//
//            if (cladeA.equals(leaf)) {
//                nodeAinTs.add(nodeA);
//            } else {
//                NetNode temp = getNodeFromClade(_subnetworks.get(i), cladeA);
//                nodeAinTs.add(temp);
//            }
//
//            if (cladeB.equals(leaf)) {
//                nodeBinTs.add(nodeB);
//            } else {
//                NetNode temp = getNodeFromClade(_subnetworks.get(i), cladeB);
//                nodeBinTs.add(temp);
//            }
//        }
//
//        List<Boolean> nAinTs = new ArrayList<>();
//        for (NetNode nodeAinT : nodeAinTs) {
//            nAinTs.add(nodeAinT != null);
//        }
//
//        List<Boolean> nBinTs = new ArrayList<>();
//        for (NetNode nodeBinT : nodeBinTs) {
//            nBinTs.add(nodeBinT != null);
//        }
//
//        int i =0;
//        boolean nAinT1 = nAinTs.get(i);
//        boolean nBinT1 = nBinTs.get(i);
//        NetNode nodeAinT1 = nodeAinTs.get(i);
//        NetNode nodeBinT1 = nodeBinTs.get(i);
//
//        int j = 1;
//        boolean nAinT2 = nAinTs.get(j);
//        boolean nBinT2 = nBinTs.get(j);
//        NetNode nodeAinT2 = nodeAinTs.get(j);
//        NetNode nodeBinT2 = nodeBinTs.get(j);
//        Network n1 = _subnetworks.get(i).clone();
//        Network n2 = _subnetworks.get(j).clone();
//
//
//
//        List<Boolean> edits = Stream.generate(() -> false).limit(_subnetworks.size()).collect(Collectors.toList());
//
//
//        NetNode<NetNodeInfo> newNodeToAdd = new BniNetNode();
//        newNodeToAdd.setName("I"+_nodeCnt);
//
//        if ((nAinT1 || nAinT2) && (nBinT1 || nBinT2)) {
//            if (nAinT1 && nAinT2) {
//                //nodeA in *both* T1 and T2
//                if (nBinT1 && nBinT2) {
//                    // Case 1: nodeB in *both* T1 and T2
//                    // Valid if nodeA and nodeB are siblings in both T1 & T2
//                    NetNode node1 = getNodeFromClade(_subnetworks.get(i), cladeAB);
//                    NetNode node2 = getNodeFromClade(_subnetworks.get(j), cladeAB);
//                    nodeAinT.setName(nodeA.getName());
//                    nodeBinT.setName(nodeB.getName());
//                    Set<String> cladeUnion = new HashSet<>();
//                    cladeUnion.addAll(cladeA);
//                    cladeUnion.addAll(cladeB);
//                    newNodeToAdd = getNodeFromClade(_subnetworks.get(i), cladeUnion);
//                    newNodeToAdd.setName("I"+_nodeCnt);
//                } else if (nBinT1) {
//                    //Case 2: node B in T1 only
//                    // Valid if nodeA and nodeB are siblings in T1
//                    NetNode node = getNodeFromClade(_subnetworks.get(i), cladeAB);
//                    if (node == null) {
//                        violates = true;
//                    }
//                } else if (nBinT2) {
//                    // Case 3: Node B in T2 only
//                    // Valid if nodeA and nodeB are siblings in T2
//                    NetNode node = getNodeFromClade(_subnetworks.get(j), cladeAB);
//                    if (node == null) {
//                        violates = true;
//                    }
//                } else {
//                    throw new Exception("Node B not found in either tree!");
//                }
//            } else if (nAinT1) {
//                // nodeA in T1 only
//                if(nBinT1 && nBinT2){
//                    // Case 4: nodeB in *both* T1 and T2
//                    NetNode node = getNodeFromClade(_subnetworks.get(i), cladeAB);
//                    if (node == null) {
//                        violates = true;
//                    }
//                }
//                else if (nBinT1){
//                    // Case 5: Node B in T1 only
//                    // Valid if nodeA and nodeB are siblings in T1
//                    NetNode node = getNodeFromClade(_subnetworks.get(i), cladeAB);
//                    if (node == null){
//                        violates = true;
//                    }
//                }
//                else if (nBinT2){
//                    // Case 6: Node B in T2 only
//                    // Do join in both trees and test for compatibility
////                            Network n1 = _subnetworks.get(i).clone();
////                            Network n2 = _subnetworks.get(j).clone();
//                    NetNode nA = getNodeFromClade(n1, cladeA);
//                    NetNode nB = getNodeFromClade(n2, cladeB);
//
//                    List<Network> networklist = joinNodesInBothNets(n1, nA, new ArrayList<>(cladeA), n2, nB, new ArrayList<>(cladeB));
//                    n1 = networklist.get(0);
//                    n2 = networklist.get(1);
//                    if (n1 !=  null){
//                        violates = (!areNetsCompatible(n1, n2));//TODO
//                    }
//                }
//               
//
//            } else if (nAinT2) {
//                //nodeA in T2 only
//                if (nBinT1 && nBinT2){
//                    //Case 7: nodeB in and T2
//                    NetNode node = getNodeFromClade(_subnetworks.get(j), cladeAB);
//                    if (node ==  null){
//                        violates = true;
//                    }
//
//                }
//                else if (nBinT1){
//                    //Case 8 (reverse of Case 6): Node B in T1 only
////                    Network n1 = _subnetworks.get(i).clone();
////                    Network n2 = _subnetworks.get(j).clone();
//
//                    NetNode nB = getNodeFromClade(n1, cladeB);
//                    NetNode nA = getNodeFromClade(n2, cladeA);
//                    List<Network> networklist = joinNodesInBothNets(n1, nB, new ArrayList<>(cladeB), n2, nA, new ArrayList<>(cladeA));
//                    n1 = networklist.get(0);
//                    n2 = networklist.get(1);
//                    if (n1 !=  null){
//                        violates = (!areNetsCompatible(n1, n2));//TODO
//                    }
//                }
//                else if (nBinT2){
//                    // Case 9: Node B in T2 only
//                    // Only valid if (nodeA, nodeB) are siblings in T2
//                    NetNode node = getNodeFromClade(_subnetworks.get(j), cladeAB);
//                    if (node ==  null){
//                        violates = true;
//                    }
//                }
//
//            }
//        }
//
//
//
//        if ((cladeA.equals(_leaves.get(0)) && cladeB.equals(_leaves.get(1))) || (cladeA.equals(_leaves.get(1)) && cladeB.equals(_leaves.get(0)))){
//            edits.set(0, true);
//            edits.set(1, true);
//
//            newNodeToAdd.adoptChild(_subnetworks.get(0).getRoot(), NetNode.NO_DISTANCE);
//            newNodeToAdd.adoptChild(_subnetworks.get(1).getRoot(), NetNode.NO_DISTANCE);
//            newNodeToAdd.setData(new NetNodeInfo());
//            _subnetworks.set(0, new BniNetwork((BniNetNode) newNodeToAdd));
//            _subnetworks.set(1, _subnetworks.get(0).clone());
//
//        }
//        else{
//            for (int i = 0; i < _subnetworks.size(); i++) {
//                newNodeToAdd = new BniNetNode();
//                newNodeToAdd.setName("I"+_nodeCnt);
//                Set<String> leaf = _leaves.get(i);
//                NetNode<NetNodeInfo> nodeAinT;
//                NetNode<NetNodeInfo> nodeBinT;
//
//                if (cladeA.equals(leaf)) {
//                    nodeAinT = _subnetworks.get(i).getRoot();
//
//                } else {
//                    nodeAinT = getNodeFromClade(_subnetworks.get(i), cladeA);
//                }
//
//                if (cladeB.equals(leaf)) {
//                    nodeBinT = _subnetworks.get(i).getRoot();
//                } else {
//                    nodeBinT = getNodeFromClade(_subnetworks.get(i), cladeB);
//                }
//
//                boolean nAinT = nodeAinT != null;
//                boolean nBinT = nodeBinT != null;
//
//                if (nAinT && nBinT) {
//                    // NetNode A and node B are both in T, do nothing!
//                    nodeAinT.setName(nodeA.getName());
//                    nodeBinT.setName(nodeB.getName());
//                    Set<String> cladeUnion = new HashSet<>();
//                    cladeUnion.addAll(cladeA);
//                    cladeUnion.addAll(cladeB);
//                    newNodeToAdd = getNodeFromClade(_subnetworks.get(i), cladeUnion);
//                    newNodeToAdd.setName("I"+_nodeCnt);
//                    break;
//
//                } else if (nAinT) {
//                    // Add node B to T
//                    // T: (A,X) -> ((A,B),X)), need to check if B and X have overlapped species
//                    if (_debug) {
//                        System.out.println("Adding node B "+nodeB.getName()+" to T"+i);
//                    }
//                    edits.set(i, true);
//                    Map<NetNode, NetNode> old2new = new HashMap<>();
//                    NetNode<NetNodeInfo> nodeBcopy = Utils.copyNodeSubnet(nodeB, old2new);
//                    nodeAinT.setName(nodeA.getName());
//                    nodeBcopy.setName(nodeB.getName());
//
//                    //check if set of species in network - A is overlapped with clade B
//                    if (leaf.equals(cladeA)) {
//                        newNodeToAdd.adoptChild(nodeAinT, NetNode.NO_DISTANCE);
//                        newNodeToAdd.setData(new NetNodeInfo());
//                        newNodeToAdd.adoptChild(nodeBcopy, NetNode.NO_DISTANCE);
//                        _subnetworks.set(i, new BniNetwork((BniNetNode) newNodeToAdd));
//                    }
//                    else {
//
//                        Set<String> cladeX = new HashSet<>(leaf);
//                        cladeX.removeAll(cladeA);
//
//                        if (!Collections.disjoint(cladeX, cladeB)) {
//                            // X and clade B overlaps
//                            // node under nodeB to node in network
//                            NetNode parentBX = null;
//                            Map<NetNode, NetNode> redundantNodeMap = getRedundantNode(nodeAinT, nodeBcopy, _subnetworks.get(i));
//                            for (NetNode nodeunderB : redundantNodeMap.keySet()) {
//                                parentBX =  (NetNode) nodeunderB.getParents().iterator().next();
//                                parentBX.removeChild(nodeunderB);
//
//                            }
//                            // use mrca of intersection of (X,B) and A as the parent of nodeBcopy
//                            Set<String> intersection = new HashSet<>(cladeX);
//                            intersection.retainAll(cladeB);
//                            intersection.addAll(cladeA);
//                            System.out.println("intersection:"+intersection);
//                            NetNode mrcaXA = getNodeFromClade(_subnetworks.get(i), intersection);
//                            // the child of mrcaXA that has cladeX is the child of nodeBcopy
//                            for(Object o: mrcaXA.getChildren()){
//                                NetNode child = (NetNode) o;
//                                if (Utils.getLeafSet(child).containsAll(cladeX)){
//                                    mrcaXA.removeChild(child);
//                                    mrcaXA.adoptChild(nodeBcopy, NetNode.NO_DISTANCE);
//                                    parentBX.adoptChild(child, NetNode.NO_DISTANCE);
//                                    break;
//                                }
//                            }
//                       }
//                        else{
//                            newNodeToAdd.adoptChild(nodeAinT, NetNode.NO_DISTANCE);
//                            newNodeToAdd.setData(new NetNodeInfo());
//                            newNodeToAdd.adoptChild(nodeBcopy, NetNode.NO_DISTANCE);
//                            NetNode parentA = nodeAinT.getParents().iterator().next();
//                            parentA.removeChild(nodeAinT);
//                            parentA.adoptChild(newNodeToAdd, NetNode.NO_DISTANCE);
//                        }
//                    }
//
//
//                } else if (nBinT) {
//                    // Add node A to T
//                    // (B,X) -> ((A,B),X)), need to check if A and X have overlapped species
//                    if (_debug) {
//                        System.out.println("Adding node A "+nodeA.getName()+" to T"+i);
//                    }
//                    edits.set(i, true);
//                    Map<NetNode, NetNode> old2new = new HashMap<>();
//                    NetNode<NetNodeInfo> nodeAcopy = Utils.copyNodeSubnet(nodeA, old2new);
//                    nodeAcopy.setName(nodeA.getName());
//                    nodeBinT.setName(nodeB.getName());
//                    if (leaf.equals(cladeB)) {
//                        newNodeToAdd.adoptChild(nodeBinT, NetNode.NO_DISTANCE);
//                        newNodeToAdd.setData(new NetNodeInfo());
//                        newNodeToAdd.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
//                        _subnetworks.set(i, new BniNetwork((BniNetNode) newNodeToAdd));
//                    } else {//there could be several networks
//                        Set<String> cladeX = new HashSet<>(leaf);
//                        cladeX.removeAll(cladeB);
//
//                        if (!Collections.disjoint(cladeX, cladeA)) {
//                          // node under nodeB to node in network
//                            NetNode parentAX = null;
//                            Map<NetNode, NetNode> redundantNodeMap = getRedundantNode(nodeBinT, nodeAcopy, _subnetworks.get(i));
//                            for (NetNode nodeunderA : redundantNodeMap.keySet()) {
//                                parentAX =  (NetNode) nodeunderA.getParents().iterator().next();
//                                parentAX.removeChild(nodeunderA);
//
//                            }
//                            Set<String> intersection = new HashSet<>(cladeX);
//                            intersection.retainAll(cladeA);
//                            intersection.addAll(cladeB);
//                            NetNode mrcaXB = getNodeFromClade(_subnetworks.get(i), intersection);
//                            // the child of mrcaXA that has cladeX is the child of nodeBcopy
//
//                            for(Object o: mrcaXB.getChildren()){
//                                NetNode child = (NetNode) o;
//                                if (Utils.getLeafSet(child).containsAll(cladeX)){
//                                    mrcaXB.removeChild(child);
//                                    mrcaXB.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
//                                    mrcaXB.adoptChild(child, NetNode.NO_DISTANCE);
//                                    break;
//                                }
//                            }
//
//                        }
//                        else{
//                            NetNodeInfo ni  = new NetNodeInfo();
//                            newNodeToAdd.setData(ni);
//                            NetNode<NetNodeInfo> parentB = nodeBinT.getParents().iterator().next();
//                            parentB.removeChild(nodeBinT);
//                            parentB.adoptChild(newNodeToAdd, NetNode.NO_DISTANCE);
//                            newNodeToAdd.adoptChild(nodeBinT, NetNode.NO_DISTANCE);
//                            newNodeToAdd.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
//                        }
//
//                    }
//                }
//            }
//        }
//        Map<NetNode, NetNode> old2new = new HashMap<>();
//        NetNode newNode = Utils.copyNodeSubnet(newNodeToAdd, old2new);
//        return new Tuple<>(edits, newNode);
//    }


    // this is the original joinnodes
    public Tuple<List<Boolean>, NetNode<NetNodeInfo>> joinNodes(NetNode<NetNodeInfo> nodeA, NetNode<NetNodeInfo> nodeB) {
//        Set<String> cladeA = Utils.getLeafSetButReti(nodeA);
//        Set<String> cladeB = Utils.getLeafSetButReti(nodeB);
        Set<String> cladeA = Utils.getLeafSet(nodeA);
        Set<String> cladeB = Utils.getLeafSet(nodeB);
        if (nodeA.getName().equals("Z")  && nodeB.getName().equals("I8")){
            System.out.println("checking point");
        }



        List<Boolean> edits = Stream.generate(() -> false).limit(_subnetworks.size()).collect(Collectors.toList());



        NetNode<NetNodeInfo> newNodeToAdd = new BniNetNode();
        newNodeToAdd.setName("I"+_nodeCnt);

        if ((cladeA.equals(_leaves.get(0)) && cladeB.equals(_leaves.get(1))) || (cladeA.equals(_leaves.get(1)) && cladeB.equals(_leaves.get(0)))){
            edits.set(0, true);
            edits.set(1, true);

            newNodeToAdd.adoptChild(_subnetworks.get(0).getRoot(), NetNode.NO_DISTANCE);
            newNodeToAdd.adoptChild(_subnetworks.get(1).getRoot(), NetNode.NO_DISTANCE);
            newNodeToAdd.setData(new NetNodeInfo());
            _subnetworks.set(0, new BniNetwork((BniNetNode) newNodeToAdd));
            _subnetworks.set(1, _subnetworks.get(0).clone());

        }
        else{
            for (int i = 0; i < _subnetworks.size(); i++) {
                newNodeToAdd = new BniNetNode();
                newNodeToAdd.setName("I"+_nodeCnt);
                Set<String> leaf = _leaves.get(i);
                NetNode<NetNodeInfo> nodeAinT;
                NetNode<NetNodeInfo> nodeBinT;

                if (cladeA.equals(leaf)) {
                    nodeAinT = _subnetworks.get(i).getRoot();

                } else {
                    nodeAinT = getNodeFromClade(_subnetworks.get(i), cladeA);
                }

                if (cladeB.equals(leaf)) {
                    nodeBinT = _subnetworks.get(i).getRoot();
                } else {
                    nodeBinT = getNodeFromClade(_subnetworks.get(i), cladeB);
                }

                boolean nAinT = nodeAinT != null;
                boolean nBinT = nodeBinT != null;

                if (nAinT && nBinT) {
                    // NetNode A and node B are both in T, do nothing!
                    nodeAinT.setName(nodeA.getName());
                    nodeBinT.setName(nodeB.getName());
                    Set<String> cladeUnion = new HashSet<>();
                    cladeUnion.addAll(cladeA);
                    cladeUnion.addAll(cladeB);
                    newNodeToAdd = getMRCA(_subnetworks.get(i), cladeUnion);
                    newNodeToAdd.setName("I"+_nodeCnt);
                    break;

                } else if (nAinT) {
                    // Add node B to T
                    // T: (A,X) -> ((A,B),X)), need to check if B and X have overlapped species
                    if (_debug) {
                        System.out.println("Adding node B "+nodeB.getName()+" to T"+i);
                    }
                    edits.set(i, true);
                    Map<NetNode, NetNode> old2new = new HashMap<>();
                    NetNode<NetNodeInfo> nodeBcopy = Utils.copyNodeSubnet(nodeB, old2new);
                    nodeAinT.setName(nodeA.getName());
                    nodeBcopy.setName(nodeB.getName());

                    //check if set of species in network - A is overlapped with clade B
                    if (leaf.equals(cladeA)) {
                        newNodeToAdd.adoptChild(nodeAinT, NetNode.NO_DISTANCE);
                        newNodeToAdd.setData(new NetNodeInfo());
                        newNodeToAdd.adoptChild(nodeBcopy, NetNode.NO_DISTANCE);
                        _subnetworks.set(i, new BniNetwork((BniNetNode) newNodeToAdd));
                    }
                    else {

                        Set<String> cladeX = new HashSet<>(leaf);
                        cladeX.removeAll(cladeA);
                        checkOverlap(cladeX, cladeB, cladeA, _subnetworks.get(i), nodeAinT, nodeBcopy, newNodeToAdd);

//                        if (!Collections.disjoint(cladeX, cladeB)) {
//                            // X and clade B overlaps
//                            // node under nodeB to node in network
//                            NetNode parentBX = null;
//                            Map<NetNode, NetNode> redundantNodeMap = getRedundantNode(nodeAinT, nodeBcopy, _subnetworks.get(i));
//                            for (NetNode nodeunderB : redundantNodeMap.keySet()) {
//                                parentBX =  (NetNode) nodeunderB.getParents().iterator().next();
//                                parentBX.removeChild(nodeunderB);
//
//                            }
//                            // use mrca of intersection of (X,B) and A as the parent of nodeBcopy
//                            Set<String> intersection = new HashSet<>(cladeX);
//                            intersection.retainAll(cladeB);
//                            intersection.addAll(cladeA);
//                            System.out.println("intersection:"+intersection);
//                            NetNode mrcaXA = getNodeFromClade(_subnetworks.get(i), intersection);
//                            // the child of mrcaXA that has cladeX is the child of nodeBcopy
//                            for(Object o: mrcaXA.getChildren()){
//                                NetNode child = (NetNode) o;
//                                if (Utils.getLeafSet(child).containsAll(cladeX)){
//                                    mrcaXA.removeChild(child);
//                                    mrcaXA.adoptChild(nodeBcopy, NetNode.NO_DISTANCE);
//                                    parentBX.adoptChild(child, NetNode.NO_DISTANCE);
//                                    break;
//                                }
//                            }
//                        }
//                        else{
//                            newNodeToAdd.adoptChild(nodeAinT, NetNode.NO_DISTANCE);
//                            newNodeToAdd.setData(new NetNodeInfo());
//                            newNodeToAdd.adoptChild(nodeBcopy, NetNode.NO_DISTANCE);
//                            NetNode parentA = nodeAinT.getParents().iterator().next();
//                            parentA.removeChild(nodeAinT);
//                            parentA.adoptChild(newNodeToAdd, NetNode.NO_DISTANCE);
//                        }
                    }

                } else if (nBinT) {
                    // Add node A to T
                    // (B,X) -> ((A,B),X)), need to check if A and X have overlapped species
                    if (_debug) {
                        System.out.println("Adding node A "+nodeA.getName()+" to T"+i);
                    }
                    edits.set(i, true);
                    Map<NetNode, NetNode> old2new = new HashMap<>();
                    NetNode<NetNodeInfo> nodeAcopy = Utils.copyNodeSubnet(nodeA, old2new);
                    nodeAcopy.setName(nodeA.getName());
                    nodeBinT.setName(nodeB.getName());
                    if (leaf.equals(cladeB)) {
                        newNodeToAdd.adoptChild(nodeBinT, NetNode.NO_DISTANCE);
                        newNodeToAdd.setData(new NetNodeInfo());
                        newNodeToAdd.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
                        _subnetworks.set(i, new BniNetwork((BniNetNode) newNodeToAdd));
                    } else {//there could be several networks
                        Set<String> cladeX = new HashSet<>(leaf);
                        cladeX.removeAll(cladeB);
                        checkOverlap(cladeX, cladeA, cladeB, _subnetworks.get(i), nodeBinT, nodeAcopy, newNodeToAdd);
//                        if (!Collections.disjoint(cladeX, cladeA)) {
//                            // node under nodeB to node in network
//                            NetNode parentAX = null;
//                            Map<NetNode, NetNode> redundantNodeMap = getRedundantNode(nodeBinT, nodeAcopy, _subnetworks.get(i));
//                            for (NetNode nodeunderA : redundantNodeMap.keySet()) {
//                                parentAX =  (NetNode) nodeunderA.getParents().iterator().next();
//                                parentAX.removeChild(nodeunderA);
//
//                            }
//                            Set<String> intersection = new HashSet<>(cladeX);
//                            intersection.retainAll(cladeA);
//                            intersection.addAll(cladeB);
//                            NetNode mrcaXB = getNodeFromClade(_subnetworks.get(i), intersection);
//                            // the child of mrcaXA that has cladeX is the child of nodeBcopy
//
//                            for(Object o: mrcaXB.getChildren()){
//                                NetNode child = (NetNode) o;
//                                if (Utils.getLeafSet(child).containsAll(cladeX)){
//                                    mrcaXB.removeChild(child);
//                                    mrcaXB.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
//                                    mrcaXB.adoptChild(child, NetNode.NO_DISTANCE);
//                                    break;
//                                }
//                            }
//
//                        }
//                        else{
//                            NetNodeInfo ni  = new NetNodeInfo();
//                            newNodeToAdd.setData(ni);
//                            NetNode<NetNodeInfo> parentB = nodeBinT.getParents().iterator().next();
//                            parentB.removeChild(nodeBinT);
//                            parentB.adoptChild(newNodeToAdd, NetNode.NO_DISTANCE);
//                            newNodeToAdd.adoptChild(nodeBinT, NetNode.NO_DISTANCE);
//                            newNodeToAdd.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
//                        }

                    }
                }
            }
        }
        Map<NetNode, NetNode> old2new = new HashMap<>();
        if(newNodeToAdd.getName().equals("I9") && _debug){
            System.out.println("debug start: copyNodeSubnet");
        }
        NetNode newNode = Utils.copyNodeSubnet(newNodeToAdd, old2new);
        return new Tuple<>(edits, newNode);
    }

    public List<Network> joinNodesInOneNet(Network net1, NetNode nodeAinT1, Set<String> cladeA, Network net2, NetNode nodeBinT2, Set<String> cladeB){
        /*
        Join clade A and clade B in just one tree

        1. tree 1 as (A,X); and extract (A,X);
        2. tree 2 as (B,Y); and extract (B,Y);
        3. attach A to new tree 2 as (A,(B,Y));
        */
        Map<NetNode, NetNode> old2new = new HashMap<>();
        NetNode<NetNodeInfo> nodeAcopy = Utils.copyNodeSubnet(nodeAinT1, old2new);
        nodeAcopy.setName(nodeAinT1.getName());

        Set<String> leaf = Utils.getLeafSet(net2.getRoot());
        NetNode<NetNodeInfo> newNodeToAdd = new BniNetNode();
        newNodeToAdd.setName("I"+_nodeCnt);

        if (leaf.equals(cladeB)) {
            newNodeToAdd.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
            newNodeToAdd.setData(new NetNodeInfo());
            newNodeToAdd.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
            net2 = new BniNetwork((BniNetNode) newNodeToAdd);
        } else {
            NetNodeInfo ni  = new NetNodeInfo();
            newNodeToAdd.setData(ni);
            NetNode<NetNodeInfo> parentB = (NetNode) nodeBinT2.getParents().iterator().next();
            parentB.removeChild(nodeBinT2);
            parentB.adoptChild(newNodeToAdd, NetNode.NO_DISTANCE);
            newNodeToAdd.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
            newNodeToAdd.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);

        }
        List<Network> netlist = new ArrayList<>();
        netlist.add(net1);
        netlist.add(net2);
        return netlist;
    }


//    public Tuple<List<Boolean>, NetNode<NetNodeInfo>> joinNodes(Network net1, Network net2, NetNode<NetNodeInfo> nodeA, NetNode<NetNodeInfo> nodeB) {
//        Set<String> cladeA = Utils.getLeafSet(nodeA);
//        Set<String> cladeB = Utils.getLeafSet(nodeB);
//
////        if (!Collections.disjoint(cladeA, cladeB)) {
////            throw new RuntimeException("Nodes are not disjoint on their leaf sets!");
////        }
//
//        List<Boolean> edits = Stream.generate(() -> false).limit(_subnetworks.size()).collect(Collectors.toList());
//
//        NetNode<NetNodeInfo> newNodeToAdd = new BniNetNode();
//        newNodeToAdd.setName("I"+_nodeCnt);
//        if ((cladeA.equals(_leaves.get(0)) && cladeB.equals(_leaves.get(1))) || (cladeA.equals(_leaves.get(1)) && cladeB.equals(_leaves.get(0)))){
//            edits.set(0, true);
//            edits.set(1, true);
//
//            newNodeToAdd.adoptChild(net1.getRoot(), NetNode.NO_DISTANCE);
//            newNodeToAdd.adoptChild(net2.getRoot(), NetNode.NO_DISTANCE);
//            newNodeToAdd.setData(new NetNodeInfo());
//            _subnetworks.set(0, new BniNetwork((BniNetNode) newNodeToAdd));
//            _subnetworks.set(1, null);
//
//        }
//        else{
//            newNodeToAdd = new BniNetNode();
//            newNodeToAdd.setName("I"+_nodeCnt);
//            Set<String> leaf1 = _leaves.get(0);
//            Set<String> leaf2 = _leaves.get(1);
//            NetNode<NetNodeInfo> nodeAinT1, nodeAinT2;
//            NetNode<NetNodeInfo> nodeBinT1, nodeBinT2;
//
//            if (cladeA.equals(leaf1)) {
//                nodeAinT1 = net1.getRoot();
//
//            } else {
//                nodeAinT1 = getNodeFromClade(net1, cladeA);
//            }
//
//            if (cladeB.equals(leaf1)) {
//                nodeBinT1 = net1.getRoot();
//            } else {
//                nodeBinT1 = getNodeFromClade(net1, cladeB);
//            }
//
//            if (cladeA.equals(leaf2)) {
//                nodeAinT2 = net2.getRoot();
//
//            } else {
//                nodeAinT2 = getNodeFromClade(net2, cladeA);
//            }
//
//            if (cladeB.equals(leaf2)) {
//                nodeBinT2 = net2.getRoot();
//            } else {
//                nodeBinT2 = getNodeFromClade(net2, cladeB);
//            }
//
//            List<Network> tuplenet = null;
//            boolean nAinT1 = nodeAinT1 != null;
//            boolean nAinT2 = nodeAinT2 != null;
//            boolean nBinT1 = nodeBinT1 != null;
//            boolean nBinT2 = nodeBinT2 != null;
//
//
//            if (nAinT1 && nAinT2) {
//                // NetNode nodeA in *both* net1 and net2, do nothing!
//                if(nBinT1 && nBinT2){
//                    // Case 1: Node A and node B in *both* net1 and net2
//                    // Only valid if (nodeA, nodeB) are siblings in *both* net1 and net2
//                    //Do nothing (except update the node set)
//                    nodeAinT1.setName(nodeA.getName());
//                    nodeAinT2.setName(nodeA.getName());
//                    nodeBinT1.setName(nodeB.getName());
//                    nodeBinT2.setName(nodeB.getName());
//
//                    Set<String> cladeUnion = new HashSet<>();
//                    cladeUnion.addAll(cladeA);
//                    cladeUnion.addAll(cladeB);
//                    newNodeToAdd = getNodeFromClade(net1, cladeUnion);
//                    newNodeToAdd.setName("I"+_nodeCnt);
//                    newNodeToAdd = getNodeFromClade(net2, cladeUnion);
//                    newNodeToAdd.setName("I"+_nodeCnt);
//                }
//                else if (nBinT1){
//                    // Case 2: Node A in both net1 and net2, node B in net1
//                    // Add node B to net2
//                    edits.set(1, true);
//                    tuplenet = joinNodesInOneNet(net1, nodeBinT1, cladeB, net2, nodeAinT2, cladeA);
//
//                }
//                else if (nBinT2){
//                    // Case 3: Node A in both net1 and net2, node B in net2
//                    // Add node B to net1
//                    edits.set(0, true);
//                    tuplenet = joinNodesInOneNet(net2, nodeBinT2, cladeB, net1, nodeAinT1, cladeA);
//                }
//                else{
//                    throw new RuntimeException("Node B was not found in either nets");
//                }
//
//
//            } else if (nAinT1) {
//                //nodeA in net1 only
//                if (nBinT1 && nBinT2){
//                    // Case 4: Node A in net1, node B in *both* net1 and net2
//                    // Add node A to net2
//                    edits.set(1, true);
//                    tuplenet = joinNodesInOneNet(net1, nodeAinT1, cladeA, net2, nodeBinT2, cladeB);
//                }
//                else if (nBinT1){
//                    // Case 5: Node A in net1, node B in net1
//                    // do nothing
//                    nodeAinT1.setName(nodeA.getName());
//                    nodeBinT1.setName(nodeB.getName());
//
//                    Set<String> cladeUnion = new HashSet<>();
//                    cladeUnion.addAll(cladeA);
//                    cladeUnion.addAll(cladeB);
//                    newNodeToAdd = getNodeFromClade(net1, cladeUnion);
//                    newNodeToAdd.setName("I"+_nodeCnt);
//
//                }
//                else if (nBinT2){
//                    // Case 6: Node A in net1, node B in net2
//                    // Add node A to net2 and node B to net1
//                    edits.set(0, true);
//                    edits.set(1, true);
//                    tuplenet = joinNodesInBothNets(net1, nodeAinT1, cladeA, net2, nodeBinT2, cladeB);
//                }
//                else{
//                    throw new RuntimeException("Node B was not found in either nets");
//                }
//            }
//            else if (nAinT2){
//                //nodeA in net2 only
//                if (nBinT1 && nBinT2){
//                    // Case 7: Node A in net2, node B in *both* net1 and net2
//                    // Add node A to net1
//                    edits.set(0, true);
//                    tuplenet = joinNodesInOneNet(net2, nodeAinT2, cladeA, net1, nodeBinT1, cladeB);
//                }
//                else if (nBinT1){
//                    // Case 8: Node A in net2, node B in net1
//                    // Add node A to net1 and node B to net2
//                    edits.set(0, true);
//                    edits.set(1, true);
//                    tuplenet = joinNodesInBothNets(net1, nodeBinT1, cladeB, net2, nodeAinT2, cladeA);
//                }
//                else if (nBinT2){
//                    // Case 9: Node A in net2, node B in net2
//                    // do nothing
//                    nodeAinT2.setName(nodeA.getName());
//                    nodeBinT2.setName(nodeB.getName());
//
//                    Set<String> cladeUnion = new HashSet<>();
//                    cladeUnion.addAll(cladeA);
//                    cladeUnion.addAll(cladeB);
//                    newNodeToAdd = getNodeFromClade(net2, cladeUnion);
//                    newNodeToAdd.setName("I"+_nodeCnt);
//
//                }
//                else{
//                    throw new RuntimeException("Node B was not found in either nets");
//                }
//            }
//            else{
//                throw new RuntimeException("Node A was not found in either nets");
//
//            }
//            if (tuplenet!=null) {
//                _subnetworks.set(0, net1);
//                _subnetworks.set(1, net2);
//
//            }
//
//
//        }
//        Map<NetNode, NetNode> old2new = new HashMap<>();
//
//
//        NetNode newNode = Utils.copyNodeSubnet(newNodeToAdd, old2new);
//        return new Tuple<>(edits, newNode);
//    }

    //todo
    public static Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> getParentAboveHeightBFS(NetNode<NetNodeInfo> leafnode, double height) {
        Queue<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> queue = new LinkedList<Tuple<NetNode<NetNodeInfo>,NetNode<NetNodeInfo>>>();
        for (NetNode<NetNodeInfo> parent: leafnode.getParents()){
            queue.add(new Tuple<>(parent, leafnode));
        }

        Set<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> visited = new HashSet<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>>();
        while(!queue.isEmpty()) {
            Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge = queue.poll();
            NetNode<NetNodeInfo> parent = edge.Item1;
            NetNode<NetNodeInfo> node = edge.Item2;
            if(visited.contains(edge)) continue;
            visited.add(edge);

            if (node.isNetworkNode() && (node.getParentProbability(parent) != NetNode.NO_PROBABILITY) && node.getParentProbability(parent) < 0.5){
                continue;
            }
//            if (parent.getData().getHeight() > height && node.getData().getHeight() < height){
//                return new Tuple<>(parent, node);
//            }
            for(NetNode<NetNodeInfo> parentParent : parent.getParents()) {
                queue.add(new Tuple<>(parentParent, parent));
            }
        }
        return null;
    }

    //todo get parent
//    private Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> getParentAboveHeight(NetNode<NetNodeInfo> node, double height){
//        for (NetNode<NetNodeInfo> parent: node.getParents()){
//            if (node.isNetworkNode() && (node.getParentProbability(parent) != NetNode.NO_PROBABILITY) &&  node.getParentProbability(parent) < 0.5){
//                continue;
//            }
//            if (parent.getData().getHeight() > height){
//                return new Tuple<>(parent, node);
//            }
//            else {
//                return getParentAboveHeight(parent, height);
//            }
//
//        }
//        return null;
//    }

//    public List<Boolean> joinNodes(Network net1, Network net2, NetNode nodeA, NetNode nodeB) {
//        Set<String> cladeA = getLeafSet(nodeA);
//        Set<String> cladeB = getLeafSet(nodeB);
//
//        if (!Collections.disjoint(cladeA, cladeB)) {
//            throw new RuntimeException("Nodes are not disjoint on their leaf sets!");
//        }
//
//        List<Boolean> edits = Stream.generate(() -> false).limit(_subnetworks.size()).collect(Collectors.toList());
//
//        for (int i = 0; i < _subnetworks.size(); i++) {
//            Set<String> leaf = _leaves.get(i);
//            NetNode nodeAinT;
//            NetNode nodeBinT;
//
//            if (cladeA.equals(leaf)) {
//                nodeAinT = nodeA;
//            } else {
//                nodeAinT = getNodeFromClade(_subnetworks.get(i), cladeA);
//            }
//
//            if (cladeB.equals(leaf)) {
//                nodeBinT = nodeB;
//            } else {
//                nodeBinT = getNodeFromClade(_subnetworks.get(i), cladeB);
//            }
//
//            boolean nAinT = nodeAinT != null;
//            boolean nBinT = nodeBinT != null;
//
//            if (nAinT && nBinT) {
//                // NetNode A and node B are both in T, do nothing!
//                continue;
//            } else if (nAinT) {
//                // Add node B to T
//                edits.set(i, true);
//                NetNode nodeBcopy = copyNodeSubnet(nodeB);
//                if (leaf.equals(cladeA)) {
//                    NetNode root = new BniNetNode();
//                    root.adoptChild(nodeAinT, NetNode.NO_DISTANCE);
//                    root.adoptChild(nodeBcopy, NetNode.NO_DISTANCE);
//                    _subnetworks.set(i, new BniNetwork((BniNetNode) root));
//                }
//                else{
//                    NetNode parentA = (NetNode) nodeAinT.getParents().iterator().next();
//                    NetNode newParentNodeT = new BniNetNode();
//                    parentA.removeChild(nodeAinT);
//                    parentA.adoptChild(newParentNodeT, NetNode.NO_DISTANCE);
//                    newParentNodeT.adoptChild(nodeAinT, NetNode.NO_DISTANCE);
//                    newParentNodeT.adoptChild(nodeBcopy, NetNode.NO_DISTANCE);
//                }
//
//            } else if (nBinT) {
//                // Add node A to T
//                edits.set(i, true);
//                NetNode nodeAcopy = copyNodeSubnet(nodeA);
//                if (leaf.equals(cladeB)) {
//                    NetNode root = new BniNetNode();
//                    root.adoptChild(nodeB, NetNode.NO_DISTANCE);
//                    root.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
//                    _subnetworks.set(i, new BniNetwork((BniNetNode) root));
//                } else {
//                    NetNode parentB = (NetNode) nodeBinT.getParents().iterator().next();
//                    NetNode newParentNodeT = new BniNetNode();
//                    parentB.removeChild(nodeBinT);
//                    parentB.adoptChild(newParentNodeT, NetNode.NO_DISTANCE);
//                    newParentNodeT.adoptChild(nodeBinT, NetNode.NO_DISTANCE);
//                    newParentNodeT.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
//                }
//            }
//        }
//        //Todo: make the topology difference of shared subset due to reticulations the same
//
//        return edits;
//    }

    public static void dfs_parents(NetNode node, Set<NetNode> visited, Set<String> clade, List<NetNode> result){
        Set<String> leaves = Utils.getLeafSet(node);
        if (leaves.containsAll(clade)){
            if (result.isEmpty()){
                result.add(node);
            }
            return;
        }
        visited.add(node);

        for (Object parent: node.getParents()){
            NetNode parentnode = (NetNode) parent;
            if (!visited.contains(parentnode)) {
                dfs_parents(parentnode, visited, clade, result);
            }

        }
    }

    public static NetNode bfs_parents(NetNode node, Set<String> clade){
//        Set<String> leaves = Utils.getLeafSet(node);
//        if (leaves.containsAll(clade)){
//            if (result.isEmpty()){
//                result.add(node);
//            }
//            return;
//        }
//        visited.add(node);
        Queue<NetNode> queue = new LinkedList<>();
        Set<NetNode> visited = new HashSet<>();
        queue.add(node);
        visited.add(node);

        while (!queue.isEmpty()){
            NetNode curnode = queue.poll();
            Set<String> leaves = Utils.getLeafSet(curnode);
            if (leaves.containsAll(clade)){
                return curnode;
            }
            for (Object parent: curnode.getParents()){
                NetNode parentnode = (NetNode) parent;
                if (!visited.contains(parentnode)) {
                    queue.add(parentnode);
                    visited.add(parentnode);
                }

            }
        }
        return null;

    }

    //TODO: debug here
    public static NetNode getMRCA(Network net, Set<String> clade){
        NetNode mrca = null;
        Iterator<String> it = clade.iterator();
        String leaf = it.next();
        mrca = bfs_parents(net.findNode(leaf), clade);

        return mrca;
     }


    public static NetNode getNodeFromClade(Network net, Set<String> clade) {
        // Returns the MRCA node of the clade igoring the branch lengths

        Set<String> leaves = Utils.getLeafSet(net.getRoot());

        // Check if the clade is the whole tree!
        if (leaves.equals(new HashSet<>(clade))) {
            return net.getRoot();
        }

        // Check if the clade contains leaves not in the tree itself
        if (!leaves.containsAll(clade)) {
            return null;
        }

        NetNode node = getMRCA(net, clade);

        Set<String> nodeLeaves = getLeafSetButReti(node);
        if(clade.containsAll(nodeLeaves)){
            return node;
        }
        else{
            return null;
        }
    }


    // Todo: how compatible: identical topology?
    public boolean areNetsCompatible(Network net1, Network net2){
        Set<String> net1Leaves = new HashSet<>();
        Set<String> net2Leaves = new HashSet<>();
        for (Object o: net1.getLeaves()){
            NetNode leaf = (NetNode) o;
            net1Leaves.add(leaf.getName());
        }

        for (Object o: net2.getLeaves()){
            NetNode leaf = (NetNode) o;
            net2Leaves.add(leaf.getName());
        }

        Set<String> intersection = new HashSet<>();
        intersection.addAll(net1Leaves);;
        intersection.retainAll(net2Leaves);

        System.out.println("areNetsCompatible:");
        System.out.println(net1.toString());
        System.out.println(net2.toString());
        if (_debug){
            System.out.println(intersection);
        }

        if (intersection.size() <= 1){
            return true;
        }

        if (intersection.contains(this._outgroup)){
            net1.resetRoot(this._outgroup);
            net2.resetRoot(this._outgroup);
        }
        Network subnet1 = SuperNetwork3.getSubNetwork(net1, new ArrayList<>(intersection), true).Item1;
        Network subnet2 = SuperNetwork3.getSubNetwork(net2, new ArrayList<>(intersection), true).Item1;


        Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(subnet1, subnet2);
//        if (Networks.hasTheSameTopology(subnet1.Item1, subnet2.Item1)){
//            return true;
//        }
        if (closest.Item3 < 0.01){
            System.out.println("Compatible");
            return true;
        }
        else{
            if (_debug){
                System.out.println("Incompatible");
                System.out.println(net1.toString());
                System.out.println(net2.toString());
            }

            return false;
        }
    }

    /*
        Test whether joining cladeA and cladeB in one
        or both networks causes the two networks to be incompatible
    */
    public boolean testJoin(NetNode nodeA, NetNode nodeB) throws Exception {
        System.out.println("testjoin:"+nodeA.getName()+","+nodeB.getName());
        if (nodeA.getName().equals("I10")  && nodeB.getName().equals("I14")){
            System.out.println("checking point");
        }

//        Set<String> cladeA = getLeafSetButReti(nodeA);
//        Set<String> cladeB = getLeafSetButReti(nodeB);
        Set<String> cladeA = Utils.getLeafSet(nodeA);
        Set<String> cladeB = Utils.getLeafSet(nodeB);
        Set<String> cladeAB = new HashSet<>(cladeA);
        cladeAB.addAll(cladeB);

//        cladeA = getLeafSetButReti(nodeAinTs.get(0));
//        cladeB = getLeafSetButReti(nodeBinTs.get(0));
//        if (!Collections.disjoint(cladeA, cladeB)) {
//            throw new Exception("Nodes are not disjoint on their leaf sets! "+nodeA.getName()+":"+cladeA+";"+nodeB.getName()+":"+cladeB);
//        }

        List<NetNode> nodeAinTs = new ArrayList<>();
        List<NetNode> nodeBinTs = new ArrayList<>();

        for (int i = 0; i < _subnetworks.size(); i++) {
            Set<String> leaf = this._leaves.get(i);

            if (cladeA.equals(leaf)) {
                nodeAinTs.add(nodeA);
            } else {
                NetNode temp = getNodeFromClade(_subnetworks.get(i), cladeA);
                nodeAinTs.add(temp);
            }

            if (cladeB.equals(leaf)) {
                nodeBinTs.add(nodeB);
            } else {
                NetNode temp = getNodeFromClade(_subnetworks.get(i), cladeB);
                nodeBinTs.add(temp);
            }
        }

        List<Boolean> nAinTs = new ArrayList<>();
        for (NetNode nodeAinT : nodeAinTs) {
            nAinTs.add(nodeAinT != null);
        }
        if (Collections.frequency(nAinTs, true) < 1) {
            throw new Exception("Node A: "+nodeA.getName()+" was not found in any tree!");
        }

        List<Boolean> nBinTs = new ArrayList<>();
        for (NetNode nodeBinT : nodeBinTs) {
            nBinTs.add(nodeBinT != null);
        }


        if (Collections.frequency(nBinTs, true) < 1) {
            throw new Exception("Node B: "+nodeB.getName()+" was not found in any tree!--"+cladeB);
        }



        boolean violates = false;
        int i =0;
        boolean nAinT1 = nAinTs.get(i);
        boolean nBinT1 = nBinTs.get(i);
        NetNode nodeAinT1 = nodeAinTs.get(i);
        NetNode nodeBinT1 = nodeBinTs.get(i);

        int j = 1;
        boolean nAinT2 = nAinTs.get(j);
        boolean nBinT2 = nBinTs.get(j);
        NetNode nodeAinT2 = nodeAinTs.get(j);
        NetNode nodeBinT2 = nodeBinTs.get(j);
        Network n1 = _subnetworks.get(i).clone();
        Network n2 = _subnetworks.get(j).clone();

        // TODO: add cases
        if ((nAinT1 || nAinT2) && (nBinT1 || nBinT2)) {
            if (nAinT1 && nAinT2) {
                //nodeA in *both* T1 and T2
                if (nBinT1 && nBinT2) {
                    // Case 1: nodeB in *both* T1 and T2
                    // Valid if nodeA and nodeB are siblings in both T1 & T2
                    NetNode node1 = getNodeFromClade(_subnetworks.get(i), cladeAB);
                    NetNode node2 = getNodeFromClade(_subnetworks.get(j), cladeAB);
                    if (node1 == null || node2 == null) {
                        violates = true;
                    }
                } else if (nBinT1) {
                    //Case 2: node B in T1 only
                    // Valid if nodeA and nodeB are siblings in T1
                    NetNode node = getNodeFromClade(_subnetworks.get(i), cladeAB);
                    if (node == null) {
                        violates = true;
                    }
                } else if (nBinT2) {
                    // Case 3: Node B in T2 only
                    // Valid if nodeA and nodeB are siblings in T2
                    NetNode node = getNodeFromClade(_subnetworks.get(j), cladeAB);
                    if (node == null) {
                        violates = true;
                    }
                } else {
                    throw new Exception("Node B not found in either tree!");
                }
            } else if (nAinT1) {
                // nodeA in T1 only
                if(nBinT1 && nBinT2){
                    // Case 4: nodeB in *both* T1 and T2
                    NetNode node = getNodeFromClade(_subnetworks.get(i), cladeAB);
                    if (node == null) {
                        violates = true;
                    }
                }
                else if (nBinT1){
                    // Case 5: Node B in T1 only
                    // Valid if nodeA and nodeB are siblings in T1
                    NetNode node = getNodeFromClade(_subnetworks.get(i), cladeAB);
                    if (node == null){
                        violates = true;
                    }
                }
                else if (nBinT2){
                    // Case 6: Node B in T2 only
                    // Do join in both trees and test for compatibility
//                            Network n1 = _subnetworks.get(i).clone();
//                            Network n2 = _subnetworks.get(j).clone();
                    NetNode nA = getNodeFromClade(n1, cladeA);
                    NetNode nB = getNodeFromClade(n2, cladeB);

                    List<Network> networklist = joinNodesInBothNets(n1, nA, cladeA, n2, nB, cladeB);
                    n1 = networklist.get(0);
                    n2 = networklist.get(1);
                    if (n1 == null && n2 == null){
                        violates = true;
                    }
                    else if (n1 !=  null){
                        violates = (!areNetsCompatible(n1, n2));//TODO
                    }
                }
                else{
                    throw new Exception("Node B not found in either network!\n");
                }

            } else if (nAinT2) {
                //nodeA in T2 only
                if (nBinT1 && nBinT2){
                    //Case 7: nodeB in and T2
                    NetNode node = getNodeFromClade(_subnetworks.get(j), cladeAB);
                    if (node ==  null){
                        violates = true;
                    }

                }
                else if (nBinT1){
                    //Case 8 (reverse of Case 6): Node B in T1 only
//                    Network n1 = _subnetworks.get(i).clone();
//                    Network n2 = _subnetworks.get(j).clone();

                    NetNode nB = getNodeFromClade(n1, cladeB);
                    NetNode nA = getNodeFromClade(n2, cladeA);
                    List<Network> networklist = joinNodesInBothNets(n1, nB, cladeB, n2, nA, cladeA);
                    n1 = networklist.get(0);
                    n2 = networklist.get(1);
                    if (n1 == null && n2 == null){
                        violates = true;
                    }
                    else if (n1 !=  null){
                        violates = (!areNetsCompatible(n1, n2));//TODO
                    }
                }
                else if (nBinT2){
                    // Case 9: Node B in T2 only
                    // Only valid if (nodeA, nodeB) are siblings in T2
                    NetNode node = getNodeFromClade(_subnetworks.get(j), cladeAB);
                    if (node ==  null){
                        violates = true;
                    }
                }

                else{
                    throw new Exception("Node B not found in either network!\n");
                }

            } else {
                throw new Exception("Node A not found in either tree!");
            }

            if (violates) {
                return violates;
            }
        }
        if(violates) {
            System.out.println("violates");
        }
        return violates;
    }

//    public static void calculateDistances(Network network1, Network network2,  List<NetNode<NetNodeInfo>> nodePool){
//        Set<NetNode> nodeSet = new HashSet<>(nodePool);
//
//        for (NetNode node: node2integers0.keySet()){
//            if (nodeSet.contains(node)){
//
//            }
//        }
//    }


    public List<Network> mergeNetsViaNJ() {

        // Check trees are on disjoint leaf sets
        for (int i = 0; i < _leaves.size() - 1; i++) {
            for (int j = i + 1; j < _leaves.size(); j++) {
                Set<String> intersection = new HashSet<String>(_leaves.get(i));
                intersection.retainAll(_leaves.get(j));
                if (!intersection.isEmpty()) {
                    System.out.println(_subnetworks.get(i).toString());
                    System.out.println(_subnetworks.get(j).toString());
                    throw new RuntimeException("Input trees are not on disjoint leaf sets!");
                }
            }
        }

        //    TODO: Check distance matrix and trees have matching leaf sets
        //    Remove some extra nonsense
        //    Map splits to nodes
        //    Initialize node pool
        List<NetNode<NetNodeInfo>> nodePool = new ArrayList<>();

        for (String ss: this._taxonList) {
            NetNodeInfo ni = new NetNodeInfo(0.0);
            NetNode<NetNodeInfo> nd = new BniNetNode();
            nd.setData(ni);
            nd.setName(ss);
            nodePool.add(nd);
        }

        for (int i = 0; i < this._taxonList.size(); i++) {
            NetNode<NetNodeInfo> nd1 = nodePool.get(i);
            for (int j = 0; j < this._taxonList.size(); j++) {
                if (i == j){
                    continue;
                }
                double d = _matrix[i][j];
                nd1.getData().setNJDistances(_taxonList.get(j), d);
                nd1.getData().incremetNJXSub(d);
            }
        }

        int n = _num_taxa;
        int mergedNetworkID = -1;

        while(n > 1){
            if (_debug){
                System.out.printf("%d joins to go!\n", n);
            }


            // Sort the Q-matrix
            List<Tuple<Integer, Integer>> pairs = new ArrayList<>();
            List<Double> qvalues = new ArrayList<>();
            List<Double> subnetDist = new ArrayList<>();
//            for(int i = 0; i < nodePool.size(); i++){
//                subnetDist.add(0.0);
//            }


//            Map<String, Integer> node2levels0 = Utils.calculateLevels(_subnetworks.get(0));
//            Map<String, Integer> node2levels1 = Utils.calculateLevels(_subnetworks.get(1));
            Set<String> nodeset0 = Utils.getNodeSet(_subnetworks.get(0));
            Set<String> nodeset1 = Utils.getNodeSet(_subnetworks.get(1));
            for (int idx1 = 0; idx1 < nodePool.size() - 1; idx1++) {
                NetNode<NetNodeInfo> nd1 = nodePool.get(idx1);
                for (int idx2 = idx1 + 1; idx2 < nodePool.size(); idx2++ ) {
                    NetNode<NetNodeInfo> nd2 = nodePool.get(idx2);
                    double v1 = (n - 2) * nd1.getData().getNjDistances().get(nd2.getName());
                    double qvalue = v1 - nd1.getData().getNJXsub() - nd2.getData().getNJXsub();
                    pairs.add(new Tuple<>(idx1, idx2));
                    qvalues.add(qvalue);
                    if (nodeset0.contains(nd1.getName()) && nodeset0.contains(nd2.getName())){
                        subnetDist.add(Utils.getDistanceBetweenTwoNodes(_subnetworks.get(0), nd1.getName(), nd2.getName()));
                    }
                    else if (nodeset1.contains(nd1.getName()) && nodeset1.contains(nd2.getName())){
                        subnetDist.add(Utils.getDistanceBetweenTwoNodes(_subnetworks.get(1), nd1.getName(), nd2.getName()));
                    }
                    else{
                        subnetDist.add(Double.MAX_VALUE);
                    }

                }
            }

            if (_debug){
                System.out.println("before:");
                System.out.println(_subnetworks.get(0).toString());
                System.out.println(_subnetworks.get(1).toString());
            }


            // Test for constraint violations
            // TODO: Use multi-threading in test_join function!
            Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> nodesToJoin = null;

            // Todo: we sort the pairwise distances, when multiple pairs have the same distance, we choose the one that are closest in the subnetworks
//            int[] sortedIndices = IntStream.range(0, qvalues.size())
//                    .boxed().sorted((i, j) -> Double.compare(qvalues.get(i), qvalues.get(j)))
//                    .mapToInt(ele -> ele).toArray();
            int[] sortedIndices = IntStream.range(0, qvalues.size())
                    .boxed()
                    .sorted((i, j) -> {
                        int qCompare = Double.compare(qvalues.get(i), qvalues.get(j));
                        if (qCompare == 0) {
                            // if qvalues are equal, compare subnetDist
                            return Double.compare(subnetDist.get(i), subnetDist.get(j));
                        }
                        return qCompare;
                    })
                    .mapToInt(ele -> ele)
                    .toArray();

//            if(_debug){
//                for (int idxq : sortedIndices) {
//                    Tuple<Integer, Integer> pair = pairs.get(idxq);
//                    NetNode nd1 = nodePool.get(pair.Item1);
//                    NetNode nd2 = nodePool.get(pair.Item2);
//                    System.out.println(nd1.getName()+","+nd2.getName()+";"+qvalues.get(idxq)+";"+subnetDist.get(idxq));
//                }
//            }

            for (int idxq : sortedIndices) {
                Tuple<Integer, Integer> pair = pairs.get(idxq);
                NetNode nd1 = nodePool.get(pair.Item1);
                NetNode nd2 = nodePool.get(pair.Item2);

                // Check join does not violate a constraint subnetwork!
                try{
//                    System.out.println(nd1.getName()+";"+nd2.getName());
                    boolean violates = testJoin(nd1, nd2);
                    if (!violates || n <= 3) {
                        nodesToJoin = new Tuple<>(nd1, nd2);
                        break;
                    }
                }catch (Exception e){
                    e.printStackTrace();
                    System.err.println("Error testJoin(nd1, nd2)!");
                }

            }

            if (nodesToJoin == null) {
                throw new RuntimeException("Unable to find valid siblinghood!");
            }
            if (_debug){
                System.out.println(nodesToJoin.Item1.getName()+";"+nodesToJoin.Item2.getName());
            }


            //Compute distance from the pair members to the new node
            double d_A_new = nodesToJoin.Item1.getData().getNjDistances().get(nodesToJoin.Item2.getName());
            double d_sum = 0;
            for (NetNode<NetNodeInfo> node : nodePool){
                if (node.equals(nodesToJoin.Item1) || node.equals(nodesToJoin.Item2)){
                    continue;
                }
                d_sum += node.getData().getNjDistances().get(nodesToJoin.Item1.getName());
                d_sum -= node.getData().getNjDistances().get(nodesToJoin.Item2.getName());
            }

            double sig_A = 0;
            double sig_B = 0;
            if (n > 2){
                d_sum /= ((n-2)*2);
                sig_A = d_A_new/2 + d_sum;
                sig_B = d_A_new - sig_A;
                if(_debug){
                    System.out.println(sig_A);
                }



            }


            Tuple<List<Boolean>, NetNode<NetNodeInfo>> tuple = joinNodes(nodesToJoin.Item1, nodesToJoin.Item2); // TODO: Add required parameters
//            Tuple<List<Boolean>, NetNode<NetNodeInfo>> tuple = joinNodes(_subnetworks.get(0), _subnetworks.get(1), nodesToJoin.Item1, nodesToJoin.Item2);
//
            if (_debug){
                System.out.println("after:");
                System.out.println(_subnetworks.get(0).toString());
                System.out.println(_subnetworks.get(1).toString());
                System.out.println("--------------");
            }

            List<Boolean> edits = tuple.Item1;
            NetNode<NetNodeInfo> newNode = tuple.Item2;
            if (edits.stream().anyMatch(Boolean::booleanValue)) {
                int i = 0;

                for (int index = 0; index < _subnetworks.size(); index++) {
                    Network net = _subnetworks.get(index);
                    boolean e = edits.get(index);

                    if (e) {
                        // Check to see if you can quit early
                        Set<String> leafSet = new HashSet<>();
                        for(Object o: net.getLeaves()){
                            NetNode node = (NetNode) o;
                            leafSet.add(node.getName());

                        }
                        _leaves.set(i, leafSet);

                        if (leafSet.equals(new HashSet<>(_taxonList))) {
                            //todo something
                            // System.out.println("Able to quit early!");
                            mergedNetworkID = index;
                            if (_debug)
                            {
                                System.out.println("equal leaf set:"+net.toString());
                            }
                            break;

                        }

                    }
                    i++;
                }
            }

            if (mergedNetworkID >= 0){
                break;
            }
            // Create the new node
//            NetNode<NetNodeInfo> newNode = new BniNetNode();

            nodePool.remove(nodesToJoin.Item1);
            nodePool.remove(nodesToJoin.Item2);

            // Calculate the distances for the new node
            double d_fu = nodesToJoin.Item1.getData().getNjDistances().get(nodesToJoin.Item2.getName());
            double d = 0;
            for (NetNode<NetNodeInfo> node : nodePool){
                d += node.getData().getNjDistances().get(nodesToJoin.Item1.getName());
                d -= node.getData().getNjDistances().get(nodesToJoin.Item2.getName());
            }
            d /= ((n-2)*2);
            double sig_fu = d_fu/2 + d;

            // Attach it to the tree
//            newNode.setName("["+nodesToJoin.Item1.getName()+"],["+nodesToJoin.Item2.getName()+"]");
//            newNode.setName("I"+_nodeCnt);
            _nodeCnt += 1;
//            newNode.adoptChild(nodesToJoin.Item1, sig_fu);
//            newNode.adoptChild(nodesToJoin.Item2, d_fu - sig_fu);


//            NetNodeInfo ni = new NetNodeInfo();
//            ni.zeroNJXSub();
            newNode.getData().zeroNJXSub();
            newNode.getData().emptyNJDistances();
            for (NetNode<NetNodeInfo> node : nodePool) {
                // actual node-to-node distances
                double v1 = node.getData().getNjDistances().get(nodesToJoin.Item1.getName());
                v1 += node.getData().getNjDistances().get(nodesToJoin.Item2.getName());

                double v3 = nodesToJoin.Item1.getData().getNjDistances().get(nodesToJoin.Item2.getName());
                double dist = 0.5 * (v1 - v3);
                newNode.getData().getNjDistances().put(node.getName(), dist);
                node.getData().getNjDistances().put(newNode.getName(), dist);//TODO: add name for each node

                // Adjust/recalculate the values needed for the Q-matrix calculations
                newNode.getData().incremetNJXSub(dist);
                node.getData().incremetNJXSub(dist);
                node.getData().incremetNJXSub(-nodesToJoin.Item1.getData().getNjDistances().get(node.getName()));
                node.getData().incremetNJXSub(-nodesToJoin.Item2.getData().getNjDistances().get(node.getName()));
            }
            // Clean up
            nodesToJoin.Item1.getData().emptyNJDistances();
            nodesToJoin.Item1.getData().setNJXsub(0.0);
            nodesToJoin.Item2.getData().emptyNJDistances();
            nodesToJoin.Item2.getData().setNJXsub(0.0);

            // Add the new node to the pool of nodes
            nodePool.add(newNode);
            if (_debug) {
                for (NetNode<NetNodeInfo> node : nodePool) {
                    System.out.print(node.getName()+",");
                }
                System.out.println("\n");
            }
            // Adjust count


            n--;
        }
        //todo
        List<Network> candidateNetworks = new ArrayList<>();
        if (Networks.hasTheSameTopology(_subnetworks.get(0), _subnetworks.get(1)) ){
            //|| _subnetworks.get(0).getReticulationCount() == 0
            candidateNetworks.add(_subnetworks.get(0));
        }
        else{
            for (int i = 0; i < _subnetworks.size(); i++) {
                if ( i == mergedNetworkID) {
                    continue;
                }
                Network net2 = _subnetworks.get(i);
                if (candidateNetworks.size() == 0){
                    Network net1 = _subnetworks.get(mergedNetworkID);
                    if (_debug) {
                        System.out.println("match networks:");
                        System.out.println(net2.toString());
                        System.out.println(net1.toString());
                    }
//                    System.out.println(net2.toString());
//                    System.out.println(net1.toString());
//                    net2.resetRoot(this._outgroup);
//                    net2.resetRoot(this._outgroup);
                    candidateNetworks = matchNetworks(net2, net1);
                    candidateNetworks.addAll(matchNetworks(_originalSubNetworks.get(i).clone(), net1));
                    if (_debug) {
                        System.out.println(candidateNetworks);
                    }

                }

            }
        }
        if (_debug){
            System.out.println("node pool size:"+nodePool.size());
            for (Network net: candidateNetworks){
                System.out.println(net.toString());
            }
        }

        List<Double> scorelist = new ArrayList<>();
        for(Network cn: candidateNetworks){
            scorelist.add(scoreNetworks(cn));
        }
        List<Network> resultnetwork = new ArrayList<>();
        int[] scoreIndexes = Utils.argsort(scorelist, true);


        for (int i = 0; i < scoreIndexes.length; i++){
            boolean alreadyAdd = false;
            for (Network resnet: resultnetwork){
                if (Networks.hasTheSameTopology(resnet, candidateNetworks.get(scoreIndexes[i]))){
                    alreadyAdd = true;
                    break;
                }
            }
            if (!alreadyAdd){
                resultnetwork.add(candidateNetworks.get(scoreIndexes[i]));
            }

            if (i + 1 < scorelist.size()-1 && scorelist.get(scoreIndexes[i]) < scorelist.get(scoreIndexes[i+1])){
                break;
            }
        }
        System.out.println("start checking");
        for (Network resnet: resultnetwork){
            System.out.println(resnet.toString());
        }

        return resultnetwork;
    }

    public double scoreNetworks(Network network){
        double score = 0;
        for (Network originalsubnet: _originalSubNetworks){
            List<String> clade = new ArrayList<>(Utils.getLeafSet(originalsubnet.getRoot()));
            Tuple<Network, Map<NetNode, NetNode>> tuple = SuperNetwork3.getSubNetwork(network, clade, true);
            score += Networks.computeDistanceBetweenTwoNetworks(originalsubnet, tuple.Item1);

            }
        return score;
    }






    public static void testCase1(){
        Network net1 = Networks.readNetwork("((((B:1.0)#H1:1,C:2.0):1,#H1:2):1,D:4);");
        Network net2 = Networks.readNetwork("((A:4,(E:2)#H1:2):1,(#H1:1,F:3):2);");
        double [][] matrix = {{0,2,3,4,4,5},
                {2,0,3,4,4,5},
                {3,3,0,4,4,5},
                {4,4,4,0,1,5},
                {4,4,4,1,0,5},
                {5,5,5,5,5,0}};
        for (int i = 0; i < matrix.length; i++){
            for (int j =0; j < matrix.length; j++){
                matrix[i][j] *= 2;
            }
        }
        List<Network> netlist = new ArrayList<>();
        netlist.add(net1);
        netlist.add(net2);

        List<String> taxonlist = new ArrayList<>();
        taxonlist.add("A");
        taxonlist.add("B");
        taxonlist.add("C");
        taxonlist.add("D");
        taxonlist.add("E");
        taxonlist.add("F");

        NJMergeTopology merge = new NJMergeTopology(netlist, matrix, taxonlist);
        List<Network> net = merge.mergeNetsViaNJ();
        System.out.println(net.toString());
    }


    public static void testCase2(){
        Network net1 = Networks.readNetwork("((A:2,(B:1)#H1:1):1,(C:2,#H1:1):1);");
        Network net2 = Networks.readNetwork("((((D:1,E:1):1)#H1:1,F:3):2,#H1:3);");
        double [][] matrix = {{0,2,3,4,4,5},
                {2,0,3,4,4,5},
                {3,3,0,4,4,5},
                {4,4,4,0,1,5},
                {4,4,4,1,0,5},
                {5,5,5,5,5,0}};
        for (int i = 0; i < matrix.length; i++){
            for (int j =0; j < matrix.length; j++){
                matrix[i][j] *= 2;
            }
        }

        List<Network> netlist = new ArrayList<>();
        netlist.add(net1);
        netlist.add(net2);

        List<String> taxonlist = new ArrayList<>();
        taxonlist.add("A");
        taxonlist.add("B");
        taxonlist.add("C");
        taxonlist.add("D");
        taxonlist.add("E");
        taxonlist.add("F");

        NJMergeTopology merge = new NJMergeTopology(netlist, matrix, taxonlist);
        List<Network> net = merge.mergeNetsViaNJ();
        System.out.println(net.toString());
    }


    public static void testNodeMap(){
        Network net1 = Networks.readNetwork("((F:3.0,((D,E))#H1:1.0):2.0,((C,(B,A)),#H1:2.0):1.0);");
        Network net2 = Networks.readNetwork("((F,(E,D)),(((A,B))#H1:2.0,(C:2.0,#H1:1.0):1.0):1.0);");
        if (net1.getReticulationCount() <= net2.getReticulationCount()){
            Map<NetNode, NetNode> nodemap = NJMerge.mapSubNets(net1, net2);
            System.out.println(net1.toString());
            System.out.println(net2.toString());
            for(NetNode n1: nodemap.keySet()){
                System.out.println(n1+"--"+nodemap.get(n1));

            }
        }

    }

    public static void testGetNodeLeafMap(){
        Map<NetNode, Set<String>> node2leaf = new HashMap<>();
        Network net1 = Networks.readNetwork("((t31,((t20:1.305398272)#H1:0.0::0.5935451058,(t34:0.7979988744,(t26:0.2254407615,t32:0.2254407615)I2:0.5725581129)I3:0.5073993979)I4)I5,(((t13,t22)I0,t2)I1,#H1:0.0::0.4064548942):5.828052169100001);");
        Networks.autoLabelNodes(net1);
        System.out.println(net1.toString());
        getNodeLeafMap(net1.getRoot(), node2leaf, null);
        for(NetNode node: node2leaf.keySet()){
            System.out.println(node.getName()+":"+node2leaf.get(node));
        }
    }

    public static void testGetParent(){
        Network net1 = Networks.readNetwork("((E:1.0,D:1.0):3.0,((B:1.0)#H1:2.0,(C:2.0,#H1:1.0):1.0):1.0);");
        Networks.autoLabelNodes(net1);
        System.out.println(net1);
        NetNode leaf = net1.findNode("B");

        double height = 2.5;
        Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge = getParentAboveHeightBFS(leaf, height);
        System.out.println(edge.Item1.getName());
        System.out.println(edge.Item2.getName());

    }

    public static void testJoinTwoSubnets(){
        Network net1 = Networks.readNetwork("((A:2,(B:1)#H1:1):1,(C:2,#H1:1):1);");
        Network net2 = Networks.readNetwork("((((D:1,E:1):1)#H1:1,F:3):2,#H1:3);");
        NetNode root = new BniNetNode();
//                    root.adoptChild(nodeAinT, NetNode.NO_DISTANCE);
        root.adoptChild(net1.getRoot(), 1);
        root.adoptChild(net2.getRoot(), 1);
//                    root.setData(new NetNodeInfo(_subnetworks.get(i).getRoot()));
        System.out.println();
    }

    public static void testMatchNetworks(){
        Network net2 = Networks.readNetwork("((F:3.0,((D:1.0,E:1.0):1.0)#H1:1.0):2.0,((C:3.0,(B:2.0,A:2.0):1.0):1.0,#H1:2.0):1.0);");
//        Network net1 = Networks.readNetwork("((E:1.0,D:1.0):3.0,((A:2.0,(B:1.0)#H1:1.0):1.0,(C:2.0,#H1:1.0):1.0):1.0);");
        Network net1 = Networks.readNetwork("(((E:1.0,D:1.0):3.0,((A:2.0,(B:1.0)#H1:1.0):1.0,(C:2.0,#H1:1.0):1.0):1.0):1,F:5);");
//        initNodeHeightMap(net1);
//        initNodeHeightMap(net2);
        Networks.autoLabelNodes(net1);
        Networks.autoLabelNodes(net2);
//        matchNetworks(net1, net2);
        matchNetworks(net2, net1);//todo debug
    }

    public static void testNetworkRooting(){
        Network network = Networks.readNetwork("(((E:1.0,D:1.0):3.0,F:6.0):0.1,(((B:1.0)#H1:1.0,C:2.0):1.0,(#H1:1.0,A:2.0):1.0):0.1);");
        network.resetRoot(network.findNode("F"));
        System.out.println(network.toString());
    }

    public static void testRooting(){
//        try{
//            STITree stiTree = new STITree("(((A,B),C),((D,E),F));");
//            stiTree.rerootTreeAtEdge("F");
//            System.out.println(stiTree.toString());
//        }catch (Exception e){
//            e.printStackTrace();
//        }

        Network net = Networks.readNetwork("(((A:2,(B:1)#H1:1):1,(C:2,#H1:1):1),((((D:1,E:1):1)#H2:1,F:3):2,#H2:3));");
        Networks.autoLabelNodes(net);
        System.out.println(net.toString());
        net.resetRoot("F");
        System.out.println(net.toString());

    }

//    public static void testCopyNode(){
//        Network net = Networks.readNetwork("(((A:2,(B:1)#H1:1):1,(C:2,#H1:1):1):3,((((D:1,E:1):1)#H2:1,F:3):2,#H2:3):1);");
//        Networks.autoLabelNodes(net);
//        System.out.println(net.toString());
//        Map<NetNode, NetNode> old2new = new HashMap<>();
//        NetNode copy = copyNodeSubnet(net.getRoot(), old2new);
//        System.out.println((new BniNetwork<NetNodeInfo>((BniNetNode<NetNodeInfo>) copy)).toString());
//        Set<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> candidateEdgeList = findEdgeInRange(net, 1, "F", false);
//        for (Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge:candidateEdgeList){
//            System.out.println(edge.Item1.getName()+","+edge.Item2.getName());
//        }
//        System.out.println("------");
//        candidateEdgeList = findEdgeInRange(net, 1, "F", true);
//        for (Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge:candidateEdgeList){
//            System.out.println(edge.Item1.getName()+","+edge.Item2.getName());
//        }
//    }


    public static void testnetwork(){
        Network net = Networks.readNetwork("(((((t1:1.8,#H1:1.4::0.5):2.9464,((t2:0.3052,t3:0.3052):3.658,((t4:0.4)#H1:1.0948::0.5,(t5:0.5832,t6:0.5832):0.9116):2.4684):0.7832):4.7008,(((t7:1.4428,t8:1.4428):4.4308,((t9:0.8492,t10:0.8492):4.3948,(((((t11:1.3768,t12:1.3768):0.4)#H3:0.638::0.6,(t13:0.4624,t14:0.4624):1.9524):0.5004,(t15:2.7952,(t16:2.0,#H3:0.2232::0.4):0.7952):0.12):0.2084,t17:3.1236):2.1204):0.6296):1.3796,((((t18:1.2424,t19:1.2424):0.898,t20:2.1404):2.1012,((t21:0.6644,t22:0.6644):2.2932,(t23:0.5328,t24:0.5328):2.4248):1.284):0.7936,(t25:2.3128,t26:2.3128):2.7224):2.218):2.194):1.1436,(t27:5.3544,((t28:0.6,#H2:0.2::0.4):0.8936,(t29:1.0344,(t30:0.4)#H2:0.6344::0.6):0.4592):3.8608):5.2364):89.4092,Z:100.0);");
        Set<String> clade = new HashSet<>();
        clade.add("t1");
        clade.add("t2");
        clade.add("t3");
        clade.add("t4");
        clade.add("t5");
        clade.add("t6");
        NetNode node = getNodeFromClade(net, clade);
        System.out.println(node);
    }

    public static void testCase3(){
        Network net_original = Networks.readNetwork("(((((t1:1.8,#H1:1.4::0.5):2.9464,((t2:0.3052,t3:0.3052):3.658,((t4:0.4)#H1:1.0948::0.5,(t5:0.5832,t6:0.5832):0.9116):2.4684):0.7832):4.7008,(((t7:1.4428,t8:1.4428):4.4308,((t9:0.8492,t10:0.8492):4.3948,(((((t11:1.3768,t12:1.3768):0.4)#H3:0.638::0.6,(t13:0.4624,t14:0.4624):1.9524):0.5004,(t15:2.7952,(t16:2.0,#H3:0.2232::0.4):0.7952):0.12):0.2084,t17:3.1236):2.1204):0.6296):1.3796,((((t18:1.2424,t19:1.2424):0.898,t20:2.1404):2.1012,((t21:0.6644,t22:0.6644):2.2932,(t23:0.5328,t24:0.5328):2.4248):1.284):0.7936,(t25:2.3128,t26:2.3128):2.7224):2.218):2.194):1.1436,(t27:5.3544,((t28:0.6,#H2:0.2::0.4):0.8936,(t29:1.0344,(t30:0.4)#H2:0.6344::0.6):0.4592):3.8608):5.2364):89.4092,Z:100.0);");
        List<List<String>> cladeslist = new ArrayList<>();
        List<String> taxonlist = new ArrayList();

        for(Object o: net_original.getLeaves()){
            NetNode n = (NetNode) o;
            taxonlist.add(n.getName());
        }
        List<Network> netlist = new ArrayList<>();

//        "[['t1', 't2', 't3', 't4', 't5', 't6'], ['Z', 't27', 't28', 't29', 't30'], ['t10', 't11', 't12', 't13', 't14', 't15', 't16', 't17', 't7', 't8', 't9'], ['t18', 't19', 't20', 't21', 't22', 't23', 't24', 't25', 't26']]"
        cladeslist.add(Arrays.asList("t1", "t2", "t3", "t4", "t5", "t6"));
        cladeslist.add(Arrays.asList("Z", "t27", "t28", "t29", "t30"));
        cladeslist.add(Arrays.asList("t10", "t11", "t12", "t13", "t14", "t15", "t16", "t17", "t7", "t8", "t9"));
        cladeslist.add(Arrays.asList("t18", "t19", "t20", "t21", "t22", "t23", "t24", "t25", "t26"));
        for (List<String> clade: cladeslist){
            Tuple<Network, Map<NetNode, NetNode>> subnet1 = SuperNetwork3.getSubNetwork(net_original, clade, true);
            netlist.add(subnet1.Item1);
        }
        System.out.println(netlist);

        double [][] matrix = {{0, 3, 6, 6, 6, 7, 9, 9, 9, 10, 10, 7, 7, 8, 9, 9, 9, 9, 9, 9, 7, 9, 9, 9, 10, 10, 7, 8, 9, 10, 10},
                {3, 0, 5, 5, 5, 6, 8, 8, 8, 9, 9, 6, 6, 7, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 9, 9, 6, 7, 8, 9, 9},
                {6, 5, 0, 2, 4, 5, 7, 7, 7, 8, 8, 7, 7, 8, 9, 9, 9, 9, 9, 9, 7, 9, 9, 9, 10, 10, 7, 8, 9, 10, 10},
                {6, 5, 2, 0, 4, 5, 7, 7, 7, 8, 8, 7, 7, 8, 9, 9, 9, 9, 9, 9, 7, 9, 9, 9, 10, 10, 7, 8, 9, 10, 10},
                {6, 5, 4, 4, 0, 3, 5, 5, 5, 6, 6, 7, 7, 8, 9, 9, 9, 9, 9, 9, 7, 9, 9, 9, 10, 10, 7, 8, 9, 10, 10},
                {7, 6, 5, 5, 3, 0, 4, 4, 4, 5, 5, 8, 8, 9, 10, 10, 10, 10, 10, 10, 8, 10, 10, 10, 11, 11, 8, 9, 10, 11, 11}, {9, 8, 7, 7, 5, 4, 0, 2, 4, 5, 5, 10, 10, 11, 12, 12, 12, 12, 12, 12, 10, 12, 12, 12, 13, 13, 10, 11, 12, 13, 13}, {9, 8, 7, 7, 5, 4, 2, 0, 4, 5, 5, 10, 10, 11, 12, 12, 12, 12, 12, 12, 10, 12, 12, 12, 13, 13, 10, 11, 12, 13, 13}, {9, 8, 7, 7, 5, 4, 4, 4, 0, 3, 3, 10, 10, 11, 12, 12, 12, 12, 12, 12, 10, 12, 12, 12, 13, 13, 10, 11, 12, 13, 13}, {10, 9, 8, 8, 6, 5, 5, 5, 3, 0, 2, 11, 11, 12, 13, 13, 13, 13, 13, 13, 11, 13, 13, 13, 14, 14, 11, 12, 13, 14, 14}, {10, 9, 8, 8, 6, 5, 5, 5, 3, 2, 0, 11, 11, 12, 13, 13, 13, 13, 13, 13, 11, 13, 13, 13, 14, 14, 11, 12, 13, 14, 14}, {7, 6, 7, 7, 7, 8, 10, 10, 10, 11, 11, 0, 2, 5, 6, 6, 6, 6, 6, 6, 6, 8, 8, 8, 9, 9, 6, 7, 8, 9, 9}, {7, 6, 7, 7, 7, 8, 10, 10, 10, 11, 11, 2, 0, 5, 6, 6, 6, 6, 6, 6, 6, 8, 8, 8, 9, 9, 6, 7, 8, 9, 9}, {8, 7, 8, 8, 8, 9, 11, 11, 11, 12, 12, 5, 5, 0, 3, 3, 5, 5, 5, 5, 7, 9, 9, 9, 10, 10, 7, 8, 9, 10, 10}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 6, 6, 3, 0, 2, 6, 6, 6, 6, 8, 10, 10, 10, 11, 11, 8, 9, 10, 11, 11}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 6, 6, 3, 2, 0, 6, 6, 6, 6, 8, 10, 10, 10, 11, 11, 8, 9, 10, 11, 11}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 6, 6, 5, 6, 6, 0, 2, 4, 4, 8, 10, 10, 10, 11, 11, 8, 9, 10, 11, 11}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 6, 6, 5, 6, 6, 2, 0, 4, 4, 8, 10, 10, 10, 11, 11, 8, 9, 10, 11, 11}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 6, 6, 5, 6, 6, 4, 4, 0, 2, 8, 10, 10, 10, 11, 11, 8, 9, 10, 11, 11}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 6, 6, 5, 6, 6, 4, 4, 2, 0, 8, 10, 10, 10, 11, 11, 8, 9, 10, 11, 11}, {7, 6, 7, 7, 7, 8, 10, 10, 10, 11, 11, 6, 6, 7, 8, 8, 8, 8, 8, 8, 0, 4, 4, 4, 5, 5, 4, 5, 6, 7, 7}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 8, 8, 9, 10, 10, 10, 10, 10, 10, 4, 0, 2, 4, 5, 5, 6, 7, 8, 9, 9}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 8, 8, 9, 10, 10, 10, 10, 10, 10, 4, 2, 0, 4, 5, 5, 6, 7, 8, 9, 9}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 8, 8, 9, 10, 10, 10, 10, 10, 10, 4, 4, 4, 0, 3, 3, 6, 7, 8, 9, 9}, {10, 9, 10, 10, 10, 11, 13, 13, 13, 14, 14, 9, 9, 10, 11, 11, 11, 11, 11, 11, 5, 5, 5, 3, 0, 2, 7, 8, 9, 10, 10}, {10, 9, 10, 10, 10, 11, 13, 13, 13, 14, 14, 9, 9, 10, 11, 11, 11, 11, 11, 11, 5, 5, 5, 3, 2, 0, 7, 8, 9, 10, 10}, {7, 6, 7, 7, 7, 8, 10, 10, 10, 11, 11, 6, 6, 7, 8, 8, 8, 8, 8, 8, 4, 6, 6, 6, 7, 7, 0, 3, 4, 5, 5}, {8, 7, 8, 8, 8, 9, 11, 11, 11, 12, 12, 7, 7, 8, 9, 9, 9, 9, 9, 9, 5, 7, 7, 7, 8, 8, 3, 0, 3, 4, 4}, {9, 8, 9, 9, 9, 10, 12, 12, 12, 13, 13, 8, 8, 9, 10, 10, 10, 10, 10, 10, 6, 8, 8, 8, 9, 9, 4, 3, 0, 3, 3}, {10, 9, 10, 10, 10, 11, 13, 13, 13, 14, 14, 9, 9, 10, 11, 11, 11, 11, 11, 11, 7, 9, 9, 9, 10, 10, 5, 4, 3, 0, 2}, {10, 9, 10, 10, 10, 11, 13, 13, 13, 14, 14, 9, 9, 10, 11, 11, 11, 11, 11, 11, 7, 9, 9, 9, 10, 10, 5, 4, 3, 2, 0}};
        for (int i = 0; i < matrix.length; i++){
            for (int j =0; j < matrix.length; j++){
                matrix[i][j] *= 2;
            }
        }

//        NetNJMerge merge = new NetNJMerge(matrix, netlist);
//        List<Network> net = merge.mergePairs();
//        System.out.println(net.toString());
    }

    public static void testMatchnetworks(){
        Network net1 = Networks.readNetwork("(t1,(((((t30:0.4)#H1:0.6344::0.6,t29:1.0344):0.4592,(#H1:0.2::0.4,t28:0.6)I0:0.8936)I1:3.8608,t27:5.3544)I2:94.6456,Z:100.0)I3)I5;");
        Network net2 = Networks.readNetwork("((((t6:0.5832,t5:0.5832)I4:0.9116,(t4:0.4)#H1:1.0948::0.5):2.4684,(t3:0.3052,t2:0.3052):3.658):0.7832,(((((((t30)#H2,t29),(#H2,t28)I0)I1,t27)I2,Z)I3,t1)I5,#H1:1.4::0.5):2.9464);");
        List<Network> netlist = matchNetworks(net1, net2);
        System.out.println(netlist);

    }

    public static void testRedundant(){
        Network net1 = Networks.readNetwork("(((t15,((t37,(t25,t14)I3)I10,(t26,(t35,(t3,t23)I0)I2)I9)I11)I14,(t13,(t22,(t20:0.123798681,t21:0.123798681)I4)I12)I15)I16,(t7:4.8857598977,(t38:4.329287445,t32:4.329287445)I7:0.5564724532)I13:0.5007753974)I17;");
        Network net2 = Networks.readNetwork("(Z:100.0,(((((t20,t21)I4,t22)I12,t13:1.700226431)I15:2.92079142,(t15:3.567157584,((t37:3.125263496,(t25:1.836351846,t14:1.836351846)I3:1.28891165)I10:0.1913379973,(t26:3.175296509,(t35:1.833662434,(t3:0.5177220855,t23:0.5177220855)I0:1.315940349)I2:1.341634075)I9:0.1413049845)I11:0.250556091)I14:1.053860267)I16:10.653153017400001,((t39:0.892152881,t12:0.892152881)I6:1.611409191,(t17:2.214562199,(t19:0.1319581023,t30:0.1319581023)I1:2.082604097)I5:0.2889998727)I8:12.7706088):84.7258291316);");
        NetNode nodeA = net1.findNode("I17");
        NetNode nodeB = net2.findNode("Z");
        Map<NetNode, NetNode> redundantList= getRedundantNode(nodeB, nodeA, net2);
        for (NetNode node: redundantList.keySet()){
            System.out.println(node.getName());
        }

    }
    public static void testgetMRCA(){
        Network net1 = Networks.readNetwork("(((t15,((t37,(t25,t14)I3)I10,(t26,(t35,(t3,t23)I0)I2)I9)I11)I14,(t13,(t22,(t20:0.123798681,t21:0.123798681)I4)I12)I15)I16,(t7:4.8857598977,(t38:4.329287445,t32:4.329287445)I7:0.5564724532)I13:0.5007753974)I17;");
        Set<String> set = new HashSet<>();
        set.add("t15");
        set.add("t37");
        set.add("t3");
        NetNode mrca = getMRCA(net1, set);
        System.out.println(mrca.getName());
    }

    public static void testgetMRCAReti(){
        Network net1 = Networks.readNetwork("(((t2,((t63:0.08872438214,t10:0.08872438214)I1:0.2531653096,t7:0.3418896917)I2)I3,(t55:0.3094094558)#H1:0.8276345935::0.4763545006):6.6052797512999994,(#H1:0.06255430952::0.5236454994,t59:0.3719637653)I13:7.3703600349);");
        Networks.autoLabelNodes(net1);
        System.out.println(net1.toString());
        Set<String> set = new HashSet<>();
        set.add("t55");
        set.add("t59");
        NetNode mrca = getMRCA(net1, set);
        System.out.println(mrca.getName());
    }

    public static void testgetNodeFromClade(){
        Set<String> clade = new HashSet<>();
        clade.add("t18");
        clade.add("t27");
        Network net = Networks.readNetwork("(((t2,t17)I1,(((t28:0.7278426847,t16:0.7278426847)I9:0.2254566606)#H1:0.1945802142::0.5341611349,(t11:0.42338387,t34:0.42338387)I8:0.7244956895)I10)I13,(((t38,(t22:2.937865176,(t33:1.479894037)#H2:1.45797114::0.4464888709)I11)I12,((t18:2.434368494,#H1:1.481069149::0.4658388651):0.7537124848,t27:3.188080979):0.2072552739):0.5577594616,#H2:2.473201678::0.5535111291):1.4304966762);");
        Networks.autoLabelNodes(net);
        System.out.println(net.toString());
        NetNode n = getNodeFromClade(net, clade);
        System.out.println(n);
    }

    public static void testgetLeafSetButReti(){
        Network n = Networks.readNetwork("((t16:0.895)#H1:1.76::0.595,((#H1:0.408::0.405,t8:1.303):0.163,(t11:1.333,(t12:0.123,t20:0.123)I2:1.21)I10:0.133)I11:1.189);");
        Set<String> clade = getLeafSetButReti(n.findNode("I11"));
        for (String s: clade){
            System.out.println(s);
        }
    }

    public static void testJoinNodes(){
        Network net1 = Networks.readNetwork("(((((t30:0.1319581023,t19:0.1319581023)I0:2.082604097,t17:2.214562199)I2:0.2889998727,(t12:0.892152881,t39:0.892152881)I1:1.611409191)I3:12.7706088,((((((t23:0.5177220855,t3:0.5177220855)I7:1.315940349,t35:1.833662434)I8:1.341634075,t26:3.175296509)I13:0.1413049845,((t14:1.836351846,t25:1.836351846)I9:1.28891165,t37:3.125263496)I12:0.1913379973)I14:0.250556091,t15:3.567157584)I15:1.053860267,(((t21,t20)I4,t22)I5,t13:1.700226431)I6:2.92079142)I16:10.653153017400001):84.7258291316,Z:100.0);");
    }

    public static void testRerootCompatible(){
        Network net1 = Networks.readNetwork("((t39,t75)I1,((((t15:1.105517021,t16:1.105517021)I3:2.396737825,t25:3.502254847)I5:1.947558176,(t20:4.818472964,(t85:3.88346984,((t29:0.2769542752,t72:0.2769542752)I4,((t80),((t83,t8)I6,(t23,t31)I7)I8)I9)I10)I11:0.9350031231)I12:0.6313400592)I13:94.55018697878,(Z,((t50),t87)I0)I2)I14);");
        Network net2 = Networks.readNetwork("((((t50:0.7284880026)#H1:0.0::0.3735040297,(t39:0.1397969294,t75:0.1397969294)I1:0.5886910732):20.33144073,(((t80:1.348108508)#H2:0.3479510206::0.5465845262,((t83:0.1083187504,t8:0.1083187504)I6:0.6998155581,(t23:0.5923585598,t31:0.5923585598)I7:0.2157757488)I8:0.8879252204)I9,(t29,t72)I4)I10):0.702792199,(#H2:5.752299178::0.4534154738,((#H1:0.0::0.6264959703,t87:0.7284880026)I0,Z)I2):14.662313241);");
        net1.resetRoot("Z");
        net2.resetRoot("Z");
        Set<String> leafset1 = new HashSet<>();
        Set<String> leafset2 = new HashSet<>();
        Set<String> intersection = new HashSet<>();
        for (Object o: net1.getLeaves()){
            NetNode node = (NetNode) o;
            leafset1.add(node.getName());
        }
        for (Object o: net2.getLeaves()){
            NetNode node = (NetNode) o;
            leafset2.add(node.getName());
        }
        intersection.addAll(leafset1);
        intersection.retainAll(leafset2);
        Network subnet1 = SuperNetwork3.getSubNetwork(Networks.readNetwork(net1.toString()), new ArrayList<>(intersection), true).Item1;
        Network subnet2 = SuperNetwork3.getSubNetwork(Networks.readNetwork(net2.toString()), new ArrayList<>(intersection), true).Item1;
        Tuple3<Network, Network, Double> tuple3 = Pipeline.CheckWithTrueNetwork(subnet1, subnet2);
        System.out.println(tuple3.Item1);
        System.out.println(tuple3.Item2);
        System.out.println(tuple3.Item3);

    }


    public static void main(String[] args) {
//        testGetParent();
//        testMapSubNets();
//        testCase1();
//        testMatchNetworks();
//        testCase2();
//        findEdge();
//        testGetNodeLeafMap();
//        testNetworkRooting();
//        testRooting();
//        testCopyNode();
//        testnetwork();
//        testCase3();
//        testMatchnetworks();
//        testRedundant();
        testRerootCompatible();

//        testgetMRCAReti();
//        testgetNodeFromClade();
//        testGetNodeLeafMap();
//        testgetLeafSetButReti();
    }
}
