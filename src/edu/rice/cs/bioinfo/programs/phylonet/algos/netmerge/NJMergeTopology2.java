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
import java.util.stream.StreamSupport;

import com.google.common.collect.Lists;
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


import static edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.util.Utils.*;
import static edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.util.Utils.getAbandonList;
import static edu.rice.cs.bioinfo.programs.phylonet.algos.summarize.majorTree.removeReti;


public class NJMergeTopology2 {
    private List<Network> _subnetworks = new ArrayList<>();
    private final List<Network> _originalSubNetworks = new ArrayList<>();
    private Network _backbone = null;
    private Network _subnet = null;
    private int _num_taxa;
    private double[][] _matrix = new double[_num_taxa][_num_taxa];
    private List<String> _taxonList = new ArrayList<>();
    private List<Set<String>> _leaves = null;
    private String _outgroup = "";
    public static double _epsilon = 0.001;
    private int _nodeCnt = 0;
//    private double _scale = 1;
    private boolean _debug = Utils._debug;

    private Network _backboneTree = null;
    private Network _truenetwork = null;

    public Set<NetNode> _abandonList = new HashSet<>();

    private boolean _addedOutgroup = false;

    /* Constructor */
    public NJMergeTopology2(List<Network> subnetworks, double[][] matrix, List<String> taxonList, String outgroup) {
        setOutgroup(outgroup);
        for (Network net: subnetworks){
            Utils.blankSubNetInternalNodeNames(net);
            this._subnetworks.add(net);
            this._originalSubNetworks.add(net.clone());
        }
        this._matrix = matrix;
        this._taxonList = taxonList;
        _num_taxa = _taxonList.size();
        _leaves = new ArrayList<>();


        if (_subnetworks.get(0).getReticulationCount() > _subnetworks.get(1).getReticulationCount()){
            _backbone = _subnetworks.get(0);
            _subnet = _subnetworks.get(1);
        }
        else{
            _backbone = _subnetworks.get(1);
            _subnet = _subnetworks.get(0);
        }

        checkRoot();

        _subnetworks.set(0, _backbone);
        _subnetworks.set(1, _subnet);
        for (Network net : _subnetworks) {
            Set<String> set = new HashSet<String>();
            for (Object o: net.getLeaves()){
                set.add(((NetNode) o).getName());
            }
            _leaves.add(set);
        }
        List<NetNode> netnodeList = Lists.newArrayList(_subnet.getNetworkNodes());
        removeReti(_subnet, netnodeList, 0);

    }

    public void checkRoot(){
        NetNode rootInT1 = _backbone.findNode(_outgroup);
        NetNode rootInT2 = _subnet.findNode(_outgroup);
        if(rootInT1 == null && rootInT2 == null){

            NetNode<NetNodeInfo> outgroupNode = new BniNetNode();
            outgroupNode.setName(_outgroup);
            NetNode<NetNodeInfo> root = new BniNetNode();
            root.adoptChild(outgroupNode, NetNode.NO_DISTANCE);
            root.adoptChild(_subnet.getRoot(), NetNode.NO_DISTANCE);
            _subnet = new BniNetwork((BniNetNode)root);
            _addedOutgroup = true;
        }
        else{
            _addedOutgroup = false;
        }
    }

    public void setOutgroup(String outgroup){
        this._outgroup = outgroup;
    }


    // A in T1 and B in T2,
    Tuple<List<NetNode>, List<Network>> joinNodesInBothNets(Network net1, NetNode nodeAinT1, Set<String> cladeA,
                                      Network net2, NetNode nodeBinT2, Set<String> cladeB
    ) {
        System.out.println(net1.toString());
        System.out.println(net2.toString());
        if (_debug) {
            System.out.println("joinNodesInBothNets");
        }
//        if (nodeAinT1.getName().equals("I17") && nodeBinT2.getName().equals("Z")){
//            System.out.println("debug point");
//
//        }
        Network net1copy = net1.clone();

        Set<String> leaves1 = Utils.getLeafSet(net1.getRoot());
        Set<String> leaves2 = Utils.getLeafSet(net2.getRoot());
        Set<String> cladeAset = new HashSet<>(cladeA);
        Set<String> cladeBset = new HashSet<>(cladeB);

        List<NetNode> res = new ArrayList<>();
        List<Network> networklist = new ArrayList<>();
        
        NetNode<NetNodeInfo> newNodeToAddinT1 = new BniNetNode();
        newNodeToAddinT1.setName("I"+_nodeCnt);
        if ((cladeA.equals(leaves1) && cladeB.equals(leaves2)) || (cladeA.equals(leaves2) && cladeB.equals(leaves1))) {

            newNodeToAddinT1.adoptChild(net1.getRoot(), NetNode.NO_DISTANCE);
            newNodeToAddinT1.adoptChild(net2.getRoot(), NetNode.NO_DISTANCE);
            newNodeToAddinT1.setData(new NetNodeInfo());
            net1 = new BniNetwork((BniNetNode) newNodeToAddinT1);
            net2 = net1.clone();
            res.add(newNodeToAddinT1);
            res.add(net2.getRoot());
            networklist = Arrays.asList(net1, net2);
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
                res.add(newNodeToAddinT1);

            }
            else {

                Set<String> cladeX = new HashSet<>(leaves1);
                cladeX.removeAll(cladeAset);
                if (_debug) {
                    System.out.println("check point before checkOverlap:"+nodeAinT1.getName()+","+nodeBinT2.getName());
                    System.out.println("cladeX: " + cladeX);
                }
                NetNode newnode = checkOverlap(cladeX, cladeBset, cladeAset, net1, nodeAinT1, net2, nodeBcopy, newNodeToAddinT1, true);
                if(newnode == null){
                    res = Arrays.asList(null, null);
                }
                else{
                    res.add(newnode);
                }

            }
            networklist.add(net1);
            if(_debug){
                System.out.println("after adding node B to T1, net1:"+ net1.toString());
                System.out.println("after adding node B to T1, net2:"+ net2.toString());
            }


            //
            // Add node A to T2
            // (B,X) -> ((A,B),X)), need to check if A and X have overlapped species
            NetNode<NetNodeInfo> newNodeToAddinT2 = new BniNetNode();
            newNodeToAddinT2.setName("I"+_nodeCnt);
            old2new = new HashMap<>();
            NetNode<NetNodeInfo> nodeAcopy = Utils.copyNodeSubnet(nodeAinT1, old2new);
            nodeAcopy.setName(nodeAinT1.getName());
            nodeBinT2.setName(nodeBinT2.getName());
            if (leaves2.equals(cladeBset)) {
                newNodeToAddinT2.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
                newNodeToAddinT2.setData(new NetNodeInfo());
                newNodeToAddinT2.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
                net2 = new BniNetwork((BniNetNode) newNodeToAddinT2);
                System.out.println(newNodeToAddinT2+","+ res);
                if (res.get(0) != null){
                    res.add(newNodeToAddinT2);
                }

            }
            else {//there could be several networks
                Set<String> cladeX = new HashSet<>(leaves2);
                cladeX.removeAll(cladeBset);
                if (_debug) {
                    System.out.println("check point before checkOverlap:"+nodeAinT1.getName()+","+nodeBinT2.getName());
                    System.out.println("cladeX: " + cladeX);
                }
                NetNode newnode  = checkOverlap(cladeX, cladeAset, cladeBset, net2, nodeBinT2, net1copy, nodeAcopy, newNodeToAddinT2, false);
                if(newnode == null){
                    res = Arrays.asList(null, null);
                }
                else if (res.get(0) != null){
                    res.add(newnode);
                }
                else{
                    res = Arrays.asList(null, null);
                }

            }
            networklist.add(net2);
        }
        if (_debug) {
            System.out.println("after adding node A to T2, net1:"+ net1.toString());
            System.out.println("after adding node A to T2, net2:"+ net2.toString());
            System.out.println("finish checkoverlap-----------------");
        }

        return new Tuple<>(res, networklist);
    }


    // node B in net2 (B,X), add node A to net2 ((A,B),X)
    public static NetNode checkOverlap(Set<String> cladeX, Set<String> cladeA, Set<String> cladeB, Network net2, NetNode nodeBinT2, Network net1, NetNode<NetNodeInfo> nodeAcopy, NetNode newNodeToAddinT2, boolean isNetBackbone){
        if (!Collections.disjoint(cladeX, cladeA)) {
            System.out.println("checkOverlap overlap");
            // node under nodeB to node in network


            // intersection of A and X
            Set<String> intersection = new HashSet<>(cladeX);
            intersection.retainAll(cladeA);
            Set<String> unionIntersectionB = new HashSet<>(cladeB);
            unionIntersectionB.addAll(intersection);
            System.out.println("testjoin intersectionB:"+intersection);
            System.out.println("testjoin unionIntersectionB:"+unionIntersectionB);


            Map<NetNode, NetNode> redundantNodeMap = null;
            redundantNodeMap = getRedundantNode(nodeBinT2, nodeAcopy, net2);
//            if(isNetBackbone){
//                redundantNodeMap = getRedundantNode(nodeBinT2, nodeAcopy, net2);
//            }
//            else{//todo
//                redundantNodeMap = getRedundantNode(nodeAcopy, nodeBinT2, net1);
//            }


            NetNode parentAX = null;
            for (NetNode nodeunderA : redundantNodeMap.keySet()) {
//                if(isNetBackbone){
                    parentAX =  (NetNode) nodeunderA.getParents().iterator().next();
                    parentAX.removeChild(nodeunderA);
//                }
//                else{
//
//
//                    NetNode nodeinNet2 = redundantNodeMap.get(nodeunderA);
//                    parentAX =  (NetNode) nodeinNet2.getParents().iterator().next();
//                    parentAX.removeChild(nodeinNet2);
//
//                }

            }

            if(redundantNodeMap.size() == 1){
                NetNode<NetNodeInfo> parentB = (NetNode) nodeBinT2.getParents().iterator().next();
                if (parentB.getName().isEmpty()){
                    for(Object o: parentB.getChildren()){
                        NetNode child = (NetNode) o;
                        if (!child.equals(nodeBinT2)){
                            parentB.removeChild(child);
                            parentB.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
                            if (!parentAX.equals(child)){
                                parentAX.adoptChild(child, NetNode.NO_DISTANCE);
                            }
                            else{
                                System.out.println("problem outlier");
                                NetNode parentparentB = (NetNode) parentB.getParents().iterator().next();
                                parentparentB.removeChild(parentB);
                                NetNode newParentB = new BniNetNode();
                                parentparentB.adoptChild(newParentB, NetNode.NO_DISTANCE);
                                newParentB.adoptChild(child, NetNode.NO_DISTANCE);
                                newParentB.adoptChild(parentB, NetNode.NO_DISTANCE);
                            }

                            parentB.setName(newNodeToAddinT2.getName());
                            return parentB;
                        }
                    }
                }
            }
            else{
                //TODO:test
                System.out.println("redundantNodeMap.size() != 1,"+redundantNodeMap.size());
                for (NetNode nodeunderA : redundantNodeMap.keySet()) {
                    System.out.println("redundantNodeMap:"+nodeunderA.getName());
                }
                return null;
            }


        }
        else{
            NetNodeInfo ni = new NetNodeInfo();
            newNodeToAddinT2.setData(ni);
            newNodeToAddinT2.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
            newNodeToAddinT2.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);

            NetNode<NetNodeInfo> parentB = (NetNode) nodeBinT2.getParents().iterator().next();
            parentB.removeChild(nodeBinT2);
            parentB.adoptChild(newNodeToAddinT2, NetNode.NO_DISTANCE);
            return newNodeToAddinT2;
        }
        return null;
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
        Utils.getNodeLeafMap(net1.getRoot(), node2leaf1, null);
        Utils.getNodeLeafMap(net2.getRoot(), node2leaf2, null);
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
//        newRetic.setParentProbability(newVnode, 0.5);
//        newRetic.setParentDistance(newVnode, 1.0);
//        newRetic.setParentProbability(nodeU, 0.5);
//        newRetic.setParentDistance(nodeU, 1.0);
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

                    if(net2 != null){
//                        System.out.println("mergenet");
                        addInheritanceProb(net2);
                        Network newnet = net2.clone();
//                        System.out.println(net2.toString());
//                        if (newnet == null){
//                            System.out.println("newnet is null");
//                        }
                        candidateNetworks.add(newnet);
//                        System.out.println("added");
                    }

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
        addInheritanceProb(net1);
        addInheritanceProb(net2);
        System.out.println(net1);
        System.out.println(net2);
        if (Utils._debug){
            Networks.autoLabelNodes(net1);
            Networks.autoLabelNodes(net2);
        }


        Map<NetNode, Set<String>> node2leaf1 = new HashMap<>();
        Map<NetNode, Set<String>> node2leaf2 = new HashMap<>();

        List<Network> candidateNetworks = new ArrayList<>();

        Set<String> intersection = Utils.getLeafSet(net1.getRoot());
        intersection.retainAll(Utils.getLeafSet(net2.getRoot()));

        Utils.getNodeLeafMap(net1.getRoot(), node2leaf1, intersection);
        Utils.getNodeLeafMap(net2.getRoot(), node2leaf2, intersection);
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
//                NetNode<NetNodeInfo> nodeUParent = nodeU.getParents().iterator().next();
//                NetNode<NetNodeInfo> nodeVParent = nodeV.getParents().iterator().next();
//                Set<String> nodeVParentLeaf = node2leaf1.get(nodeVParent);
//                Set<String> nodeUParentLeaf = node2leaf1.get(nodeUParent);

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


    public static List<Network> matchNetworks(Network<NetNodeInfo> net1, Network<NetNodeInfo> net2, NetNode retic){

        Networks.removeBinaryNodes(net1);
        Networks.removeBinaryNodes(net2);
        addInheritanceProb(net1);
        addInheritanceProb(net2);
        System.out.println(net1);
        System.out.println(net2);
        if (Utils._debug){
            Networks.autoLabelNodes(net1);
            Networks.autoLabelNodes(net2);
        }


        Map<NetNode, Set<String>> node2leaf1 = new HashMap<>();
        Map<NetNode, Set<String>> node2leaf2 = new HashMap<>();

        List<Network> candidateNetworks = new ArrayList<>();

        Set<String> intersection = Utils.getLeafSet(net1.getRoot());
        intersection.retainAll(Utils.getLeafSet(net2.getRoot()));

        Utils.getNodeLeafMap(net1.getRoot(), node2leaf1, intersection);
        Utils.getNodeLeafMap(net2.getRoot(), node2leaf2, intersection);
        Network subnet1 = SuperNetwork3.getSubNetwork(Networks.readNetwork(net1.toString()), new ArrayList<>(intersection), true).Item1;
        Network subnet2 = SuperNetwork3.getSubNetwork(Networks.readNetwork(net2.toString()), new ArrayList<>(intersection), true).Item1;

        Map<Set<String>, List<NetNode>> set2node2 = Utils.invertMapUsingGroupingBy(node2leaf2);

        if( Networks.hasTheSameTopology(subnet1, subnet2)){
            candidateNetworks.add(net2);
        }
        else{

            Set<String> reticLeafSet = node2leaf1.get(retic);
            Iterator<NetNode<NetNodeInfo>> parentIter = retic.getParents().iterator();
            NetNode<NetNodeInfo> nodeU = parentIter.next();
            NetNode<NetNodeInfo> nodeV = parentIter.next();
            Set<String> nodeUleaf = node2leaf1.get(nodeU);
            Set<String> nodeVleaf = node2leaf1.get(nodeV);
//                NetNode<NetNodeInfo> nodeUParent = nodeU.getParents().iterator().next();
//                NetNode<NetNodeInfo> nodeVParent = nodeV.getParents().iterator().next();
//                Set<String> nodeVParentLeaf = node2leaf1.get(nodeVParent);
//                Set<String> nodeUParentLeaf = node2leaf1.get(nodeUParent);

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
        if (Utils._debug){
            System.out.println("candidate networks after matching two networks are below");
            for (Network net: candidateNetworks){
                System.out.println(net.toString());
            }
        }

        return candidateNetworks;
    }


    //match the reticulations in two networks, add reticulations on net1 to net2
    public static List<Network> match2Networks(Network<NetNodeInfo> net1, Network<NetNodeInfo> net2){

        Networks.removeBinaryNodes(net1);
        Networks.removeBinaryNodes(net2);

        List<NetNode> reticulationList = new ArrayList<>();
        for (NetNode<NetNodeInfo> retic: net1.getNetworkNodes()){
            reticulationList.add(retic);
        }
        //todo permutate the reticulations
        List<Network> candidatenets = new ArrayList<>();
        candidatenets.add(net2);

        for(NetNode<NetNodeInfo> retic: reticulationList){
            List<Network> nextcandidatenets = new ArrayList<>();
            for(Network net: candidatenets){
                nextcandidatenets.addAll(matchNetworks(net1, net, retic));
            }
            candidatenets = nextcandidatenets;
        }

        return candidatenets;
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



    public Tuple3<Boolean, NetNode<NetNodeInfo>, List<Boolean>> joinNodes(NetNode<NetNodeInfo> nodeA, NetNode<NetNodeInfo> nodeB) {
        if (nodeA.getName().equals("I10")&& nodeB.getName().equals("I14")) {
            System.out.println("nodeA or nodeB is check join");
        }
        Set<String> cladeA = Utils.getLeafSet(nodeA);
        Set<String> cladeB = Utils.getLeafSet(nodeB);
        Set<String> cladeAB = new HashSet<>(cladeA);
        cladeAB.addAll(cladeB);

//        List<Boolean> edits = Stream.generate(() -> false).limit(_subnetworks.size()).collect(Collectors.toList());
        boolean edited = false;
        List<Boolean> isChildOfNetNode = new ArrayList<>();
        isChildOfNetNode.add(false);
        isChildOfNetNode.add(false);
        NetNode<NetNodeInfo> newNodeToAdd = new BniNetNode();
        newNodeToAdd.setName("I"+_nodeCnt);
//        cladeA = getLeafSetButReti(nodeAinTs.get(0));
//        cladeB = getLeafSetButReti(nodeBinTs.get(0));
//        if (!Collections.disjoint(cladeA, cladeB)) {
//            throw new Exception("Nodes are not disjoint on their leaf sets! "+nodeA.getName()+":"+cladeA+";"+nodeB.getName()+":"+cladeB);
//        }

        List<NetNode> nodeAinTs = new ArrayList<>();
        List<NetNode> nodeBinTs = new ArrayList<>();
        List<Set<String>> cladeAinTs = new ArrayList<>();
        List<Set<String>> cladeBinTs = new ArrayList<>();

        List<Network> subnetworks = new ArrayList<>();
        subnetworks.add(_backbone);
        subnetworks.add(_subnet);

        for (int i = 0; i < subnetworks.size(); i++) {
            Set<String> leaf = this._leaves.get(i);

            if (cladeA.equals(leaf)) {
                nodeAinTs.add(subnetworks.get(i).getRoot());
            } else {
//                NetNode temp = getNodeFromClade(subnetworks.get(i), cladeA);
                NetNode temp = subnetworks.get(i).findNode(nodeA.getName());
                nodeAinTs.add(temp);

            }

            if (cladeB.equals(leaf)) {
                nodeBinTs.add(subnetworks.get(i).getRoot());
            } else {
//                NetNode temp = getNodeFromClade(subnetworks.get(i), cladeB);
                NetNode temp = subnetworks.get(i).findNode(nodeB.getName());
                nodeBinTs.add(temp);
            }
        }

        List<Boolean> nAinTs = new ArrayList<>();
        for (NetNode nodeAinT : nodeAinTs) {
            nAinTs.add(nodeAinT != null);
            if (nodeAinT != null && !nodeAinT.isRoot()){
                for (Object o: nodeAinT.getParents()){
                    NetNode p = (NetNode) o;
                    if (p.isNetworkNode()){
                        isChildOfNetNode.set(0, true);
                    }
                }
            }
            if (nodeAinT != null){
                cladeAinTs.add(Utils.getLeafSet(nodeAinT));
            }
            else{
                cladeAinTs.add(null);
            }
        }


        List<Boolean> nBinTs = new ArrayList<>();
        for (NetNode nodeBinT : nodeBinTs) {
            nBinTs.add(nodeBinT != null);
            if(nodeBinT != null && !nodeBinT.isRoot()){
                for (Object o: nodeBinT.getParents()){
                    NetNode p = (NetNode) o;
                    if (p.isNetworkNode()){
                        isChildOfNetNode.set(1, true);
                    }
                }
            }
            if (nodeBinT != null){
                cladeBinTs.add(Utils.getLeafSet(nodeBinT));
            }
            else{
                cladeBinTs.add(null);
            }
        }

        int i =0;
        boolean nAinT1 = nAinTs.get(i);
        boolean nBinT1 = nBinTs.get(i);
        NetNode nodeAinT1 = nodeAinTs.get(i);
        NetNode nodeBinT1 = nodeBinTs.get(i);
        Set<String> cladeAinT1 = cladeAinTs.get(i);
        Set<String> cladeBinT1 = cladeBinTs.get(i);
        Set<String> cladeABinT1 = new HashSet<>();
        if (cladeAinT1!= null && cladeBinT1!=null){
            cladeABinT1.addAll(cladeAinT1);
            cladeABinT1.addAll(cladeBinT1);
        }

        int j = 1;
        boolean nAinT2 = nAinTs.get(j);
        boolean nBinT2 = nBinTs.get(j);
        NetNode nodeAinT2 = nodeAinTs.get(j);
        NetNode nodeBinT2 = nodeBinTs.get(j);
        Set<String> cladeAinT2 = cladeAinTs.get(j);
        Set<String> cladeBinT2 = cladeBinTs.get(j);
        Set<String> cladeABinT2 = new HashSet<>();
        if (cladeAinT2!= null && cladeBinT2!=null){
            cladeABinT2.addAll(cladeAinT2);
            cladeABinT2.addAll(cladeBinT2);
        }


        Network n1 = _backbone.clone();
        Network n2 = _subnet.clone();

        if ((nAinT1 || nAinT2) && (nBinT1 || nBinT2)) {
            if (nAinT1 && nAinT2) {
                //nodeA in *both* T1 and T2
                if (nBinT1 && nBinT2) {
                    // Case 1: nodeB in *both* T1 and T2
                    // Valid if nodeA and nodeB are siblings in both T1 & T2
//                    NetNode node1 = getNodeFromClade(_backbone, cladeAB);
//                    NetNode node2 = getNodeFromClade(_subnet, cladeAB);
                    NetNode node1 = getNodeFromClade(_backbone, cladeABinT1);
                    NetNode node2 = getNodeFromClade(_subnet, cladeABinT2);

//                    NetNode node1 = commonParents(nodeAinT1, nodeBinT1, _abandonList);
//                    NetNode node2 = commonParents(nodeAinT2, nodeBinT2, _abandonList);
//                    _abandonList.addAll(getAbandonList(node1));
//                    _abandonList.addAll(getAbandonList(node2));
                    node1.setName("I"+_nodeCnt);
                    node2.setName("I"+_nodeCnt);


                }
                else if (nBinT1) {
//                    //Case 2: node B in T1 only
//                    // Valid if nodeA and nodeB are siblings in T1
                    // add nodeB to T2
//                    NetNode node = getNodeFromClade(_subnetworks.get(i), cladeAB);
//                    if (node == null) {
//                        violates = true;
//                    }
//                    NetNode node = getNodeFromClade(_backbone, cladeAB);
                    NetNode node = getNodeFromClade(_backbone, cladeABinT1);
//                    NetNode node = commonParents(nodeAinT1, nodeBinT1, _abandonList);
                    node.setName("I"+_nodeCnt);
//                    _abandonList.addAll(getAbandonList(node));
//                    joinNodesInOneNet(nodeBinT1, cladeB, _subnet, nodeAinT2, cladeA);
                    joinNodesInOneNet(nodeBinT1, cladeBinT1, _subnet, nodeAinT2, cladeAinT2, false);
                    edited = true;
                }
                else {
                    // Case 3: Node B in T2 only
                    // Valid if nodeA and nodeB are siblings in T2
                    // add node B to T1
//                    NetNode node = getNodeFromClade(_subnet, cladeAB);
                    NetNode node = getNodeFromClade(_subnet, cladeABinT2);
//                    NetNode node = commonParents(nodeAinT2, nodeBinT2, _abandonList);
                    node.setName("I"+_nodeCnt);
//                    _abandonList.addAll(getAbandonList(node));
//                    joinNodesInOneNet(nodeBinT2, cladeB, _backbone, nodeAinT1, cladeA);
                    joinNodesInOneNet(nodeBinT2, cladeBinT2, _backbone, nodeAinT1, cladeAinT1, true);
                    edited = true;

                }
            } else if (nAinT1) {
                // nodeA in T1 only
                if(nBinT1 && nBinT2){
                    // Case 4: nodeB in *both* T1 and T2
                    // add nodeA to T2
//                    NetNode node = getNodeFromClade(_subnetworks.get(i), cladeAB);
//                    if (node == null) {
//                        violates = true;
//                    }
//                    NetNode node = getNodeFromClade(_backbone, cladeAB);
                    NetNode node = getNodeFromClade(_backbone, cladeABinT1);
//                    NetNode node = commonParents(nodeAinT1, nodeBinT1, _abandonList);
                    node.setName("I"+_nodeCnt);
//                    _abandonList.addAll(getAbandonList(node));
//                    joinNodesInOneNet(nodeA, cladeA, _subnet, nodeBinT2, cladeB);
                    joinNodesInOneNet(nodeAinT1, cladeAinT1, _subnet, nodeBinT2, cladeBinT2, false);
                    edited = true;
                }
                else if (nBinT1){
                    // Case 5: Node B in T1 only
                    // Valid if nodeA and nodeB are siblings in T1
//                    NetNode node = getNodeFromClade(_backbone, cladeAB);
//                    if (node == null){
//                        violates = true;
//                    }
//                    NetNode node = getNodeFromClade(_backbone, cladeAB);
                    NetNode node = getNodeFromClade(_backbone, cladeABinT1);
//                    NetNode node = commonParents(nodeAinT1, nodeBinT1, _abandonList);
                    node.setName("I"+_nodeCnt);
//                    _abandonList.addAll(getAbandonList(node));

                }
                else if (nBinT2){
                    // Case 6: Node B in T2 only
                    // Do join in both trees and test for compatibility
//                            Network n1 = _subnetworks.get(i).clone();
//                            Network n2 = _subnetworks.get(j).clone();
//                    NetNode nA = getNodeFromClade(n1, cladeA);
//                    NetNode nB = getNodeFromClade(n2, cladeB);
                    edited = true;
//                    List<NetNode> nodelist = joinNodesInBothNets(_backbone, nodeAinT1, cladeA, _subnet, nodeBinT2, cladeB);
                    Tuple<List<NetNode>, List<Network>> tuple = joinNodesInBothNets(_backbone, nodeAinT1, cladeAinT1, _subnet, nodeBinT2, cladeBinT2);
                    List<NetNode> nodelist = tuple.Item1;
                    List<Network> networklist = tuple.Item2;
//                    _abandonList.addAll(getAbandonList(nodelist.get(0)));
//                    _abandonList.addAll(getAbandonList(nodelist.get(1)));
//                    List<Network> networklist = joinNodesInBothNets(n1, nA, cladeA, n2, nB, cladeB);
//                    n1 = networklist.get(0);
//                    n2 = networklist.get(1);
                    _backbone = networklist.get(0);
                    _subnet = networklist.get(1);
                    _subnetworks = networklist;

                }

            } else if (nAinT2) {
                //nodeA in T2 only
                if (nBinT1 && nBinT2){
                    //Case 7: nodeB in and T2
                    // add nodeA to T1
//                    NetNode node = getNodeFromClade(_subnet, cladeAB);
                    NetNode node = getNodeFromClade(_subnet, cladeABinT2);

//                    NetNode node = commonParents(nodeAinT2, nodeBinT2, _abandonList);
                    node.setName("I"+_nodeCnt);

//                    joinNodesInOneNet(nodeAinT2, cladeA, _backbone, nodeBinT1, cladeB);
                    joinNodesInOneNet(nodeAinT2, cladeAinT2, _backbone, nodeBinT1, cladeBinT1, true);
                    edited = true;
                }
                else if (nBinT1){
                    //Case 8 (reverse of Case 6): Node B in T1 only
//                    Network n1 = _subnetworks.get(i).clone();
//                    Network n2 = _subnetworks.get(j).clone();

//                    NetNode nB = getNodeFromClade(n1, cladeB);
//                    NetNode nA = getNodeFromClade(n2, cladeA);
//                    List<Network> networklist = joinNodesInBothNets(n1, nB, cladeB, n2, nA, cladeA);
//                    List<NetNode> nodelist = joinNodesInBothNets(_backbone, nodeBinT1, cladeB, _subnet, nodeAinT2, cladeA);
                    Tuple<List<NetNode>, List<Network>> tuple = joinNodesInBothNets(_backbone, nodeBinT1, cladeBinT1, _subnet, nodeAinT2, cladeAinT2);
                    List<NetNode> nodelist = tuple.Item1;
                    List<Network> networklist = tuple.Item2;

                    _backbone = networklist.get(0);
                    _subnet = networklist.get(1);
                    _subnetworks = networklist;

//                    _abandonList.addAll(getAbandonList(nodelist.get(0)));
//                    _abandonList.addAll(getAbandonList(nodelist.get(1)));
                    edited = true;

                }
                else if (nBinT2){
                    // Case 9: Node B in T2 only
                    // Only valid if (nodeA, nodeB) are siblings in T2
                    //do nothing
//                    NetNode node = getNodeFromClade(_subnet, cladeAB);
                    NetNode node = getNodeFromClade(_subnet, cladeABinT2);
//                    NetNode node = commonParents(nodeAinT2, nodeBinT2, _abandonList);
                    node.setName("I"+_nodeCnt);

                }

            }
        }
        NetNode newNode = new BniNetNode();
        newNode.setData(new NetNodeInfo());
        newNode.setName("I"+_nodeCnt);
        Map<NetNode, NetNode> old2new = new HashMap<>();
        NetNode c1 = Utils.copyNodeSubnet(nodeA, old2new);
        NetNode c2 = Utils.copyNodeSubnet(nodeB, old2new);
        newNode.adoptChild(c1, NetNode.NO_DISTANCE);
        newNode.adoptChild(c2, NetNode.NO_DISTANCE);
        return new Tuple3<>(edited, newNode, isChildOfNetNode);
    }

    public NetNode joinNodesInOneNet(NetNode nodeAinT1,Set<String> cladeA, Network net2, NetNode nodeBinT2, Set<String> cladeB, boolean isnetBackbone) {
        /*
        Join clade A and clade B in just one tree

        1. tree 1 as (A,X); and extract (A,X);
        2. tree 2 as (B,Y); and extract (B,Y);
        3. attach A to new tree 2 as ((A,B),Y);
        */
        if (_debug) {
            System.out.println("joinNodesInOneNet");
        }




        if(!Collections.disjoint(cladeA, cladeB)) {
            Set<String> intersection = new HashSet<>(cladeA);
            intersection.retainAll(cladeB);
            NetNode mrca = getNodeFromClade(_subnet, intersection);
            if(mrca!=null){
                NetNode parent  = (NetNode)mrca.getParents().iterator().next();
                parent.removeChild(mrca);
            }
            else{
                for(String leafname :intersection){
                    NetNode node2rm = _subnet.findNode(leafname);
                    NetNode parent = (NetNode) node2rm.getParents().iterator().next();
                    parent.removeChild(_subnet.findNode(leafname));
                }
            }
            if(isnetBackbone){
                cladeB.removeAll(intersection);
            }
            else {
                cladeA.removeAll(intersection);
            }


        }

        Map<NetNode, NetNode> old2new = new HashMap<>();
        NetNode<NetNodeInfo> nodeAcopy = Utils.copyNodeSubnet(nodeAinT1, old2new);
        nodeAcopy.setName(nodeAinT1.getName());

        Set<String> leaf = Utils.getLeafSet(net2.getRoot());
        NetNode<NetNodeInfo> newNodeToAdd = new BniNetNode();
        newNodeToAdd.setName("I" + _nodeCnt);
        Set<String> cladeX = new HashSet<>(leaf);
        cladeX.removeAll(cladeB);

        if (leaf.equals(cladeB)) {
            newNodeToAdd.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
            newNodeToAdd.setData(new NetNodeInfo());
            newNodeToAdd.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
            net2 = new BniNetwork((BniNetNode) newNodeToAdd);
            return newNodeToAdd;
        }
        else {

            if (_debug) {
                System.out.println("check point before checkOverlap:"+nodeAinT1.getName()+","+nodeBinT2.getName());
                System.out.println("cladeX: " + cladeX);
            }
            NetNode newNodeToAddinT2 = new BniNetNode();
            newNodeToAddinT2.setName("I" + _nodeCnt);
            NetNode newnode  = checkOverlap(cladeX, cladeA, cladeB, net2, nodeBinT2, _backbone, nodeAcopy, newNodeToAddinT2, isnetBackbone);

            return newnode;

        }


    }


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

    public static void findAncestors(NetNode node, Set<NetNode> ancestors){
        for (Object parent: node.getParents()){
            NetNode parentnode = (NetNode) parent;
            ancestors.add(parentnode);
            findAncestors(parentnode, ancestors);
        }
    }

    public static List<NetNode> findMRCA(Network net, Set<String> clade){
        Map<NetNode, Set<NetNode>> ancestorsMap = new HashMap();
        Set<NetNode> commonAncestors = new HashSet<>();
        boolean begin = false;
        for(String leafname: clade) {
            NetNode node = net.findNode(leafname);
            if(node == null) {
                throw new IllegalArgumentException("Leaf not found in the network: " + leafname);
            }
            Set<NetNode> ancestors = new HashSet<>();
            findAncestors(node, ancestors);
            ancestorsMap.put(node, ancestors);
            if (!begin) {
                commonAncestors.addAll(ancestors);
                begin = true;
            } else {
                commonAncestors.retainAll(ancestors);
            }
        }
        List<NetNode> MRCALsit = new ArrayList<>();
        for (NetNode ca: commonAncestors){
            boolean isMostRecent = true;
            for (Object o: ca.getChildren()){
                NetNode child = (NetNode) o;
                if (commonAncestors.contains(child)){
                    isMostRecent = false;
                    break;
                }
            }
            if (isMostRecent){
                MRCALsit.add(ca);
            }
        }
        return MRCALsit;
    }


    public static NetNode getNodeFromClade(Network net, Set<String> clade) {
        // Returns the MRCA node of the clade

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
    public double areNetsCompatible(Network net1, Network net2){
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

        if (_debug){
            System.out.println("areNetsCompatible:");
            System.out.println(net1.toString());
            System.out.println(net2.toString());
        }

        if (_debug){
            System.out.println(intersection);
        }

//        if (intersection.size() <= 1){
//            return true;
//        }

        if (intersection.contains(this._outgroup)){
            NetNode node1 = net1.findNode(this._outgroup);
            if(isAtReticulationNode(node1)){
                return _compatible_max_value;
            }

            NetNode node2 = net2.findNode(this._outgroup);
            if(isAtReticulationNode(node2)){
                return _compatible_max_value;
            }

            net1.resetRoot(this._outgroup);
            net2.resetRoot(this._outgroup);



        }
        Network subnet1 = SuperNetwork3.getSubNetwork(net1, new ArrayList<>(intersection), true).Item1;
        Network subnet2 = SuperNetwork3.getSubNetwork(net2, new ArrayList<>(intersection), true).Item1;

        addInheritanceProb(subnet1);
        addInheritanceProb(subnet2);
        System.out.println(subnet1.toString());
        System.out.println(subnet2.toString());
        Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(subnet1, subnet2);
//        if (Networks.hasTheSameTopology(subnet1.Item1, subnet2.Item1)){
//            return true;
//        }
//        if (closest.Item3 < 0.01){
//            if(_debug){
//                System.out.println("Compatible");
//            }
//
//            return true;
//        }
//        else{
//            if (_debug){
//                System.out.println("Incompatible:"+closest.Item3);
//                System.out.println(closest.Item1.toString());
//                System.out.println(closest.Item2.toString());
//
//            }
//
//            return false;
//        }
        return closest.Item3;
    }

    /*
        Test whether joining cladeA and cladeB in one
        or both networks causes the two networks to be incompatible
    */
    public double testJoin(NetNode nodeA, NetNode nodeB) throws Exception {
        Set<String> cladeA = Utils.getLeafSet(nodeA);
        Set<String> cladeB = Utils.getLeafSet(nodeB);
        Set<String> cladeAB = new HashSet<>(cladeA);
        cladeAB.addAll(cladeB);
        if(nodeA.getName().equals("I16")&&nodeB.getName().equals("I27")){
            System.out.println("check point");
        }

//        cladeA = getLeafSetButReti(nodeAinTs.get(0));
//        cladeB = getLeafSetButReti(nodeBinTs.get(0));
//        if (!Collections.disjoint(cladeA, cladeB)) {
//            throw new Exception("Nodes are not disjoint on their leaf sets! "+nodeA.getName()+":"+cladeA+";"+nodeB.getName()+":"+cladeB);
//        }

        List<NetNode> nodeAinTs = new ArrayList<>();
        List<NetNode> nodeBinTs = new ArrayList<>();
        List<Set<String>> cladeAinTs = new ArrayList<>();
        List<Set<String>> cladeBinTs = new ArrayList<>();

        List<Network> subnetworks = new ArrayList<>();
        subnetworks.add(_backbone);
        subnetworks.add(_subnet);
        for (int i = 0; i < subnetworks.size(); i++) {
            Set<String> leaf = this._leaves.get(i);

            if (cladeA.equals(leaf)) {
                nodeAinTs.add(nodeA);
            } else {
//                NetNode temp = getNodeFromClade(subnetworks.get(i), cladeA);
                NetNode temp = subnetworks.get(i).findNode(nodeA.getName());
                nodeAinTs.add(temp);
            }

            if (cladeB.equals(leaf)) {
                nodeBinTs.add(nodeB);
            } else {
//                NetNode temp = getNodeFromClade(subnetworks.get(i), cladeB);
                NetNode temp = subnetworks.get(i).findNode(nodeB.getName());
                nodeBinTs.add(temp);
            }
        }

        List<Boolean> nAinTs = new ArrayList<>();
        for (NetNode nodeAinT : nodeAinTs) {
            nAinTs.add(nodeAinT != null);
            if (nodeAinT != null){
                cladeAinTs.add(Utils.getLeafSet(nodeAinT));
            }
            else{
                cladeAinTs.add(null);
            }
        }
        if (Collections.frequency(nAinTs, true) < 1) {
            throw new Exception("Node A: "+nodeA.getName()+" was not found in any tree!");
        }

        List<Boolean> nBinTs = new ArrayList<>();
        for (NetNode nodeBinT : nodeBinTs) {
            nBinTs.add(nodeBinT != null);
            if (nodeBinT != null){
                cladeBinTs.add(Utils.getLeafSet(nodeBinT));
            }
            else{
                cladeBinTs.add(null);
            }
        }


        if (Collections.frequency(nBinTs, true) < 1) {
            throw new Exception("Node B: "+nodeB.getName()+" was not found in any tree!--"+cladeB);
        }



        double violates = 0.0;
        int i =0;
        boolean nAinT1 = nAinTs.get(i);
        boolean nBinT1 = nBinTs.get(i);
        NetNode nodeAinT1 = nodeAinTs.get(i);
        NetNode nodeBinT1 = nodeBinTs.get(i);
        Set<String> cladeAinT1 = cladeAinTs.get(i);
        Set<String> cladeBinT1 = cladeBinTs.get(i);
        Set<String> cladeABinT1 = new HashSet<>();
        if (cladeAinT1!= null && cladeBinT1!=null){
            cladeABinT1.addAll(cladeAinT1);
            cladeABinT1.addAll(cladeBinT1);
        }


        int j = 1;
        boolean nAinT2 = nAinTs.get(j);
        boolean nBinT2 = nBinTs.get(j);
        NetNode nodeAinT2 = nodeAinTs.get(j);
        NetNode nodeBinT2 = nodeBinTs.get(j);
        Set<String> cladeAinT2 = cladeAinTs.get(j);
        Set<String> cladeBinT2 = cladeBinTs.get(j);
        Set<String> cladeABinT2 = new HashSet<>();
        if (cladeAinT2!= null && cladeBinT2!=null){
            cladeABinT2.addAll(cladeAinT2);
            cladeABinT2.addAll(cladeBinT2);
        }



        Network n1 = _backbone.clone();
        Network n2 = _subnet.clone();

        if ((nAinT1 || nAinT2) && (nBinT1 || nBinT2)) {
            if (nAinT1 && nAinT2) {
                //nodeA in *both* T1 and T2
                if (nBinT1 && nBinT2) {
                    // Case 1: nodeB in *both* T1 and T2
                    // Valid if nodeA and nodeB are siblings in both T1 & T2
//                    NetNode node1 = getNodeFromClade(_backbone, cladeAB);
//                    NetNode node2 = getNodeFromClade(_subnet, cladeAB);

                    NetNode node1 = getNodeFromClade(_backbone, cladeABinT1);
                    NetNode node2 = getNodeFromClade(_subnet, cladeABinT2);


//                    NetNode node1 = commonParents(nodeAinT1, nodeBinT1, _abandonList);
//                    NetNode node2 = commonParents(nodeAinT2, nodeBinT2, _abandonList);

                    if (node1 == null || node2 == null) {
                        violates = _compatible_max_value;
                    }
                }
                else if (nBinT1) {
//                    //Case 2: node B in T1 only
//                    // Valid if nodeA and nodeB are siblings in T1
//                    NetNode node = getNodeFromClade(_backbone, cladeAB);
                    NetNode node = getNodeFromClade(_backbone, cladeABinT1);

//                    NetNode node = commonParents(nodeAinT1, nodeBinT1, _abandonList);
                    if (node == null) {
                        violates = _compatible_max_value;
                    }
                }
                else {
                    // Case 3: Node B in T2 only
                    // Valid if nodeA and nodeB are siblings in T2
//                    NetNode node = getNodeFromClade(_subnet, cladeAB);
                    NetNode node = getNodeFromClade(_subnet, cladeABinT2);
//                    NetNode node = commonParents(nodeAinT2, nodeBinT2, _abandonList);
                    if (node == null) {
                        violates = _compatible_max_value;
                    }
                }
            } else if (nAinT1) {
                // nodeA in T1 only
                if(nBinT1 && nBinT2){
                    // Case 4: nodeB in *both* T1 and T2
//                    NetNode node = getNodeFromClade(_backbone, cladeAB);
                    NetNode node = getNodeFromClade(_backbone, cladeABinT1);
//                    NetNode node = commonParents(nodeAinT1, nodeBinT1, _abandonList);
                    if (node == null) {
                        violates = _compatible_max_value;
                    }
                }
                else if (nBinT1){
                    // Case 5: Node B in T1 only
                    // Valid if nodeA and nodeB are siblings in T1
//                    NetNode node = getNodeFromClade(_backbone, cladeAB);
                    NetNode node = getNodeFromClade(_backbone, cladeABinT1);
//                    NetNode node = commonParents(nodeAinT1, nodeBinT1, _abandonList);
                    if (node == null){
                        violates = _compatible_max_value;
                    }
                }
                else if (nBinT2){
                    // Case 6: Node B in T2 only
                    // Do join in both trees and test for compatibility
//                            Network n1 = _subnetworks.get(i).clone();
//                            Network n2 = _subnetworks.get(j).clone();
//                    NetNode nA = getNodeFromClade(n1, cladeA);
//                    NetNode nB = getNodeFromClade(n2, cladeB);
                    NetNode nA = n1.findNode(nodeA.getName());
                    NetNode nB = n2.findNode(nodeB.getName());
//                    List<NetNode> netnodelist = joinNodesInBothNets(n1, nA, cladeA, n2, nB, cladeB);
                    Tuple<List<NetNode>, List<Network>> tuple = joinNodesInBothNets(n1, nA, cladeAinT1, n2, nB, cladeBinT2);
                    List<NetNode> netnodelist = tuple.Item1;
                    List<Network> networklist = tuple.Item2;
                    n1 = networklist.get(0);
                    n2 = networklist.get(1);
                    NetNode node1 = netnodelist.get(0);
                    NetNode node2 = netnodelist.get(1);
                    if (node1 == null && node2 == null){
                        violates = _compatible_max_value;
                    }
                    else if (node1 !=  null){
                        violates = areNetsCompatible(n1, n2);
                    }
                }
                else{
                    throw new Exception("Node B not found in either network!\n");
                }

            } else if (nAinT2) {
                //nodeA in T2 only
                if (nBinT1 && nBinT2){
                    //Case 7: nodeB in T1 and T2
//                    NetNode node = getNodeFromClade(_subnet, cladeAB);
                    NetNode node = getNodeFromClade(_subnet, cladeABinT2);
//                    NetNode node = commonParents(nodeAinT2, nodeBinT2, _abandonList);
                    if (node == null){
                        violates = _compatible_max_value;
                    }

                }
                else if (nBinT1){
                    //Case 8 (reverse of Case 6): Node B in T1 only
//                    Network n1 = _subnetworks.get(i).clone();
//                    Network n2 = _subnetworks.get(j).clone();

//                    NetNode nB = getNodeFromClade(n1, cladeB);
//                    NetNode nA = getNodeFromClade(n2, cladeA);
                    NetNode nB = n1.findNode(nodeB.getName());
                    NetNode nA = n2.findNode(nodeA.getName());

                    Tuple<List<NetNode>, List<Network>> tuple = joinNodesInBothNets(n1, nB, cladeBinT1, n2, nA, cladeAinT2);
                    List<NetNode> netnodelist = tuple.Item1;
                    List<Network> networklist = tuple.Item2;
                    n1 = networklist.get(0);
                    n2 = networklist.get(1);
                    addInheritanceProb(n1);
                    addInheritanceProb(n2);
//                    List<NetNode> netnodelist = joinNodesInBothNets(n1, nB, cladeB, n2, nA, cladeA);
                    System.out.println(netnodelist);
                    NetNode node1 = netnodelist.get(0);
                    NetNode node2 = netnodelist.get(1);
                    if (node1 == null && node2 == null){
                        violates = _compatible_max_value;
                    }
                    else if (node1 !=  null){
                        addInheritanceProb(n1);
                        addInheritanceProb(n2);
                        violates = areNetsCompatible(n1, n2);
                    }

                }
                else if (nBinT2){
                    // Case 9: Node B in T2 only
                    // Only valid if (nodeA, nodeB) are siblings in T2
//                    NetNode node = getNodeFromClade(_subnet, cladeAB);
                    NetNode node = getNodeFromClade(_subnet, cladeABinT2);
//                    NetNode node = commonParents(nodeAinT2, nodeBinT2, _abandonList);
                    if (node ==  null){
                        violates = _compatible_max_value;
                    }
                }

                else{
                    throw new Exception("Node B not found in either network!\n");
                }

            } else {
                throw new Exception("Node A not found in either tree!");
            }


        }

        return violates;
    }



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
//            if (n == 7){
//                System.out.println("debug");
//            }

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
                    pairs.add(new Tuple<>(idx1, idx2));
                    if(Utils._USEQ){
                        double v1 = (n - 2) * nd1.getData().getNjDistances().get(nd2.getName());
                        double qvalue = v1 - nd1.getData().getNJXsub() - nd2.getData().getNJXsub();

                        qvalues.add(qvalue);
                    }
                    else{
                        qvalues.add(nd1.getData().getNjDistances().get(nd2.getName()));
                    }

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
                System.out.println(_backbone.toString());
                System.out.println(_subnet.toString());
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
            double min_violates = Double.MAX_VALUE;
            for (int idxq : sortedIndices) {
                Tuple<Integer, Integer> pair = pairs.get(idxq);
                NetNode nd1 = nodePool.get(pair.Item1);
                NetNode nd2 = nodePool.get(pair.Item2);

                // Check join does not violate a constraint subnetwork!
                try{
                    if(_debug){
                        System.out.println("testjoin:"+nd1.getName()+";"+nd2.getName());
                    }
                    double violates = testJoin(nd1, nd2);
//                    System.out.println();
//                    System.out.println();
                    if (violates < 0.01) {
                        nodesToJoin = new Tuple<>(nd1, nd2);
                        break;
                    }
                    else if (violates < min_violates){
                        min_violates = violates;
                        nodesToJoin = new Tuple<>(nd1, nd2);
                    }
                }catch (Exception e){
                    e.printStackTrace();
                    System.err.println("Error testJoin(nd1, nd2)!");
                    System.exit(-1);
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


//            Tuple<Boolean, NetNode<NetNodeInfo>> tuple = joinNodes(nodesToJoin.Item1, nodesToJoin.Item2); // TODO: Add required parameters
//            Tuple<List<Boolean>, NetNode<NetNodeInfo>> tuple = joinNodes(_subnetworks.get(0), _subnetworks.get(1), nodesToJoin.Item1, nodesToJoin.Item2);
            Tuple3<Boolean, NetNode<NetNodeInfo>, List<Boolean>> tuple = joinNodes(nodesToJoin.Item1, nodesToJoin.Item2);
            if (_debug){
                System.out.println("after:");
                System.out.println(_backbone.toString());
                System.out.println(_subnet.toString());
                System.out.println("--------------");
            }

            boolean edited = tuple.Item1;
            NetNode<NetNodeInfo> newNode = tuple.Item2;
            List<Boolean> isChildOfNetworkNode = tuple.Item3;

            if (edited) {
                removeInternalLeaf(_subnetworks);
                for (int index = 0; index < _subnetworks.size(); index++) {
                    Network net = _subnetworks.get(index);

                    // Check to see if you can quit early
                    Set<String> leafSet = new HashSet<>();
                    for(Object o: net.getLeaves()){
                        NetNode node = (NetNode) o;
                        if(node.getName().startsWith("I")){
                            System.out.println("node name error");
                            System.out.println(net.toString());
                        }
                        leafSet.add(node.getName());

                    }
                    _leaves.set(index, leafSet);

                    //todo: resetroot
//                    NetNode outGroupNode = net.findNode(_outgroup);
//                    if(outGroupNode != null){
//                        net.resetRoot(_outgroup);
//                    }

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
            }
            if (_debug){
                System.out.println("leaves:"+edited);
                System.out.println(_leaves.get(0));
                System.out.println(_leaves.get(1));
            }


            if (mergedNetworkID >= 0){
                break;
            }
            // Create the new node
//            NetNode<NetNodeInfo> newNode = new BniNetNode();
            nodePool.remove(nodesToJoin.Item1);
            nodePool.remove(nodesToJoin.Item2);
//            if(!isChildOfNetworkNode.get(0)){
//                nodePool.remove(nodesToJoin.Item1);
//
//            }
//            if (!isChildOfNetworkNode.get(0)){
//                nodePool.remove(nodesToJoin.Item2);
//
//            }

//            nodePool.remove(nodesToJoin.Item2);

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

        removeInternalLeaf(_subnetworks);
        for (Network net: _subnetworks) {
            NetNode outgroupNode = net.findNode(_outgroup);
            System.out.println(net.toString());
            if (outgroupNode != null) {
                removeRetiOverOutgroup(outgroupNode);
                net.resetRoot(_outgroup);
                System.out.println(net.toString());
                if (_addedOutgroup) {
                    net.getRoot().removeChild(outgroupNode);
                    Networks.removeBinaryNodes(net);
                }
            }
            System.out.println(net.toString());
        }

//        Queue<Tuple<Network, Integer>> netQueue = new LinkedList<>();
        List<Network> candidateNetworks = new ArrayList<>();
        if(mergedNetworkID < 0){
            System.out.println(_subnetworks.get(0).toString());
            System.out.println(_subnetworks.get(1).toString());
            System.out.println(_taxonList);
            System.out.println(_leaves.get(0));
            System.out.println(_leaves.get(1));
            Set<String> leaves = new HashSet<>(_taxonList);

            System.out.println(_leaves.get(0).equals(new HashSet<>(_taxonList)));
            System.out.println(_leaves.get(1).equals(new HashSet<>(_taxonList)));
            System.out.println(_leaves.get(1).size()+","+_leaves.get(0).size()+","+_taxonList.size());
            _leaves.get(1).removeAll(_taxonList);
            System.out.println(_leaves.get(1));

            System.out.println("leaves:");
            System.out.println(leaves);
        }
        Network candidateNetwork = _subnetworks.get(mergedNetworkID);
        addInheritanceProb(candidateNetwork);
        System.out.println(candidateNetwork.toString());
        candidateNetworks.add(candidateNetwork);
        List<Network> backboneList2 = Pipeline.getAllBackboneNets(candidateNetwork, Integer.MAX_VALUE);
        backboneList2.add(candidateNetwork);
//        int level = 0;
//        netQueue.add(new Tuple<>(candidateNetwork, level));
//        int maxLevel = Integer.max(_subnetworks.get(0).getReticulationCount(), _subnetworks.get(1).getReticulationCount());
//        while(!netQueue.isEmpty() && level <= maxLevel){
//            Tuple<Network, Integer> tuple = netQueue.poll();
//            Network net2match = tuple.Item1;
//            level = tuple.Item2;
//            candidateNetworks.add(net2match);
//            for (Network originalSubnet: _originalSubNetworks){
//                List<Network> matchedlist = matchNetworks(originalSubnet, net2match);
//
//                for (Network net: matchedlist){
//                    boolean newnet = true;
//                    for (Network net2: candidateNetworks){
//                        if (Networks.hasTheSameTopology(net, net2)){
//                            newnet = false;
//                            break;
//                        }
//                    }
//                    if (newnet){
//                        netQueue.add(new Tuple<>(net, level+1));
//                    }
//                }
////                candidateNetworks.addAll(matchedlist);
//            }
//        }



        for (Network net: backboneList2){
            for (Network originalSubnet: _originalSubNetworks){
//                System.out.println("match networks:");
//
//                System.out.println(originalSubnet.toString());
//                System.out.println(net.toString());
                List<Network> matchedlist = match2Networks(originalSubnet, net);
                candidateNetworks.addAll(matchedlist);
            }
        }
//        candidateNetworks.addAll(matchNetworks(_originalSubNetworks.get(0), _subnetworks.get(mergedNetworkID)));
//        candidateNetworks.addAll(matchNetworks(_originalSubNetworks.get(1), _subnetworks.get(mergedNetworkID)));
//        if (Networks.hasTheSameTopology(_backbone, _subnet)){
//            //|| _subnetworks.get(0).getReticulationCount() == 0
//            candidateNetworks.add(_subnetworks.get(0));
//        }
//        else{
//            for (int i = 0; i < _subnetworks.size(); i++) {
//                if ( i == mergedNetworkID) {
//                    continue;
//                }
//                Network net2 = _subnetworks.get(i);
//                if (candidateNetworks.size() == 0){
//                    Network net1 = _subnetworks.get(mergedNetworkID);
//                    if (_debug) {
//                        System.out.println("match networks:");
//                        System.out.println(net2.toString());
//                        System.out.println(net1.toString());
//                    }
////                    System.out.println(net2.toString());
////                    System.out.println(net1.toString());
////                    net2.resetRoot(this._outgroup);
////                    net2.resetRoot(this._outgroup);
//                    candidateNetworks = matchNetworks(net2, net1);
//                    candidateNetworks.addAll(matchNetworks(_originalSubNetworks.get(i).clone(), net1));
//                    if (_debug) {
//                        System.out.println(candidateNetworks);
//                    }
//
//                }
//
//            }
//        }
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
//        System.out.println("start checking");
        for (Network resnet: resultnetwork){
            System.out.println(resnet.toString());

        }
        System.out.println(scorelist.get(scoreIndexes[0]));

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
//    public static NetNode commonParents(NetNode nodeA, NetNode nodeB) {
//        for(Object o1: nodeA.getParents()) {
//            NetNode parentA = (NetNode) o1;
//
//            for (Object o2: nodeB.getParents()) {
//                NetNode parentB = (NetNode) o2;
//                if (parentA.equals(parentB)) {
//                    return parentA;
//                }
//                if (parentB.isNetworkNode()){
//
//                    NetNode commonParent = commonParents(nodeA, parentB);
//                    if (commonParent != null){
//                        return commonParent;
//
//                    }
//                }
//            }
//            if (parentA.isNetworkNode()){
//
//                NetNode commonParent = commonParents(parentA, nodeB);
//                if (commonParent != null){
//                    return commonParent;
//                }
//
//            }
//        }
//        return null;
//    }
    public static void test(){
        Network net1 = Networks.readNetwork("(((t38,(((t5:1.070155914,t54:1.070155914)I29:0.4144493305)#H1:0.2519861945::0.5148253302,(t127:0.7949493714,t49:0.7949493714)I28:0.9416420672)I30)I31,(((t160:0.2027410277,t46:0.2027410277)I20:1.213592991)#H2:0.7276044754::0.6194566142,t57:2.143938494)I26:4.155281518):2.130819673,(((((t6,t88)I21,(t124,t48)I22)I23,t151)I24,((#H2:0.4875205662::0.3805433858,((t20:0.03432457825,t111:0.03432457825)I17:0.7883145163,t50:0.8226390945)I19:1.081215491):2.270330631,(t153:0.7514288513,t164:0.7514288513)I18:3.422756365)I25)I27,#H1:2.854063559::0.4851746698):4.091370882);" );
        Network net2 = Networks.readNetwork("((((((((((t160,t46)I20),((t20,t111)I17,t50)I19),(t153,t164)I18)I25,(((t6:0.8751778875,t88:0.8751778875)I21:0.6139868211,(t124:0.879112813,t48:0.879112813)I22:0.6100518956)I23:0.9900317328,t151:2.479196441)I24)I27,(t30:3.423483102,(((((t5,t54)I29),(t127,t49)I28)I30,t38)I31,t136:2.453060145)I32:0.9704229574)I33:5.006556583):3.712209664,(t18:6.627727187,t19:6.627727187)I16:5.514522162):2.103999998,t114:14.246249343):3.855278046,((t102:5.2134814700000005,(((((t72:0.506068326,t132:0.506068326)I4:0.8684131795,(((t78:0.001949521357,t13:0.001949521357)I2:0.4882188782,t125:0.4901683996)I5:0.5462996732,t86:1.036468073)I7:0.3380134327)I8:0.1638615484,(t142:0.5976462728,(t63:0.2278918387,t27:0.2278918387)I3:0.3697544341)I6:0.9406967811)I9:2.049656702,(t117:2.189447962,(t113:0.4804166958,t130:0.4804166958)I1:1.709031266)I10:1.3985517943999999)I12:0.8361047715,(t34:2.536646868,t152:2.536646868)I11:1.88745766)I13:0.7893769417)I14:1.094675071,(t37:0.2528462507,t105:0.2528462507)I0:6.05531029)I15:11.79337085):81.898472607,Z:100.0);");
        System.out.println(net1.findNode("I27"));
        System.out.println(net2.findNode("I34"));
        NetNode p = net1.findNode("I27");
        NetNode c = net1.findNode("I20");

        NetNode oc = p.getTheOtherChildren(c);
        System.out.println("otherchild:");
        System.out.println(oc.getName());


        Networks.autoLabelNodes(net1);
        NetNode n1 = net1.findNode("I29");
        NetNode n2 = net1.findNode("I27");

        System.out.println(net1.toString());
//        NetNode cp = commonParents(n1, n2);
//        System.out.println(cp.getName());
        Set<String> clades = new HashSet<>();
//        clades.add("t160");
//        clades.add("t57");
//        for (NetNode x : findMRCA(net1, clades)){
//            System.out.println(x.getName());
//        }
        Set<NetNode> res = getAbandonList(net1.findNode("I2"));
        System.out.println("abandon list:");
        for (NetNode x : res){
            System.out.println(x.getName());
        }

    }

    public static void testNet(){
//        Network net = Networks.readNetwork("((((t50:0.7284880026)#H1:0.0::0.3735040297,(t39:0.1397969294,t75:0.1397969294):0.5886910732):20.33144073,((((((t80:1.348108508)#H2:0.3479510206::0.5465845262,((t83:0.1083187504,t8:0.1083187504):0.6998155581,(t23:0.5923585598,t31:0.5923585598):0.2157757488):0.8879252204),(t29,t72)),t85),t20),((t15,t16),t25))):0.702792199,(#H2:5.752299178::0.4534154738,((#H1:0.0::0.6264959703,t87:0.7284880026),Z)):14.662313241);");
//        net.resetRoot("Z");
//        System.out.println(net.toString());
//        int leafcnt = 0;
//        Set<String> leafset = new HashSet<>();
//        for(Object o: net.bfs()){
//            NetNode node = (NetNode) o;
//            if (node.isLeaf()){
//
//                String leaf = node.getName();
//                if(leafset.contains(leaf)){
//                    System.out.println("duplicate leaf:"+leaf);
//                }
//                System.out.println(leaf);
//                leafset.add(leaf);
//                leafcnt ++;
//            }
//
//        }
//        System.out.println(leafcnt);
//        System.out.println(leafset.size());
        Network net = Networks.readNetwork("((((t87:0.7284880026,(t50:0.7284880026)I20#H1:0.0::0.6264959703)I0,((t80:1.348108508)#H3::0.547)#H2::0.5),((((t25,(t16,t15)I3)I5,((t19,t20)I4,(t85,(((t28,(t10,t78)I7)I12,(t72,t29)I6)I13,(#H2::0.5,((t31:0.5923585598,t23:0.5923585598)I9:0.2157757488,(t8:0.1083187504,t83:0.1083187504)I8:0.6998155581)I10:0.8879252204)I11)I14)I15)I16)I17,((t75:0.1397969294,t39:0.1397969294)I1:0.5886910732,I20#H1:0.0::0.3735040297)I23:20.33144073)I22,#H3:5.752299178::0.4534154738)I19)I18,Z)I2;");
        System.out.println(net);

    }

    public static void main(String[] args) {
//        test();
//        testNet();
//        Network net1 = Networks.readNetwork("((((t21:0.689109324)I9#H1:1.0::0.5,(t22:0.5326000818,t10:0.5326000818)I8:0.1565092422)I5:4.114409669,((t33:1.432442835,(t16:0.4712545803,(t25:0.1626650406,t30:0.1626650406)I11:0.3085895397)I10:0.9611882551)I6:2.7344783452)I3#H2:1.0::0.5)I2:1.9346975199099998,((t11:2.87898364,(t38:0.689109324,I9#H1:1.0::0.5)I7:2.189874316)I4:1.905507974,I3#H2:1.0::0.5)I1:1.9537249)I0;(((t4:3.392681616,(((((t6:0.002301929895,t20:0.002301929895)I19:0.006815844532,t32:0.009117774426)I18:0.3733282863,t29:0.3824460607)I16:0.5613405463,(t13:0.4412615513,t24:0.4412615513)I15:0.5025250558)I14:1.897741706,(t9:0.1334509069,t1:0.1334509069)I13:2.708077406)I12:0.5511533024)I11:6.347379667,(((t38,t11)I5,(((((t30,t25)I2,t16)I3,t33)I4,(t21,(t10,t22)I1)I17)I6,t2)I7)I8,t18:8.056923686000001)I10:1.683137597)I9,Z:100.0)I0;\n");
//        Network net2 = Networks.readNetwork("(Z:100.0,((t18:8.056923686000001,((t2,(((t22,t10)I1,(t21)#H1:1.0::0.5)I17,(t33,(t16,(t25,t30)I2)I3)I4)I6)I7,(t11,(#H1:1.0::0.5,t38))I5)I8)I10:1.683137597,(((t1:0.1334509069,t9:0.1334509069)I13:2.708077406,((t24:0.4412615513,t13:0.4412615513)I15:0.5025250558,(t29:0.3824460607,(t32:0.009117774426,(t20:0.002301929895,t6:0.002301929895)I19:0.006815844532)I18:0.3733282863)I16:0.5613405463)I14:1.897741706)I12:0.5511533024,t4:3.392681616)I11:6.347379667)I9)I0;");
//        List<Network> net = matchNetworks(net1, net2);
//        System.out.println("merged");
//        for (Network n: net){
//            System.out.println(n.toString());
//        }

//        Network net1 =  Networks.readNetwork("((((((t4:0.04115588418,t10:0.04115588418):0.1650568111,t29:0.2062126953):0.4854316696,t19:0.6916443648):2.461066606,((t31:0.1486444659,t39:0.1486444659):0.4802126509,t26:0.6288571168):2.523853855):4.532398833,((t20:0.07990613658,t3:0.07990613658):6.728437711,(((((t16,t28))#H1:1.0::0.5,(t34,t11)),((((((t18,#H1:1.0::0.5),(t16,t28)),t27),((t33)#H2:1.0::0.5,t22)),#H2:1.0::0.5),t38)),(t17:0.5189419599,t2:0.5189419599):5.987387876):0.302014012):0.8767659562),Z:100.0);");
//        Network net2 = Networks.readNetwork("(Z:100.0,(((t19:0.6916443648,((t4:0.04115588418,t10:0.04115588418):0.1650568111,t29:0.2062126953):0.4854316696):2.461066606,(t26:0.6288571168,(t39:0.1486444659,t31:0.1486444659):0.4802126509):2.523853855):4.532398833,((t20:0.07990613658,t3:0.07990613658):6.728437711,((t17:0.5189419599,t2:0.5189419599):5.987387876,((t38:4.652012987,((t33:1.479894037)#H1:2.473201678::0.5535111291,((#H1:1.45797114::0.4464888709,t22:2.937865176):0.4574710767,(t27:3.188080979,(((t28:0.7278426847,t16:0.7278426847):0.2254566606)#H2:1.481069149::0.4658388651,t18:2.434368494):0.7537124848):0.2072552739):0.5577594616):0.698917272):0.7315794042,(#H2:0.1945802142::0.5341611349,(t34:0.42338387,t11:0.42338387):0.7244956895):4.235712831):1.122737445):0.302014012):0.8767659562):92.3148901963);");

        Network net1 = Networks.readNetwork("(Z:100.0,((t1:9.797942233,((t27:4.278925328,((t36:0.05663118624,t42:0.05663118624):0.6518610285,(t9:0.2405302977)#H1:0.4679619171::0.5394632656):3.570433114):4.016134627,(((((t39:0.3213490742,(t7:0.06406507032,t45:0.06406507032):0.2572840039):1.586194693,((t23:0.2621080793,t19:0.2621080793):0.1648868735)#H2:1.480548815::0.6742175817):0.8174030769,t3:2.724946844):1.041229637,((t41:0.356303595)#H3:2.941455351::0.3356585442,t29:3.297758946):0.4684175362):3.207411187,((t15:4.105503166,(#H3:0.2318478403::0.6643414558,t44:0.5881514353):3.51735173):1.341313197,(#H2:2.712915248::0.3257824183,t2:3.1399102):2.306906162):1.526771306):1.321472286):1.502882278):5.138230099,((t43:0.001440716722,t6:0.001440716722):9.375552728,(#H1:6.404664292::0.4605367344,(t8:0.5759259639,t32:0.5759259639):6.069268626):2.731798855):5.559178887):85.0638276681);");
        Network net2 = Networks.readNetwork("((((t32,t8),(t6:0.001440716722,t43:0.001440716722)),((((t9,(t42,t36)),t27),((t2,((t41,t44),t15)),(t29,(((t23,t19),((t7,t45),t39)),t3))):1.321472286):1.502882278,t1:9.797942233):5.138230099),Z:100.0);");

        double distance = Networks.computeDistanceBetweenTwoNetworks(net1, net2);
        Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(net1, net2);
        System.out.println(distance);
        System.out.println(closest.Item3);
        System.out.println(closest.Item1.toString());
        System.out.println(closest.Item2.toString());
    }

}
