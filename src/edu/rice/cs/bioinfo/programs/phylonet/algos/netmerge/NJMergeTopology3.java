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
import org.apache.commons.collections15.list.SynchronizedList;

import static edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.util.Utils.*;
import static edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.util.Utils.getAbandonList;
import static edu.rice.cs.bioinfo.programs.phylonet.algos.summarize.majorTree.removeReti;


public class NJMergeTopology3 {
    private List<Network> _subnetworks = new ArrayList<>();
    private List<Network> _originalSubNetworks = new ArrayList<>();
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

    public Set<NetNode> _abandonList = new HashSet<>();

    private boolean _addedOutgroup = false;


    /* Constructor */
    public NJMergeTopology3(List<Network> subnetworks, double[][] matrix, List<String> taxonList, String outgroup) {
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

//        checkRoot();

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
        if (_debug) {
            System.out.println("joinNodesInBothNets");
        }
        if (nodeAinT1.getName().equals("I17") && nodeBinT2.getName().equals("Z")){
            System.out.println("debug point");

        }
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

                res.add(newnode);

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
                res.add(newNodeToAddinT2);
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
                res.add(newnode);
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

            // node under nodeB to node in network


            // intersection of A and X
            Set<String> intersection = new HashSet<>(cladeX);
            intersection.retainAll(cladeA);
            Set<String> unionIntersectionB = new HashSet<>(cladeB);
            unionIntersectionB.addAll(intersection);
            if (Utils._debug){
                System.out.println("checkOverlap overlap");
                System.out.println("testjoin intersectionB:"+intersection);
                System.out.println("testjoin unionIntersectionB:"+unionIntersectionB);
            }



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
//                NetNode nodeinNet2 = redundantNodeMap.get(nodeunderA);
//                parentAX =  (NetNode) nodeinNet2.getParents().iterator().next();
//                parentAX.removeChild(nodeinNet2);

//                parentAX =  (NetNode) nodeunderA.getParents().iterator().next();
//                parentAX.removeChild(nodeunderA);

                if(isNetBackbone){
                    parentAX =  (NetNode) nodeunderA.getParents().iterator().next();
                    parentAX.removeChild(nodeunderA);
                }
                else{


                    NetNode nodeinNet2 = redundantNodeMap.get(nodeunderA);
                    parentAX =  (NetNode) nodeinNet2.getParents().iterator().next();
                    parentAX.removeChild(nodeinNet2);

                }

            }

            if(redundantNodeMap.size() == 1){
                NetNode<NetNodeInfo> parentB = (NetNode) nodeBinT2.getParents().iterator().next();
                if (parentB.getName().isEmpty()){
                    for(Object o: parentB.getChildren()){
                        NetNode child = (NetNode) o;
                        if (!child.equals(nodeBinT2)){
                            parentB.removeChild(child);
                            parentB.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);
                            parentAX.adoptChild(child, NetNode.NO_DISTANCE);
                            parentB.setName(newNodeToAddinT2.getName());
                            return parentB;
                        }
                    }
                }
            }
            else{
                //TODO:test
                if (Utils._debug){
                    System.out.println("redundantNodeMap.size() != 1,"+redundantNodeMap.size());
                    for (NetNode nodeunderA : redundantNodeMap.keySet()) {
                        System.out.println("redundantNodeMap:"+nodeunderA.getName());
                    }
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
//                    joinNodesInOneNet(nodeBinT1, cladeBinT1, _subnet, nodeAinT2, cladeAinT2, false);
//                    edited = true;
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
//                    joinNodesInOneNet(nodeBinT2, cladeBinT2, _backbone, nodeAinT1, cladeAinT1, true);
                    edited = true;
                    NetNode nA = n1.findNode(nodeA.getName());
                    NetNode nB = n2.findNode(nodeB.getName());
                    Tuple<NetNode, Network> tuple = joinNodesInOneNet(n2, n1, nB, cladeBinT2, nA, cladeAinT1, true);
                    _backbone = tuple.Item2;

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
//                    joinNodesInOneNet(nodeA, cladeAinT1, _subnet, nodeBinT2, cladeBinT2, false);
//                    edited = true;
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
                    NetNode nA = n1.findNode(nodeA.getName());
                    NetNode nB = n2.findNode(nodeB.getName());
                    Tuple<NetNode, Network> tuple = joinNodesInOneNet(n2, n1, nB, cladeBinT2, nA, cladeAinT1, true);
                    _backbone = tuple.Item2;
//                    List<NetNode> nodelist = tuple.Item1;
//                    List<Network> networklist = tuple.Item2;
//                    _abandonList.addAll(getAbandonList(nodelist.get(0)));
//                    _abandonList.addAll(getAbandonList(nodelist.get(1)));
//                    List<Network> networklist = joinNodesInBothNets(n1, nA, cladeA, n2, nB, cladeB);
//                    n1 = networklist.get(0);
//                    n2 = networklist.get(1);

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
                    NetNode nA = n2.findNode(nodeA.getName());
                    NetNode nB = n1.findNode(nodeB.getName());
                    Tuple<NetNode, Network> tuple = joinNodesInOneNet(n2, n1, nA, cladeAinT2, nB, cladeBinT1, true);
                    _backbone = tuple.Item2;
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
//                    Tuple<List<NetNode>, List<Network>> tuple = joinNodesInBothNets(_backbone, nodeBinT1, cladeBinT1, _subnet, nodeAinT2, cladeAinT2);
                    NetNode nA = n2.findNode(nodeA.getName());
                    NetNode nB = n1.findNode(nodeB.getName());
                    Tuple<NetNode, Network> tuple =  joinNodesInOneNet(n2, n1, nA, cladeAinT2, nB, cladeBinT1, true);
                    _backbone = tuple.Item2;
//                    List<NetNode> nodelist = tuple.Item1;

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

    public Tuple<NetNode, Network> joinNodesInOneNet(Network net1, Network net2,  NetNode nodeAinT1,Set<String> cladeA, NetNode nodeBinT2, Set<String> cladeB, boolean isnet2Backbone) {
        /*
        net2: backbone
        Join clade A and clade B in just one tree

        1. tree 1 as (A,X); and extract (A,X);
        2. tree 2 as (B,Y); and extract (B,Y);
        3. attach A to tree 2 as ((A,B),Y);
        */
        if (_debug) {
            System.out.println("joinNodesInOneNet");
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
            return new Tuple<>(newNodeToAdd, net2);
        }
        else {

            if (_debug) {
                System.out.println("check point before checkOverlap:"+nodeAinT1.getName()+","+nodeBinT2.getName());
                System.out.println("cladeX: " + cladeX);
            }
            NetNode newNodeToAddinT2 = new BniNetNode();
            newNodeToAddinT2.setName("I" + _nodeCnt);
            NetNode newnode  = checkOverlap(cladeX, cladeA, cladeB, net2, nodeBinT2, net1, nodeAcopy, newNodeToAddinT2, isnet2Backbone);

            return new Tuple<>(newnode, net2);

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
        return: violates score: 0.0: compatible; 100.0: incompatible
    */
    public double testJoin(NetNode nodeA, NetNode nodeB) throws Exception {
        Set<String> cladeA = Utils.getLeafSet(nodeA);
        Set<String> cladeB = Utils.getLeafSet(nodeB);
        Set<String> cladeAB = new HashSet<>(cladeA);
        cladeAB.addAll(cladeB);
        if(nodeA.getName().equals("I16")&&nodeB.getName().equals("I27")){
            System.out.println("check point");
        }


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
                    Tuple<NetNode, Network> tuple = joinNodesInOneNet(n2, n1, nB, cladeBinT2, nA, cladeAinT1, true);

                    //                    List<NetNode> netnodelist = joinNodesInBothNets(n1, nA, cladeA, n2, nB, cladeB);
//                    Tuple<List<NetNode>, List<Network>> tuple = joinNodesInBothNets(n1, nA, cladeAinT1, n2, nB, cladeBinT2);
//                    List<NetNode> netnodelist = tuple.Item1;
//                    List<Network> networklist = tuple.Item2;
//                    n1 = networklist.get(0);
//                    n2 = networklist.get(1);
//                    NetNode node1 = netnodelist.get(0);
//                    NetNode node2= netnodelist.get(1);
                    if (tuple.Item1 == null){
                        violates = _compatible_max_value;
                    }
                    else{
                        n1 = tuple.Item2;
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
                    if (node ==  null){
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
                    Tuple<NetNode, Network> tuple = joinNodesInOneNet(n2, n1, nA, cladeAinT2, nB, cladeBinT1, true);
//                    Tuple<List<NetNode>, List<Network>> tuple = joinNodesInBothNets(n1, nB, cladeBinT1, n2, nA, cladeAinT2);
//                    List<NetNode> netnodelist = tuple.Item1;
//                    List<Network> networklist = tuple.Item2;
//                    n1 = networklist.get(0);
//                    n2 = networklist.get(1);
//                    addInheritanceProb(n1);
//                    addInheritanceProb(n2);
//                    List<NetNode> netnodelist = joinNodesInBothNets(n1, nB, cladeB, n2, nA, cladeA);
//                    NetNode node1 = netnodelist.get(0);
//                    NetNode node2 = netnodelist.get(1);
                    if(tuple.Item1 == null ){
                        violates = _compatible_max_value;
                    }
                    else {
                        n1 = tuple.Item2;
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
//        if(violates > 0.01) {
//            System.out.println("violates");
//        }
        return violates;
    }



    public List<Network> mergeNetsViaNJ() {

        // Check trees are on disjoint leaf sets
        for (int i = 0; i < _leaves.size() - 1; i++) {
            for (int j = i + 1; j < _leaves.size(); j++) {
                Set<String> intersection = new HashSet<String>(_leaves.get(i));
                intersection.retainAll(_leaves.get(j));
                if (!intersection.isEmpty()) {
                    if (_debug){
                        System.out.println(_subnetworks.get(i).toString());
                        System.out.println(_subnetworks.get(j).toString());
                    }

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
            if (n == 8){
                System.out.println("debug");
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
                System.out.println("Unable to find valid siblinghood!");
                return null;

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
                _subnetworks.set(0, _backbone);
                _leaves.set(0, Utils.getLeafSet(_backbone.getRoot()));

                if (_leaves.get(0).equals(new HashSet<>(_taxonList))) {
                    //todo something
                    // System.out.println("Able to quit early!");
                    mergedNetworkID = 0;
                    if (_debug)
                    {
                        System.out.println("equal leaf set:"+_backbone.toString());
                    }

                }
            }


            if (mergedNetworkID >= 0){
                break;
            }
            // Create the new node

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

        removeInternalLeaf(_subnetworks);
        for (Network net: _subnetworks) {
            NetNode outgroupNode = net.findNode(_outgroup);
            if (outgroupNode != null) {
                net.resetRoot(_outgroup);
                if (_addedOutgroup) {
                    net.getRoot().removeChild(outgroupNode);
                    Networks.removeBinaryNodes(net);
                }
            }
        }

        Queue<Tuple<Network, Integer>> netQueue = new LinkedList<>();
        List<Network> candidateNetworks = new ArrayList<>();
        Network candidateNetwork = _subnetworks.get(mergedNetworkID);
        candidateNetworks.add(candidateNetwork);
        List<Network> backboneList2 = new ArrayList<>();
        backboneList2.addAll(Pipeline.getAllBackboneNets(candidateNetwork, Integer.MAX_VALUE));
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
                if(_debug){
                    System.out.println("match networks:");

                    System.out.println(originalSubnet.toString());
                    System.out.println(net.toString());
                }

                List<Network> matchedlist = matchNetworks(originalSubnet, net);
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
        if(_debug){
            for (Network resnet: resultnetwork){
                System.out.println(resnet.toString());

            }
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

    public static void main(String[] args) {
        Network net1 = Networks.readNetwork("(((((t31,(((t23:0.05570223274,(t54)#H1:1.0::0.5),t37),t55)),t75),((t73,t76),((t20,t79),(t33,(t71,(#H1:1.0::0.5,t14)))))),((((t5,t88),t22),(t10,t68)),((t44,(t24,t86)),((((t29,t83),(t39,((t67,(t34,(t47,t52))),(t3,t50)))),(t78,t8)),(((t70,(t49,t6)),((t56,t42),t66)),(t69,t11)))))),Z);");
        Network net2 = Networks.readNetwork("(Z:100.0,(((((t13:0.7862574784,t20:0.7862574784):0.2110213342,(t101:0.2697379102)#H1:0.7275409024::0.617067135):3.191420771,((((t95:0.7678981129,t25:0.7678981129):0.4870647091)#H2:1.0::0.4485395048,(t17:0.5847538844,t107:0.5847538844):0.6702089375):0.173015734,(t119:0.4655055949,(t51:0.2645597783,(t92:0.07271986976,t30:0.07271986976):0.1918399086):0.2009458166):0.9624729611):2.760721027):4.241777967,((t5:1.16584051)#H3:2.91462107::0.3677499908,(t38:0.7504702136,(t83:0.2691405969,(t104:0.1865110569,t63:0.1865110569):0.08262954):0.4813296166):3.329991367):4.35001597):9.870800433,(t113:8.036259856,((((t9:0.634638641,(t47:0.1726339417,t80:0.1726339417):0.4620046993):0.7202467625,(t7:0.4481092484)#H4:0.906776155::0.6218269879):4.291666011,((#H1:0.07959503738::0.382932865,t111:0.3493329476):2.652980824,#H3:1.836473261::0.6322500092):2.644237643):0.6120568492,(((t94:2.709348407,((t58:0.8959147328)#H5:0.9631728146::0.6644778975,t11:1.859087547):0.8502608598):1.955251247,((#H4:2.439989918::0.3781730121,t84:2.888099166):0.7407316471,((t78:0.4867129676)#H6:1.528097656::0.4894336475,t10:2.014810624):1.61402019):1.035768841):0.108053052,(((#H5:0.01075035333::0.3355221025,t1:0.9066650862):2.675132952,(((t90:0.7597486022,t56:0.7597486022):1.754998093,(#H2:0.5234249969::0.5514604952,((t36:0.08657109689)#H7:0.6749242745::0.5050840506,((#H7:1.0::0.4949159494,(t45:0.05570802959,t112:0.05570802959):0.0308630673):0.5575321296,(t42:0.3380560209,t75:0.3380560209):0.3060472055):0.1173921449):1.016892447):0.7363588763):0.7511832825,(t6:0.0862916006,t52:0.0862916006):3.179638377):0.315868061):0.5683914573,((t105:0.4421176959,t49:0.4421176959):0.04459527164,#H6:1.0::0.5105663525):3.663476528):0.6224632102):1.485955557):1.777651592):10.26501813):81.69872201519999);");
        System.out.println(Pipeline.CheckWithTrueNetwork(net1, net2).Item3);
    }

}
