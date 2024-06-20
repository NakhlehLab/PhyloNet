package edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge;
/*
 * @ClassName:   Merge
 * @Description:
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
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.Pipeline;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.SuperNetwork3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;


public class NJMerge {
    private List<Network> _subnetworks = new ArrayList<>();
    private int _num_taxa;
    private double[][] _matrix = new double[_num_taxa][_num_taxa];
    private List<String> _taxonList = new ArrayList<>();
    private List<Set<String>> _leaves = null;
    private String _outgroup = "";
    public static double _epsilon = 0.001;
    private int _nodeCnt = 0;
    private double _scale = 1;

    /* Constructor */
    public NJMerge(List<Network> subnetworks, double[][] matrix, List<String> taxonList) {
        this._subnetworks = subnetworks;
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

    public static NetNode<NetNodeInfo> copyNodeSubnet(NetNode<NetNodeInfo> node, Map<NetNode, NetNode> old2new){

        NetNode<NetNodeInfo> root = new BniNetNode();
        if (old2new.keySet().contains(node)){
            return old2new.get(node);
        }
        if (node.isLeaf()){
            root.setName(node.getName());
            root.setData(new NetNodeInfo(0.0));
        }
        else{
            root.setName(node.getName());
            root.setData(new NetNodeInfo( node.getData().getHeight()));
            for (Object o: node.getChildren()){
                NetNode<NetNodeInfo> child = copyNodeSubnet((NetNode) o, old2new);
                double distance = ((NetNode<NetNodeInfo>) o).getParentDistance(node);
                root.adoptChild(child, distance);
//                root.setData(new NetNodeInfo(distance+child.getData().getHeight()));
            }
        }
        old2new.put(node, root);
        return root;
    }

    /*Join clade A and clade B in both nets
    */
    public List<Network> joinNodesInBothNets(Network net1, NetNode nodeAinT1, List<String> cladeAList,
                                                  Network net2, NetNode nodeBinT2, List<String> cladeBList,
                                                  boolean test) {

        Set<String> cladeA = new HashSet<>(cladeAList);
        Set<String> cladeB = new HashSet<>(cladeBList);
        Set<String> leaves1 = getLeafSet(net1.getRoot());
        Set<String> leaves2 = getLeafSet(net2.getRoot());

        boolean cladeAisT1 = leaves1.equals(cladeA);
        boolean cladeBisT2 = leaves2.equals(cladeB);

        if (cladeAisT1 && cladeBisT2) {
            // Join trees 1 and 2 together, i.e.,
            // tree 1 will equal tree 2
            if (test) {
                return Arrays.asList(null, null);
            }
            //todo distances
            NetNode root = new BniNetNode();
            root.adoptChild(nodeAinT1, NetNode.NO_DISTANCE);
            root.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
            net1 = new BniNetwork((BniNetNode) root);
            net2 = null;
        } else if (cladeAisT1) {
            if (test) {
                return Arrays.asList(null, null);
            }
            NetNode root = new BniNetNode();
            root.adoptChild(nodeAinT1, NetNode.NO_DISTANCE);
            root.adoptChild(net2.getRoot(), NetNode.NO_DISTANCE);
            net1 = new BniNetwork((BniNetNode) root);
            net2 = null;
        } else if (cladeBisT2) {
            if (test) {
                return Arrays.asList(null, null);
            }
            NetNode root = new BniNetNode();
            root.adoptChild(net1.getRoot(), NetNode.NO_DISTANCE);
            root.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
            net1 = new BniNetwork((BniNetNode) root);
            net2 = null;
        } else {
            Map<NetNode, NetNode> old2newT1 = new HashMap<>();
            Map<NetNode, NetNode> old2newT2 = new HashMap<>();
            NetNode nodeAinT1copy = copyNodeSubnet(nodeAinT1, old2newT1);
            NetNode nodeBinT2copy = copyNodeSubnet(nodeBinT2, old2newT2);

            NetNode parentA = (NetNode) nodeAinT1.getParents().iterator().next();
            NetNode newParentNodeT1 = new BniNetNode();
            parentA.removeChild(nodeAinT1);
            parentA.adoptChild(newParentNodeT1, NetNode.NO_DISTANCE);
            newParentNodeT1.adoptChild(nodeAinT1, NetNode.NO_DISTANCE);
            newParentNodeT1.adoptChild(nodeBinT2copy, NetNode.NO_DISTANCE);

            NetNode parentB = (NetNode) nodeBinT2.getParents().iterator().next();
            NetNode newParentNodeT2 = new BniNetNode();
            parentB.removeChild(nodeBinT2);
            parentB.adoptChild(newParentNodeT2, NetNode.NO_DISTANCE);
            newParentNodeT2.adoptChild(nodeAinT1copy, NetNode.NO_DISTANCE);
            newParentNodeT2.adoptChild(nodeBinT2, NetNode.NO_DISTANCE);
        }

        return Arrays.asList(net1, net2);
    }

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
        Tuple<Network, Map<NetNode, NetNode>> tuple = cloneNetwork(net1);
        Network curnet = tuple.Item1;

        for(Object nodeObj : curnet.dfs()) {
            NetNode node = (NetNode) nodeObj;
            if(node.isNetworkNode()) {
                curReticulations.add(node);
            }
        }
        System.out.println(curnet.toString());
        List<Tuple<NetNode, NetNode>> removedRetiList = new ArrayList<>();

        getAllBackbonesDfs(0, curnet, curReticulations, nodeCur2Net2, closest.Item2, removedRetiList);

        //add reticulation edges to net2
        //todo debug not the same node
        Map<NetNode, Set<String>> node2leaf1 = new HashMap<>();
        Map<NetNode, Set<String>> node2leaf2 = new HashMap<>();
        getNodeLeafMap(net1.getRoot(), node2leaf1, null);
        getNodeLeafMap(net2.getRoot(), node2leaf2, null);
        System.out.println("reticulations");
        for (Tuple<NetNode, NetNode> edge: removedRetiList){
            NetNode head = nodeCur2Net2.get(edge.Item1);
            NetNode tail = nodeCur2Net2.get(edge.Item2);
            System.out.println(head+","+ tail);
        }

        return nodeCur2Net2;
    }

//    public static Map<Set<String>, List<NetNode>> getMapReverse(Map<NetNode, Set<String>> map){
//
//
//    }

//    private static <V, K> Map<V, List<K>> invertMapUsingGroupingBy(Map<K, V> map) {
//        Map<V, List<K>> inversedMap = map.entrySet()
//                .stream()
//                .collect(Collectors.groupingBy(Map.Entry::getValue, Collectors.mapping(Map.Entry::getKey, Collectors.toList())));
//        return inversedMap;
//    }

    private static <V, K> Map<V, List<K>> invertMapUsingGroupingBy(Map<K, V> map) {
        Map<V, List<K>> inversedMap = new HashMap<>();
        for (K key: map.keySet()){
            boolean found = false;
            for (V value: inversedMap.keySet()){
                if (value.equals(map.get(key))){
                    inversedMap.get(value).add(key);
                    found = true;
                }
            }
            if (! found){
                V value = map.get(key);
                inversedMap.put(value, new ArrayList<>());
                inversedMap.get(value).add(key);
            }

        }
        return inversedMap;
    }

    // Add newVnode to edge(nodeVParent, nodeV), newRetic to edge(nodeU, child), and edge(newVnode, newRetic)
    public static void addReticulation(NetNode<NetNodeInfo> nodeV, NetNode<NetNodeInfo> nodeVParent, NetNode<NetNodeInfo> newVnode, NetNode<NetNodeInfo> child, NetNode<NetNodeInfo> nodeU, NetNode<NetNodeInfo> newRetic, double newVHeight, double retiHeight){
//        NetNode<NetNodeInfo> newVnode = new BniNetNode();
//        NetNode<NetNodeInfo> newRetic = new BniNetNode();
        NetNodeInfo niRetic = new NetNodeInfo(retiHeight);
        newRetic.setData(niRetic);
        NetNodeInfo niV = new NetNodeInfo(newVHeight);
        newVnode.setData(niV);
        newVnode.adoptChild(newRetic, newVnode.getData().getHeight() - newRetic.getData().getHeight());
        newVnode.adoptChild(nodeV, newVnode.getData().getHeight() - nodeV.getData().getHeight());
        nodeVParent.adoptChild(newVnode, nodeVParent.getData().getHeight()-newVnode.getData().getHeight());
        nodeVParent.removeChild(nodeV);
        nodeU.adoptChild(newRetic, nodeU.getData().getHeight() - newRetic.getData().getHeight());
        newRetic.adoptChild(child, newRetic.getData().getHeight() - child.getData().getHeight());
        nodeU.removeChild(child);
    }

    public static void removeReticulation(NetNode<NetNodeInfo> nodeV, NetNode<NetNodeInfo> nodeVParent, NetNode<NetNodeInfo> newVnode, NetNode<NetNodeInfo> child, NetNode<NetNodeInfo> nodeU, NetNode<NetNodeInfo> newRetic){
        newVnode.removeChild(nodeV);
        newVnode.removeChild(newRetic);
        nodeVParent.removeChild(newVnode);
        nodeVParent.adoptChild(nodeV, nodeVParent.getData().getHeight() - nodeV.getData().getHeight());
        nodeU.removeChild(newRetic);
        nodeU.adoptChild(child, nodeU.getData().getHeight() - child.getData().getHeight());
        newRetic.removeChild(child);

    }

    //TODO: finish this function
    //match the reticulations in two networks, add reticulations on net1 to net2
    public static List<Network> matchNetworks(Network<NetNodeInfo> net1, Network<NetNodeInfo> net2){
//        Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(net1, net2);
        Network net1copy = Networks.readNetwork(net1.toString());
        Network net2copy = Networks.readNetwork(net2.toString());
        Networks.autoLabelNodes(net1copy);
        Networks.autoLabelNodes(net2copy);
        System.out.println(net1copy.toString());
        System.out.println(net2copy.toString());
        Map<NetNode, Set<String>> node2leaf1 = new HashMap<>();
        Map<NetNode, Set<String>> node2leaf2 = new HashMap<>();

        List<Network> candidateNetworks = new ArrayList<>();

        Set<String> intersection = getLeafSet(net1.getRoot());
        intersection.retainAll(getLeafSet(net2.getRoot()));

        getNodeLeafMap(net1.getRoot(), node2leaf1, intersection);
        getNodeLeafMap(net2.getRoot(), node2leaf2, intersection);


        Map<Set<String>, List<NetNode>> set2node2 = invertMapUsingGroupingBy(node2leaf2);

        for (NetNode<NetNodeInfo> retic: net1.getNetworkNodes()){
            double reticHeight = retic.getData().getHeight();
            Set<String> reticLeafSet = node2leaf1.get(retic);
//            List<Set<String>> parentLeafList = new ArrayList<>();
//            Set<String> ulist = new HashSet<>();
            Iterator<NetNode<NetNodeInfo>> parentIter = retic.getParents().iterator();
            NetNode<NetNodeInfo> nodeU = parentIter.next();
            NetNode<NetNodeInfo> nodeV = parentIter.next();
            Set<String> nodeUleaf = node2leaf1.get(nodeU);
            Set<String> nodeVleaf = node2leaf1.get(nodeV);
            NetNode<NetNodeInfo> nodeUParent = nodeU.getParents().iterator().next();
            NetNode<NetNodeInfo> nodeVParent = nodeV.getParents().iterator().next();
            Set<String> nodeVParentLeaf = node2leaf1.get(nodeVParent);
            Set<String> nodeUParentLeaf = node2leaf1.get(nodeUParent);


            if (set2node2.containsKey(reticLeafSet)){
                for (NetNode<NetNodeInfo> child: set2node2.get(reticLeafSet)){
                    for (Object op: child.getParents()){
                        NetNode<NetNodeInfo> parent = (NetNode<NetNodeInfo>) op;
                        if (node2leaf2.get(parent).equals(nodeUleaf)){
                            List<NetNode> nodeVParentList = new ArrayList<>();
                            for (Set<String> temp: set2node2.keySet()){
                                if(temp.equals(nodeVParentLeaf)){
                                    nodeVParentList = set2node2.get(temp);
                                    break;
                                }
                            }
                            if (nodeVParentList.size() > 0){
                                for (NetNode<NetNodeInfo> nodeVParent2: nodeVParentList){
//                                    if (nodeVParent2.getData().getHeight() > retic.getData().getHeight()) {
//                                        continue;
//                                    }
                                    List<NetNode> nodeVParent2ChildList = StreamSupport.stream(nodeVParent2.getChildren().spliterator(), false)
                                            .collect(Collectors.toList());

                                    for (NetNode<NetNodeInfo> nodeV2: nodeVParent2ChildList){
                                        if (nodeV2.equals(parent)){
                                            continue;
                                        }
                                        double newVHeight = nodeV.getData().getHeight();
                                        NetNode<NetNodeInfo> newVnode = new BniNetNode();
                                        NetNode<NetNodeInfo> newRetic = new BniNetNode();
                                        addReticulation(nodeV2, nodeVParent2, newVnode, child, parent, newRetic, newVHeight, reticHeight);
                                        Network newnet = net2.clone();
                                        candidateNetworks.add(newnet);
                                        removeReticulation(nodeV2, nodeVParent2, newVnode, child, parent, newRetic);
                                    }
                                }
                            }
                        }

                        else if (node2leaf2.get(parent).equals(nodeVleaf)){
                            List<NetNode> nodeUParentlist = new ArrayList<>();
                            for (Set<String> temp: set2node2.keySet()) {
                                if (temp.equals(nodeUParentLeaf)) {
                                    nodeUParentlist = set2node2.get(temp);
                                    break;
                                }
                            }

                            if (nodeUParentlist.size() > 0){
                                for (NetNode<NetNodeInfo> nodeUParent2: nodeUParentlist){
//                                    if (nodeU2.getData().getHeight() > retic.getData().getHeight()) {
//                                        continue;
//                                    }
                                    List<NetNode> nodeUParent2ChildList = StreamSupport.stream(nodeUParent2.getChildren().spliterator(), false)
                                            .collect(Collectors.toList());

                                    for (NetNode<NetNodeInfo> nodeU2: nodeUParent2ChildList) {
                                        if (nodeU2.equals(parent)){
                                            continue;
                                        }
//                                        if (nodeUparent.getData().getHeight() < retic.getData().getHeight()) {
//                                            continue;
//                                        }
                                        double newUHeight = nodeU.getData().getHeight();
                                        NetNode<NetNodeInfo> newVnode = new BniNetNode();
                                        NetNode<NetNodeInfo> newRetic = new BniNetNode();
                                        addReticulation(nodeU2, nodeUParent2, newVnode, child, parent, newRetic, newUHeight, reticHeight);
                                        Network newnet = net2.clone();
                                        candidateNetworks.add(newnet);
                                        removeReticulation(nodeU2, nodeUParent2, newVnode, child, parent, newRetic);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        System.out.println("candidate networks after matching two networks are below");
        for (Network net: candidateNetworks){
            System.out.println(net.toString());
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

            Tuple<Network, Map<NetNode, NetNode>> tuple = cloneNetwork(cur);
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

    public static Tuple<Network, Map<NetNode, NetNode>> cloneNetwork(Network network){
        Network currentNetwork = network;

        // bottom-up build subnetwork
        Map<NetNode, BniNetNode> old2new = new HashMap<>();
        Map<NetNode, Integer> hitCount = new HashMap<>();
        for(Object leafnode : network.getLeaves()) {
            String leaf = ((NetNode) (leafnode)).getName();
            Set<NetNode> visited = new HashSet<>();
            Queue<NetNode> queue = new LinkedList<>();
            NetNode node = currentNetwork.findNode(leaf);
            BniNetNode newleaf = new BniNetNode();
            newleaf.setName(node.getName());
            old2new.put(node, newleaf);
            queue.add(node);
            while(queue.size() > 0) {
                node = queue.poll();
                if(visited.contains(node)) continue;
                if(!hitCount.containsKey(node)) {
                    hitCount.put(node, 0);
                }
                hitCount.put(node, hitCount.get(node) + 1);
                visited.add(node);
                for(Object parentObj : node.getParents()) {
                    NetNode parent = (NetNode) parentObj;
                    queue.add(parent);
                    BniNetNode newParentNode;
                    if(old2new.containsKey(parent)) {
                        newParentNode = old2new.get(parent);
                    } else {
                        newParentNode = new BniNetNode();
                        newParentNode.setName(parent.getName());
                        old2new.put(parent, newParentNode);
                    }

                    newParentNode.adoptChild(old2new.get(node), node.getParentDistance(parent));
                    old2new.get(node).setParentProbability(newParentNode, node.getParentProbability(parent));
                    old2new.get(node).setParentSupport(newParentNode, node.getParentSupport(parent));
                }
            }
        }

        Network newsubnetwork = new BniNetwork(old2new.get(currentNetwork.getRoot()));
        newsubnetwork.getRoot().setRootPopSize(currentNetwork.getRoot().getRootPopSize());

        Map<NetNode, NetNode> new2old = new HashMap<>();
        for(NetNode oldNode : old2new.keySet()) {
            new2old.put(old2new.get(oldNode), oldNode);
        }

        return new Tuple<>(newsubnetwork, new2old);
    }



    // pair subnetworks together
    public Tuple<Network, Map<NetNode, NetNode>> pairNets(Network net1, Network net2){
        Set<String> intersection = getLeafSet(net1.getRoot());
        intersection.retainAll(getLeafSet(net2.getRoot()));
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



    public static void getAllEdges(NetNode<NetNodeInfo> node, NetNode<NetNodeInfo> parent, Set<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edgeList){
        if (!node.isLeaf()){
            for (NetNode<NetNodeInfo>  child: node.getChildren()){
                getAllEdges(child, node, edgeList);
            }
        }
        if (parent != null){
            edgeList.add(new Tuple<>(parent, node));
        }
    }

    public static void getPath2Root(Network<NetNodeInfo> net, String leafName, Set<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> setEdgeInPath){
        findEdgeInPathRecursive(net.findNode(leafName), setEdgeInPath);
    }

    public static void findEdgeInPathRecursive(NetNode<NetNodeInfo> node, Set<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> setEdgeInPath){
        if (node.isRoot()){
            return ;
        }
        for (NetNode<NetNodeInfo> parent: node.getParents()){
            setEdgeInPath.add(new Tuple<>(parent, node));
        }
    }

    public static Set<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> findEdgeInRange(Network<NetNodeInfo> net, double height, String outgroup, boolean excludeOutgroupPath){
        Set<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edgeList = new HashSet<>();
        Set<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> candidateEdgeList = new HashSet<>();
        Set<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edgeInOutgroupPath = new HashSet<>();
        getAllEdges(net.getRoot(), null, edgeList);
        if (excludeOutgroupPath){
            getPath2Root(net, outgroup, edgeInOutgroupPath);
        }
        for (Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge: edgeList){
            if (height > edge.Item1.getData().getHeight() && height < edge.Item2.getData().getHeight() && !edgeInOutgroupPath.contains(edge)){
                candidateEdgeList.add(edge);
            }
        }
        return candidateEdgeList;

    }


    //TODO: add root of net1 to net0
//    public void joinNodesInOneNet(Network<NetNodeInfo> net0, Network<NetNodeInfo> net1){
//        NetNode<NetNodeInfo> net1Root = net1.getRoot();
//        Map<NetNode, NetNode> old2new  = new HashMap<>();
//        NetNode<NetNodeInfo> net1RootCopy = copyNodeSubnet(net1Root, old2new);
//        Set<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> candidateEdgeList = new HashSet<>();
//        double height = net1Root.getData().getHeight();
//        Set<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edgesInRange = findEdgeInRange(net0, height, _outgroup, true);
//        for (Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge: edgesInRange){
//            NetNode<NetNodeInfo> newNode = new BniNetNode<>();
////            newNode.setData(new NetNodeInfo(height));
////            edge.Item1.adoptChild(newNode, edge.Item1.getData().getHeight() - height);
////            newNode.adoptChild(edge.Item2, height - edge.Item2.getData().getHeight());
////            newNode.adoptChild(net1RootCopy, );
//
//        }
//    }

    //Join cladeA and cladeB in one or both trees
    public Tuple<List<Boolean>, NetNode<NetNodeInfo>> joinNodes(NetNode<NetNodeInfo> nodeA, NetNode<NetNodeInfo> nodeB, double sig_A, double sig_B) {
        Set<String> cladeA = getLeafSet(nodeA);
        Set<String> cladeB = getLeafSet(nodeB);

        if (!Collections.disjoint(cladeA, cladeB)) {
            throw new RuntimeException("Nodes are not disjoint on their leaf sets!");
        }

        List<Boolean> edits = Stream.generate(() -> false).limit(_subnetworks.size()).collect(Collectors.toList());

        //todo: add situation where two networks are to be merged together
//        if (cladeA.equals(_leaves.get(0)) && cladeB.equals(_leaves.get(1))) {
//            //cladeA in net0, cladeB in net1
//            if (cladeA.contains(_outgroup)){
//                //add netB to netA
//
//                if (sig_A < _epsilon){
//                    //if nodeheight of A > B
//                    // add netB inside netA
//
//                }
//                else{
//                    NetNode nodeBinT = getNodeFromClade(_subnetworks.get(1), cladeB);
//                    Map<NetNode, NetNode> old2new  = new HashMap<>();
//                    NetNode nodeBCopy = copyNodeSubnet(nodeBinT, old2new);
//                    NetNode root = new BniNetNode();
//                    root.adoptChild(_subnetworks.get(0).getRoot(), sig_A);
//                    root.adoptChild(nodeBCopy, sig_B);
//                    root.setData(new NetNodeInfo(((NetNode<NetNodeInfo>)_subnetworks.get(0).getRoot()).getData().getHeight()+sig_A));
//                    _subnetworks.set(0, new BniNetwork((BniNetNode) root));
//
//                }
//            }
//            else if (cladeB.contains(_outgroup)){
//                //add netA to netB
//                if (sig_B < _epsilon){
//                    // add netA inside netB
//
//                }
//                else{
//                    NetNode nodeBinT = getNodeFromClade(_subnetworks.get(1), cladeB);
//                    Map<NetNode, NetNode> old2new  = new HashMap<>();
//                    NetNode nodeBCopy = copyNodeSubnet(nodeBinT, old2new);
//                    NetNode root = new BniNetNode();
//                    root.adoptChild(_subnetworks.get(0).getRoot(), sig_A);
//                    root.adoptChild(nodeBCopy, sig_B);
//                    root.setData(new NetNodeInfo(((NetNode<NetNodeInfo>)_subnetworks.get(0).getRoot()).getData().getHeight()+sig_A));
//                    _subnetworks.set(0, new BniNetwork((BniNetNode) root));
//                }
//
//
//            }
//        }
//        else if (cladeB.equals(_leaves.get(0))&&cladeA.equals(_leaves.get(1))) {
//            if (cladeA.contains(_outgroup)){
//                //add netB to netA
//                if (sig_B < _epsilon){
//                    // add netB inside netA
//
//                }
//                else{
//
//                }
//            }
//            else if (cladeB.contains(_outgroup)){
//                //add netA to netB
//
//            }
//        }
        NetNode<NetNodeInfo> newNodeToAdd = new BniNetNode();
        newNodeToAdd.setName("I"+_nodeCnt);
        if ((cladeA.equals(_leaves.get(0)) && cladeB.equals(_leaves.get(1))) || (cladeA.equals(_leaves.get(1)) && cladeB.equals(_leaves.get(0)))){
            edits.set(0, true);
            edits.set(1, true);

            double distance = nodeA.getData().getNjDistances().get(nodeB.getName());
            double height1 = ((NetNode<NetNodeInfo>)_subnetworks.get(0).getRoot()).getData().getHeight();
            double height2 = ((NetNode<NetNodeInfo>)_subnetworks.get(1).getRoot()).getData().getHeight();
            double rootHeight = (distance * 2 + height1 + height2)/2;
            newNodeToAdd.adoptChild(_subnetworks.get(0).getRoot(), rootHeight - height1);
            newNodeToAdd.adoptChild(_subnetworks.get(1).getRoot(), rootHeight - height2);
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
                    newNodeToAdd = getNodeFromClade(_subnetworks.get(i), cladeUnion);
                    newNodeToAdd.setName("I"+_nodeCnt);
                    break;

                } else if (nAinT) {
                    // Add node B to T
                    edits.set(i, true);
                    Map<NetNode, NetNode> old2new = new HashMap<>();
                    NetNode<NetNodeInfo> nodeBcopy = copyNodeSubnet(nodeB, old2new);
                    nodeAinT.setName(nodeA.getName());
                    nodeBcopy.setName(nodeB.getName());
                    if (leaf.equals(cladeA)) {//todo add node height
                        newNodeToAdd.adoptChild(nodeAinT, sig_A);
                        newNodeToAdd.setData(new NetNodeInfo(nodeAinT.getData().getHeight()+sig_A));
                        newNodeToAdd.adoptChild(nodeBcopy, newNodeToAdd.getData().getHeight() - nodeBcopy.getData().getHeight());
                        _subnetworks.set(i, new BniNetwork((BniNetNode) newNodeToAdd));
                    }
                    else{
                        double height = sig_A + nodeAinT.getData().getHeight();
                        NetNodeInfo ni  = new NetNodeInfo(height);
                        newNodeToAdd.setData(ni);
                        Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> tuple = getParentAboveHeightBFS(nodeAinT, height);
                        if (tuple == null){
                            System.out.println("null here");
                        }
                        NetNode<NetNodeInfo> parentA = tuple.Item1;
                        parentA.removeChild(tuple.Item2);
                        parentA.adoptChild(newNodeToAdd, parentA.getData().getHeight() - height);
                        newNodeToAdd.adoptChild(tuple.Item2, height - tuple.Item2.getData().getHeight());
//                        newNodeToAdd.adoptChild(nodeBcopy, height - nodeBcopy.getData().getHeight());
                        newNodeToAdd.adoptChild(nodeBcopy, sig_B);
                    }

                } else if (nBinT) {
                    // Add node A to T
                    edits.set(i, true);
                    Map<NetNode, NetNode> old2new = new HashMap<>();
                    NetNode<NetNodeInfo> nodeAcopy = copyNodeSubnet(nodeA, old2new);
                    nodeAcopy.setName(nodeA.getName());
                    nodeBinT.setName(nodeB.getName());
                    if (leaf.equals(cladeB)) {
                        newNodeToAdd.adoptChild(nodeBinT, sig_B);
                        newNodeToAdd.setData(new NetNodeInfo(nodeBinT.getData().getHeight()+sig_B));
                        newNodeToAdd.adoptChild(nodeAcopy, newNodeToAdd.getData().getHeight() - nodeAcopy.getData().getHeight());
                        _subnetworks.set(i, new BniNetwork((BniNetNode) newNodeToAdd));
                    } else {//there could be several networks
//                        double height = sig_B + nodeBinT.getData().getHeight();
                        NetNodeInfo ni  = new NetNodeInfo();
                        newNodeToAdd.setData(ni);

//                        Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> tuple = getParentAboveHeightBFS(nodeBinT, height);
//                        if (tuple == null){
//                            System.out.println("null here");
//                        }
                        NetNode<NetNodeInfo> parentB = nodeBinT.getParents().iterator().next();
                        parentB.removeChild(nodeBinT);
                        parentB.adoptChild(newNodeToAdd, NetNode.NO_DISTANCE);
                        newNodeToAdd.adoptChild(nodeBinT, NetNode.NO_DISTANCE);
                        //todo: do we need to consider the case?
                        // case1: when the height is higher than the root
                        // case2: when the 'height' is exactly the same as any edge above the node

//                        newNodeToAdd.adoptChild(nodeAcopy, height - nodeAcopy.getData().getHeight());
                        newNodeToAdd.adoptChild(nodeAcopy, NetNode.NO_DISTANCE);

                    }
                }
            }
        }
        Map<NetNode, NetNode> old2new = new HashMap<>();
        NetNode newNode = copyNodeSubnet(newNodeToAdd, old2new);
        return new Tuple<>(edits, newNode);
    }

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
            if (parent.getData().getHeight() > height && node.getData().getHeight() < height){
                return new Tuple<>(parent, node);
            }
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


    public static Set<String> getLeafSet(NetNode node){
        Set<String> leafset = new HashSet<>();
        if (node.isLeaf()){
            leafset.add(node.getName());

        }
        else{
            for (Object o: node.getChildren()){
                NetNode child = (NetNode) o;
                if (child.isLeaf()){
                    leafset.add(child.getName());
                }
                else{
                    leafset.addAll(getLeafSet(child));
                }
            }
        }

        return leafset;

    }
    public NetNode getMRCA(Network net, Set<String> clade){
//        Map<NetNode, Set<String>> nodeLeafset = new HashMap<>();
        NetNode mrca = null;
        for (Object o : net.bfs()) {
            NetNode curnet = (NetNode) o;
            Set<String> leafset = getLeafSet(curnet);
            if (leafset.containsAll(clade)){
                mrca = curnet;
            }

        }
        return mrca;
    }

//    public void List<String> getTaxaNamesUnderReticulation(Network net) {
//        int count = 0;
//        List<String> results = new ArrayList<>();
//        for(Object leafObject : net.getLeaves()) {
//            NetNode leaf = (NetNode) leafObject;
//            NetNode node = leaf;
//            while(!node.isRoot()) {
//                if(node.isNetworkNode()) {
//                    results.add(leaf.getName());
//                    break;
//                }
//                node = (NetNode) node.getParents().iterator().next();
//            }
//        }
//        return results;
//    }

    public NetNode getNodeFromClade(Network net, Set<String> clade) {
        // Returns the MRCA node of the clade igoring the branch lengths
        //Todo:check whether it's better to use branch length

        Set<String> leaves = getLeafSet(net.getRoot());

        // Check if the clade is the whole tree!
        if (leaves.equals(new HashSet<>(clade))) {
            return net.getRoot();
        }

        // Check if the clade contains leaves not in the tree itself
        if (!leaves.containsAll(clade)) {
            return null;
        }

        // Encode labels as split (integer) and return the node or none
//        int split = net.getTaxonNamespace().taxaBitmask(clade);
//        if (splitToNodeMap.containsKey(split)) {
//            return splitToNodeMap.get(split);
//        } else {
//            return null;
//        }
        NetNode node = getMRCA(net, clade);
        Set<String> nodeLeaves = getLeafSet(node);
        if(clade.containsAll(nodeLeaves)){
            return node;
        }
        else{
            return null;
        }
    }


    // Todo: how compatible: identical topology?
    public boolean areNetsCompatible(Network net1, Network net2){
        Set<String> intersection = getLeafSet(net1.getRoot());
        intersection.retainAll(getLeafSet(net2.getRoot()));
//        System.out.println(net1.toString());
//        System.out.println(net2.toString());
        System.out.println(intersection);
        if (intersection.size() <= 1){
            return true;
        }
        Tuple<Network, Map<NetNode, NetNode>> subnet1 = SuperNetwork3.getSubNetwork(net1, new ArrayList<>(intersection), true);
        Tuple<Network, Map<NetNode, NetNode>> subnet2 = SuperNetwork3.getSubNetwork(net2, new ArrayList<>(intersection), true);

        Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(subnet1.Item1, subnet2.Item1);
//        if (Networks.hasTheSameTopology(subnet1.Item1, subnet2.Item1)){
//            return true;
//        }
        if (closest.Item3 < 0.01){
            return true;
        }
        else{
            return false;
        }
    }

    /*
        Test whether joining cladeA and cladeB in one
        or both trees causes the two trees to be incompatible
    */
    public boolean testJoin(NetNode nodeA, NetNode nodeB) throws Exception {
        Set<String> cladeA = getLeafSet(nodeA);
        Set<String> cladeB = getLeafSet(nodeB);
        Set<String> cladeAB = new HashSet<>(cladeA);
        cladeAB.addAll(cladeB);

        if (!Collections.disjoint(cladeA, cladeB)) {
            throw new Exception("Nodes are not disjoint on their leaf sets!");
        }

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
            throw new Exception("Node A was not found in any tree!");
        }

        List<Boolean> nBinTs = new ArrayList<>();
        for (NetNode nodeBinT : nodeBinTs) {
            nBinTs.add(nodeBinT != null);
        }
        if (Collections.frequency(nBinTs, true) < 1) {
            throw new Exception("Node B was not found in any tree!");
        }

        boolean violates = false;
        for (int i = 0; i < this._subnetworks.size() - 1; i++) {
            boolean nAinT1 = nAinTs.get(i);
            boolean nBinT1 = nBinTs.get(i);
            NetNode nodeAinT1 = nodeAinTs.get(i);
            NetNode nodeBinT1 = nodeBinTs.get(i);

            for (int j = i+1; j < this._subnetworks.size(); j++) {
                boolean nAinT2 = nAinTs.get(j);
                boolean nBinT2 = nBinTs.get(j);
                NetNode nodeAinT2 = nodeAinTs.get(j);
                NetNode nodeBinT2 = nodeBinTs.get(j);

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
                            Network n1 = _subnetworks.get(i).clone();
                            Network n2 = _subnetworks.get(j).clone();
                            NetNode nA = getNodeFromClade(n1, cladeA);
                            NetNode nB = getNodeFromClade(n2, cladeB);
//                            Map<NetNode, NetNode> old2new = new HashMap<>();
//                            old2new.put(nodeAinT1, nA);
//                            old2new.put(nodeBinT2, nB);
                            List<Network> networklist = joinNodesInBothNets(n1, nA, new ArrayList<>(cladeA), n2, nB, new ArrayList<>(cladeB), true);
                            if (n1 !=  null){
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
                            Network n1 = _subnetworks.get(i).clone();
                            Network n2 = _subnetworks.get(j).clone();

                            NetNode nB = getNodeFromClade(n1, cladeB);
                            NetNode nA = getNodeFromClade(n2, cladeA);
                            joinNodesInBothNets(n1, nB, new ArrayList<>(cladeB), n2, nA, new ArrayList<>(cladeA), true);
                            if (n1 !=  null){
                                violates = (!areNetsCompatible(n1, n2));//TODO
                            }
                        }
                        else if (nBinT2){
                             // Case 9: Node B in T2 only
                            // Only valid if (nodeA, nodeB) are siblings in T2
                            // Do nothing (except update node set)
                        }

                        else{
                            throw new Exception("Node B not found in either network!\n");
                        }

                    } else {
                        throw new Exception("Node A not found in either tree!");
                    }

                }

            }

            if (violates) {
                return violates;
            }
        }

        return violates;
    }


    private static void initNodeHeightMap(Network network) {
//        Map<NetNode, Double> networkHeights = new HashMap<NetNode, Double>();
        for(Object nodeO: Networks.postTraversal(network)){
            NetNode node = (NetNode)nodeO;
            double height = 0.0;
            if(node.isLeaf()){
                NetNodeInfo ni = new NetNodeInfo(height);
                node.setData(ni);
            }
            else{
                for(Object childO: node.getChildren()){
                    NetNode<NetNodeInfo> childNode = (NetNode<NetNodeInfo>) childO;
                    System.out.println(childNode.getData().getHeight());
                    double ht = childNode.getParentDistance(node) + childNode.getData().getHeight();
                    height = Math.max(height, ht);
                    NetNodeInfo ni = new NetNodeInfo(height);
                    node.setData(ni);
                }
            }


        }
    }


    public List<Network> mergeNetsViaNJ() {

        // Check trees are on disjoint leaf sets
        for (int i = 0; i < _leaves.size() - 1; i++) {
            for (int j = i + 1; j < _leaves.size(); j++) {
                Set<String> intersection = new HashSet<String>(_leaves.get(i));
                intersection.retainAll(_leaves.get(j));
                if (!intersection.isEmpty()) {
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

        for (Network<NetNodeInfo> net: _subnetworks){
            initNodeHeightMap(net);
        }

        int n = _num_taxa;
        int mergedNetworks = -1;

        while(n > 1){

            System.out.printf("%d joins to go!\n", n);
            if (n == 4){
                System.out.println("break points 2");
            }

            // Sort the Q-matrix
            List<Tuple<Integer, Integer>> pairs = new ArrayList<>();
            List<Double> qvalues = new ArrayList<>();
            for (int idx1 = 0; idx1 < nodePool.size() - 1; idx1++) {
                NetNode<NetNodeInfo> nd1 = nodePool.get(idx1);
                for (int idx2 = idx1 + 1; idx2 < nodePool.size(); idx2++ ) {
                    NetNode<NetNodeInfo> nd2 = nodePool.get(idx2);
                    double v1 = (n - 2) * nd1.getData().getNjDistances().get(nd2.getName());
                    double qvalue = v1 - nd1.getData().getNJXsub() - nd2.getData().getNJXsub();
                    pairs.add(new Tuple<>(idx1, idx2));
                    qvalues.add(qvalue);
                }
            }
            System.out.println("before:");
            System.out.println(_subnetworks.get(0).toString());
            System.out.println(_subnetworks.get(1).toString());

            // Test for constraint violations
            // TODO: Use multi-threading in test_join function!
            Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> nodesToJoin = null;
            int[] sortedIndices = IntStream.range(0, qvalues.size())
                    .boxed().sorted((i, j) -> Double.compare(qvalues.get(i), qvalues.get(j)))
                    .mapToInt(ele -> ele).toArray();
            for (int idxq : sortedIndices) {
                Tuple<Integer, Integer> pair = pairs.get(idxq);
                NetNode nd1 = nodePool.get(pair.Item1);
                NetNode nd2 = nodePool.get(pair.Item2);

                // Check join does not violate a constraint subnetwork!
                try{
                    boolean violates = testJoin(nd1, nd2); // TODO: Add required parameters
                    if (!violates) {
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
            System.out.println(nodesToJoin.Item1.getName()+";"+nodesToJoin.Item2.getName());

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
            if (n == 2){
                System.out.println("break point 1");
            }
            double sig_A = 0;
            double sig_B = 0;
            if (n > 2){
                d_sum /= ((n-2)*2);
                sig_A = d_A_new/2 + d_sum;
                sig_B = d_A_new - sig_A;
                System.out.println(sig_A);
//                if (sig_A == NetNode.NO_DISTANCE){
//                    System.out.println("sig_A no distances");
//                }
            }


            Tuple<List<Boolean>, NetNode<NetNodeInfo>> tuple = joinNodes(nodesToJoin.Item1, nodesToJoin.Item2, sig_A*_scale, sig_B*_scale); // TODO: Add required parameters
            System.out.println("after:");
            System.out.println(_subnetworks.get(0).toString());
            System.out.println(_subnetworks.get(1).toString());
            System.out.println("--------------");
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
                            mergedNetworks = index;
                            System.out.println("equal leaf set:"+net.toString());
                            break;
//                            return net;
                        }

                    }
                    i++;
                }
            }

            if (mergedNetworks >= 0){
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
            // Adjust count

            n--;
        }
        //todo
        List<Network> candidateNetworks = new ArrayList<>();
        if (Networks.hasTheSameTopology(_subnetworks.get(0), _subnetworks.get(1))){
            candidateNetworks.add(_subnetworks.get(0));
        }
        else{
            for (int i = 0; i < _subnetworks.size(); i++) {
                if ( i == mergedNetworks) {
                    continue;
                }
                Network net2 = _subnetworks.get(i);
                if (candidateNetworks.size() == 0){
                    Network net1 = _subnetworks.get(mergedNetworks);
                    System.out.println("match networks:");
                    System.out.println(net2.toString());
                    System.out.println(net1.toString());
                    candidateNetworks = matchNetworks(net2, net1);
                    System.out.println(candidateNetworks);
                }
                else{
                    List<Network> temp = new ArrayList<>();
                    for (Network net1: candidateNetworks){
                        System.out.println("match networks:");
                        System.out.println(net2.toString());
                        System.out.println(net1.toString());
                        temp.addAll(matchNetworks(net2, net1));

                    }
                    candidateNetworks = temp;
                }
            }
        }

        System.out.println("node pool size:"+nodePool.size());
        for (Network net: candidateNetworks){
            System.out.println(net.toString());
        }
        return candidateNetworks;
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

        NJMerge merge = new NJMerge(netlist, matrix, taxonlist);
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

        NJMerge merge = new NJMerge(netlist, matrix, taxonlist);
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
        Network net1 = Networks.readNetwork("((F:3.0,((D,E))#H1:1.0):2.0,((C,(B,A)),#H1:2.0):1.0);");
        getNodeLeafMap(net1.getRoot(), node2leaf, null);
        for(NetNode node: node2leaf.keySet()){
            System.out.println(node+":"+node2leaf.get(node));
        }
    }

    public static void testGetParent(){
        Network net1 = Networks.readNetwork("((E:1.0,D:1.0):3.0,((B:1.0)#H1:2.0,(C:2.0,#H1:1.0):1.0):1.0);");
        initNodeHeightMap(net1);
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
        initNodeHeightMap(net1);
        initNodeHeightMap(net2);
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

    public static void testCopyNode(){
        Network net = Networks.readNetwork("(((A:2,(B:1)#H1:1):1,(C:2,#H1:1):1):3,((((D:1,E:1):1)#H2:1,F:3):2,#H2:3):1);");
        Networks.autoLabelNodes(net);
        System.out.println(net.toString());
        initNodeHeightMap(net);
        Map<NetNode, NetNode> old2new = new HashMap<>();
        NetNode copy = copyNodeSubnet(net.getRoot(), old2new);
        System.out.println((new BniNetwork<NetNodeInfo>((BniNetNode<NetNodeInfo>) copy)).toString());
        Set<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> candidateEdgeList = findEdgeInRange(net, 1, "F", false);
        for (Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge:candidateEdgeList){
            System.out.println(edge.Item1.getName()+","+edge.Item2.getName());
        }
        System.out.println("------");
        candidateEdgeList = findEdgeInRange(net, 1, "F", true);
        for (Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge:candidateEdgeList){
            System.out.println(edge.Item1.getName()+","+edge.Item2.getName());
        }
    }
    public static void main(String[] args) {
//        testGetParent();
//        testMapSubNets();
        testCase1();
//        testMatchNetworks();
//        testCase2();
//        findEdge();
//        testGetNodeLeafMap();
//        testNetworkRooting();
//        testRooting();
//        testCopyNode();
    }
}
