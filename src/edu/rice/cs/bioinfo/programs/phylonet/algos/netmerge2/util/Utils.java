package edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge2.util;
/*
 * @ClassName:   Utils
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        9/12/23 2:17 PM
 */

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge2.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

import static edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core.NetworkPrior.getParentsList;


public class Utils {
    public static long _SEED = 12345678;
    public static boolean _debug = false;

    /* Constructor */
    public Utils() {

    }

    public static int[] argsort(final List<Double> a, final boolean ascending) {
        Integer[] indexes = new Integer[a.size()];
        for (int i = 0; i < indexes.length; i++) {
            indexes[i] = i;
        }
        Arrays.sort(indexes, new Comparator<Integer>() {
            @Override
            public int compare(final Integer i1, final Integer i2) {
                return (ascending ? 1 : -1) * Double.compare(a.get(i1), a.get(i2));
            }
        });
        return asArray(indexes);
    }

    public static <T extends Number> int[] asArray(final T... a) {
        int[] b = new int[a.length];
        for (int i = 0; i < b.length; i++) {
            b[i] = a[i].intValue();
        }
        return b;
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

    public static NetNode<NetNodeInfo> copyNodeSubnet(NetNode<NetNodeInfo> node, Map<NetNode, NetNode> old2new){

        NetNode<NetNodeInfo> root = new BniNetNode();
        if (old2new.keySet().contains(node)){
            return old2new.get(node);
        }
        if (node.isLeaf()){
            root.setName(node.getName());
            root.setData(new NetNodeInfo());
        }
        else{
            root.setName(node.getName());
            root.setData(new NetNodeInfo());
            for (Object o: node.getChildren()){
                NetNode<NetNodeInfo> child = copyNodeSubnet((NetNode) o, old2new);
                root.adoptChild(child, NetNode.NO_DISTANCE);
            }
        }
        old2new.put(node, root);
        return root;
    }

    public static <V, K> Map<V, List<K>> invertMapUsingGroupingBy(Map<K, V> map) {
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

    public static <T> void emptyNodeLabel(Network<T> net)
    {
//        Map<String, NetNode<T>> map = new Hashtable<String, NetNode<T>>();
//        for (NetNode<T> node : net.bfs()) {
//            if (!node.getName().equals(NetNode.NO_NAME)) {
//                map.put(node.getName(), node);
//            }
//        }

        for (NetNode<T> node : net.bfs()) {
            if (!node.isLeaf()) {
                String name = "";
//                do {
//                    name = "";
//                } while (map.get(name) != null);
                node.setName(name);
            }
        }
    }
//
//    public void rerootNetworkAtEdge(String nodeName){
//        rerootNetworkAtEdge((nodeName));
//    }
//
//
//    /**
//     * Reroot this Network at edge incident with node <code>node</code>
//     */
//    public void rerootNetworkAtEdge(NetNode node){
//        doRerooting(node.getParents().iterator().next());
//        List<NetNode> siblinglist = ((STINode)node).getSiblings();
//        STINode newnode = _root.createChild();
//        for(Object o: siblinglist){
//            newnode.adoptChild((TMutableNode)o);
//        }
//        Networks.removeBinaryNodes(this);
//    }
//
//
//    /**
//     * Reroot this Network at node <code>node</code>
//     */
//    public static void rerootNetworkAtNode(NetNode node){
//        if(!_node_set.contains(node)){
//            throw new RuntimeException("node " + node + " is not in the Network "+ this.toNewick());
//        }
//        if(node.isRoot()){
//            return;
//        }
//        if(node.isLeaf()){
//            rerootNetworkAtEdge(node);
//        }
//        else{
//            doRerooting(node);
//        }
//        Networks.removeBinaryNodes(this);
//    }
//
//
//    /**
//     * Reroot this Network at node <code>node</code>
//     */
//    private void doRerooting(NetNode node){
//        STINode parent = (STINode)(node.getParent());
//        if(parent == null){
//            return;
//        }
//        if(!parent.isRoot()){
//            doRerooting(node.getParent());
//        }
//        parent._children.remove(node);
//        if(parent.getChildCount()==1){
//            ((STINode)node).adoptChild(parent);
//            ((STINode)node).removeChild(parent, true);
//        }
//        else{
//            ((STINode)node).adoptChild(parent);
//        }
//
//    }

    public static Map<String, Integer> calculateLevels(Network network) {
        HashMap<String, Integer> distanceMap = new HashMap<>();
        calculateLevelsUtil(network.getRoot(), 0, distanceMap);

        return distanceMap;
    }

    private static void calculateLevelsUtil(NetNode node, int level, HashMap<String, Integer> distanceMap) {
        if (node == null || distanceMap.containsKey(node.getName())) return;

        // Update the distance for this node
        if (node.getName() != ""){
            distanceMap.put(node.getName(), level);
        }


        // Process all children
        for (Object o : node.getChildren()) {
            NetNode child = (NetNode) o;
            calculateLevelsUtil(child, level + 1, distanceMap);
        }
    }

    public static Set<String> getNodeSet(Network network){
        Set<String> nodeSet = new HashSet<>();
        for (Object o: network.bfs()){
            NetNode node = (NetNode) o;
            if(node.getName() != null && node.getName() != ""){
                nodeSet.add(node.getName());
            }

        }
        return nodeSet;
    }

    public static List<NetNode> findPathRecursive(NetNode node, List<NetNode> path){
        if (node.isRoot()) {
            path.add(node);
            return path;
        }
        path.add(node);
        List<NetNode> shortest = null;
        for (Object o : node.getParents()) {
            NetNode parent = (NetNode) o;
            List<NetNode> newPath = new ArrayList<>(path);
            List<NetNode> tmp = findPathRecursive(parent, newPath);
            if (shortest == null || (tmp.size() < shortest.size())) {
                shortest = tmp;
            }
        }
        return shortest;
    }


//    public static void getPath2Root(Network<NetNodeInfo> net, String nodename, Set<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> setEdgeInPath, Map<String, Double> distanceMap){
//        findEdgeInPathRecursive(net.findNode(nodename), setEdgeInPath, distanceMap);
//    }

    public static double getDistanceBetweenTwoNodes(Network network, String node1, String node2) {
        List<NetNode> path1 = findPathRecursive(network.findNode(node1), new ArrayList<>());

        List<NetNode> path2 = findPathRecursive(network.findNode(node2), new ArrayList<>());
//        System.out.println(path1);
//        System.out.println(path2);

        int i = path1.size() - 1;
        int j = path2.size() - 1;

        while (i >= 0 && j >= 0) {
            if (path1.get(i).equals(path2.get(j))) {
                i--;
                j--;

            } else {
                break;
            }
        }
        return i + j + 2;
    }

    // remove reticulation edge (v1, v2)
    public static void removeReticulation(NetNode v1, NetNode v2){
        Iterator itv2p = v2.getParents().iterator();
        NetNode v2p = (NetNode) itv2p.next();
        NetNode v1p = (NetNode) v1.getParents().iterator().next();
        NetNode v2c = (NetNode) v2.getChildren().iterator().next();

        Iterator itv1c = v1.getChildren().iterator();
        NetNode v1c = (NetNode) itv1c.next();
        if (v1c.equals(v2)){
            v1c = (NetNode) itv1c.next();
        }
        double b1 = v1.getParentDistance(v1p)+v1c.getParentDistance(v1);
        double b2 = v2.getParentDistance(v2p)+v2c.getParentDistance(v2);

        if (v2p.equals(v1)){
            v2p = (NetNode) itv2p.next();
        }

//        v2p.removeChild(v2);
        v1p.removeChild(v1);
        v1.removeChild(v1c);
        v1.removeChild(v2);
//        v2p.adoptChild(v2c, b2);
        v1p.adoptChild(v1c, b1);





    }

    public static void getLeavesUnderNode(NetNode node, Set<String> leaves){
        if (node.isLeaf()){
            leaves.add(node.getName());
        }
        else{
            for (Object o: node.getChildren()){
                getLeavesUnderNode((NetNode) o, leaves);
            }
        }
    }

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

    public static Set<String> getLeafSetButReti(NetNode node){
        Set<String> leafset = new HashSet<>();
        if (node.isLeaf()){
            leafset.add(node.getName());

        }
        else if (!node.isNetworkNode()){
            for (Object o: node.getChildren()){
                NetNode child = (NetNode) o;
                leafset.addAll(getLeafSetButReti(child));
            }
        }

        return leafset;

    }


    public static void blankSubNetInternalNodeNames(Network net){
        for (Object o: net.bfs()){
            NetNode node = (NetNode) o;
            if (!node.isLeaf()){
                node.setName("");
            }
        }
    }

//    public static void getOverlapNodes(Network network, Set<String> cladeX, Set<String> cladeB, NetNode nodeB){
//        List<String> leaves = new ArrayList<>();
//        Queue<Tuple<NetNode, String>> queue = new LinkedList<>();
//        queue.add(new Tuple<>(nodeB.getName());
//
//        // TODO
//    }

    public static void getAllBackbonesDfs(int index, int retiLimit, Network cur, List<NetNode> curReticulations, List<Network> results) {
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
            if(cur.getReticulationCount() <= retiLimit) {
                Network temp = cur.clone();
                Networks.removeBinaryNodes(temp);
                results.add(temp);
            }
            getAllBackbonesDfs(index + 1, retiLimit, cur, curReticulations, results);

            parent.adoptChild(curnode, distanceBackup);
        }

        getAllBackbonesDfs(index + 1, retiLimit, cur, curReticulations, results);
    }

    public static List<Network> getAllBackboneNets(Network network, int retiLimit) {
        Network clonedNetwork = network.clone();

        List<NetNode> trueReticulations = new ArrayList<>();
        for(Object nodeObj : clonedNetwork.dfs()) {
            NetNode node = (NetNode) nodeObj;
            if(node.isNetworkNode()) {
                trueReticulations.add(node);
                for(Object parentObj : node.getParents()) {
                    NetNode parent = (NetNode) parentObj;
                    node.setParentProbability(parent, NetNode.NO_PROBABILITY);
                }
            }
        }

        List<Network> results = new ArrayList<>();

        if(network.getReticulationCount() > retiLimit) {

        }

        //added by Zhen, because of reticulation dependency
        Collections.reverse(trueReticulations);

        getAllBackbonesDfs(0, retiLimit, clonedNetwork, trueReticulations, results);
        //added by Zhen, because of reticulation dependency
        removeEmptyLeaf(results);
        return results;
    }

    public static void removeEmptyLeaf(List<Network> results){
        for (Network<String> network: results){
            for (NetNode<String> node: network.getLeaves()){
                if (node.getName().isEmpty()){
                    node.removeItself();
                }
            }
            Networks.removeBinaryNodes(network);
        }
    }


    public static Tuple<List<NetNode>, NetNode> getCommonAncestorList(Network net1, List<String> intersectionABList){
        Set<NetNode> parentIntersection = new HashSet<>();
        NetNode leaf = net1.findNode(intersectionABList.get(0));
        List<NetNode> parentlist = getParentsList(net1, leaf);
        parentIntersection.addAll(parentlist);
        for (int i = 1; i < intersectionABList.size()-1; i++) {
            NetNode leaf2 = net1.findNode(intersectionABList.get(i));
            List<NetNode> parentlist2 = getParentsList(net1, leaf2);
            parentIntersection.retainAll(parentlist2);

        }
        NetNode leaflast = net1.findNode(intersectionABList.get(intersectionABList.size()-1));
        List<NetNode> parentlistlast = getParentsList(net1, leaflast);
        List<NetNode> commonAncestorList = new ArrayList<>();
        NetNode firsthybrid = null;
        for (NetNode parent: parentlistlast){
            if (parentIntersection.contains(parent)){
                commonAncestorList.add(parent);
                if(parent.isNetworkNode()){
                    firsthybrid = parent;
                }
            }
        }
        return new Tuple<>(commonAncestorList, firsthybrid);

    }
//    public static void
//    public static NetNode getMRCA(Network net, Set<String> clade){
//        NetNode mrca = null;
//        for (String leafname: clade){
//            NetNode leaf = net.findNode(leafname);
//
//        }
//        return mrca;
//    }'
    public NetNode commonParents(NetNode nodeA, NetNode nodeB) {
        for(Object o1: nodeA.getParents()) {
            NetNode parentA = (NetNode) o1;
            for (Object o2: nodeB.getParents()) {
                NetNode parentB = (NetNode) o2;
                if (parentA.equals(parentB)) {
                    return parentA;
                }
            }
        }
        return null;
    }

}
