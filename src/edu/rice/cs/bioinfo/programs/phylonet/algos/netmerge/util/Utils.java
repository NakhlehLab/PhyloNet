package edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.util;
/*
 * @ClassName:   Utils
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        9/12/23 2:17 PM
 */

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import sun.nio.ch.Net;

import java.util.*;

import static edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core.NetworkPrior.getParentsList;

public class Utils {
    public static long _SEED = 12345678;
    public static boolean _debug = false;
    public static double _compatible_max_value = 100;
    public static boolean _USEQ = true;
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
        root.setName(node.getName());
        root.setData(new NetNodeInfo());
        if (!node.isLeaf()){

            root.setName(node.getName());
            root.setData(new NetNodeInfo());
            for (Object o: node.getChildren()) {
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
//    }

    public static boolean ReticulationInDescendants(NetNode node){
        if (node.isLeaf()){
            return false;
        }
        if (node.isNetworkNode()){
            return true;
        }
        for (Object o: node.getChildren()){
            if (ReticulationInDescendants((NetNode) o)){
                return true;
            }
        }
        return false;
    }

    public static void getNodeLeafMap(NetNode node, Map<NetNode, Set<String>> nodeLeafMap, Set<String> restrictedSet){
        Set<String> leafset = new HashSet<>();
//        if (node.getName().equals("I1")) {
//            System.out.println("checkpoint");
//        }
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

    public static NetNode commonParents(NetNode nodeA, NetNode nodeB, Set<NetNode> abandonlist) {
        for(Object o1: nodeA.getParents()) {
            NetNode parentA = (NetNode) o1;

            for (Object o2: nodeB.getParents()) {
                NetNode parentB = (NetNode) o2;
                if (parentA.equals(parentB)) {
                    return parentA;
                }
                if (parentB.isNetworkNode() || abandonlist.contains(parentB)){

                    NetNode commonParent = commonParents(nodeA, parentB, abandonlist);
                    if (commonParent != null){
                        return commonParent;

                    }
                }
            }
            if (parentA.isNetworkNode() || abandonlist.contains(parentA)){
                NetNode commonParent = commonParents(parentA, nodeB, abandonlist);
                if (commonParent != null){
                    return commonParent;
                }

            }
        }
        return null;
    }

    public static Set<NetNode> getAbandonList(NetNode node){
        Set<NetNode> abandonList = new HashSet<>();
        for (Object o: node.getChildren()){
            NetNode child = (NetNode) o;
            if (child.isNetworkNode()){
                for (Object o1: child.getParents()){
                    NetNode otherparent = (NetNode) o1;
                    if (!otherparent.equals(node)){
                        abandonList.add(otherparent);
                    }
                }
                abandonList.addAll(getAbandonList(child));
            }
        }
        return abandonList;
    }

    public static void removeInternalLeaf(List<Network> results){
        for (Network<String> network: results){
            for (NetNode<String> node: network.getLeaves()){
                if (node.getName().isEmpty() || node.getName().startsWith("I")){
                    node.removeItself();
                }
            }
            Networks.removeBinaryNodes(network);
        }
    }

    public static void addInheritanceProb(Network network){
        for (Object o: network.getNetworkNodes()){
            NetNode netnode = (NetNode) o;
            for(Object o1: netnode.getParents()){
                NetNode netNodeParent = (NetNode) o1;
                netnode.setParentDistance(netNodeParent, 1.0);
                netnode.setParentProbability(netNodeParent, 0.5);
            }
        }
    }



    public static boolean isAtReticulationNode(NetNode node){
        if(node.isNetworkNode()){
            return true;
        }
        NetNode parent = (NetNode)(node.getParents().iterator().next());
        if (parent.isRoot()){
            return false;
        }
        return isAtReticulationNode(parent);
    }

    public static void removeRetiOverOutgroup(NetNode outgroup){
        for (Object o: outgroup.getParents()){
            NetNode parent = (NetNode) o;
            if (parent.isNetworkNode()){
                NetNode grandparent = (NetNode) parent.getParents().iterator().next();
                grandparent.removeChild(parent);
            }
            removeRetiOverOutgroup(parent);
        }
    }

//    public static void main(String[] args) {
//        Network net = Networks.readNetwork("((((t21:0.689109324)I9#H1:1.0::0.5,(t22:0.5326000818,t10:0.5326000818)I8:0.1565092422)I5:4.114409669,((t33:1.432442835,(t16:0.4712545803,(t25:0.1626650406,t30:0.1626650406)I11:0.3085895397)I10:0.9611882551)I6:2.7344783452)I3#H2:1.0::0.5)I2:1.9346975199099998,((t11:2.87898364,(t38:0.689109324,I9#H1:1.0::0.5)I7:2.189874316)I4:1.905507974,I3#H2:1.0::0.5)I1:1.9537249)I0;(((t4:3.392681616,(((((t6:0.002301929895,t20:0.002301929895)I19:0.006815844532,t32:0.009117774426)I18:0.3733282863,t29:0.3824460607)I16:0.5613405463,(t13:0.4412615513,t24:0.4412615513)I15:0.5025250558)I14:1.897741706,(t9:0.1334509069,t1:0.1334509069)I13:2.708077406)I12:0.5511533024)I11:6.347379667,(((t38,t11)I5,(((((t30,t25)I2,t16)I3,t33)I4,(t21,(t10,t22)I1)I17)I6,t2)I7)I8,t18:8.056923686000001)I10:1.683137597)I9,Z:100.0)I0;");
//        NetNode node = net.findNode("I1");
//        Map<NetNode, Set<String>> node2leaf1 = new HashMap<>();
//        Utils.getNodeLeafMap(net.getRoot(), node2leaf1, null);
//        System.out.println(node2leaf1.get(node));
//
//    }

    public static void main(String[] args) {
        Network net = Networks.readNetwork("(((((((t6:0.2353460826,t86:0.2353460826)I51:3.619316958,(((((t119:0.8344322293,t112:0.8344322293)I56:0.9382349602,t84:1.77266719)I57:0.7383950385,t83:2.511062228)I60:0.4692187249,((t51:1.998570766,t31:1.998570766)I58:0.3513426429,t99:2.349913409)I59:0.6303675438)I61:0.7390462469,(((t95:0.2082478511,t106:0.2082478511)I53:0.3499951833,t36:0.5582430344)I54:0.9137831323,(t40:0.0471266567,t22:0.0471266567)I52:1.42489951)I55:2.247301033)I63:0.1353358407)I64:0.9243357015,((t18:0.1862424469,t43:0.1862424469)I50:3.471568874,t62:3.657811321)I62:1.121187421)I65:3.113922607,((((t8:0.9087307318,(t110:0.3928876034,t29:0.3928876034)I43:0.5158431284)I44:1.464101747,t107:2.372832479)I46:1.149154701,(t53:1.428409981,t57:1.428409981)I45:2.093577199)I47:0.9984538019,(t39:3.573919069,t41:3.573919069)I48:0.9465219131)I49:3.372480367)I66:0.758856714,((((t79:0.6916003336,((t65:0.2456928933)#H2:1.0::0.5,t127:0.2456928933)I39:0.4459074403),Z)I69,((t91:1.636746677,t59:1.636746677)I41:0.5076824401,((t71:0.02710812164,t76:0.02710812164)I38:0.5626914415,t80:0.5897995631)I40:1.554629554)I42:0.2613761506)I70:0.4850631353)#H1:1.0::0.5):5.445982798,((((t121,t23)I20,((t118,t46)I13,t89)I21)I22,((((t17,t38)I16,t77)I17,t4)I34,((((t9:1.423452137,((t50:0.869166171,((t52:0.266739528,t78:0.266739528)I24:0.08528160935,t68:0.3520211373)I25:0.5171450336)I26:0.1137464698,t73:0.9829126408)I27:0.4405394964)I28:3.050170956,((t33:0.7122115354,t104:0.7122115354)I23:1.040462068,t61:1.752673603)I29:2.72094949)I31:0.6485129492,t67:5.122136042)I32:0.6525161646,((t74:0.4601377482,t109:0.4601377482)I19:2.395272024,t54:2.855409772)I30:2.919242434)I33)I35)I37,#H1:1.0::0.5):6.704588014)I67:1.281042679,((((((t124,t81)I4,(t100,(t72,t30)I5)I9)I10,(t15,t19)I8)I14,(t66,t42)I12)I15,((((t120,t103)I0,t34)I2,((t108,t45)I1,t75)I3)I11,(t70,(t24,t96)I6)I7)I18)I36,#H2:1.0::0.5):14.903641259999999)I68;");
        removeRetiOverOutgroup(net.findNode("Z"));
        System.out.println(net.toString());
        net.resetRoot("Z");
        System.out.println(net.toString());
    }
}
