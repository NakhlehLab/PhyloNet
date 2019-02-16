package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

import static com.sun.tools.javac.jvm.ByteCodes.swap;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 5/31/18
 * Time: 9:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class SuperNetwork2 {
    private double MCRAs_[][];
    private int n_;
    private List<Network> subnetworks_;
    private List<String> leafnames_;
    private Network backbone_;
    private Map<NetNode, Double> heightMap_ = new HashMap<>();
    private Map<NetNode, NetNode> resolved_ = new HashMap<>();

    SuperNetwork2(List<Network> subnetworks) {
        subnetworks_ = new ArrayList<>();
        for(Network net : subnetworks) {
            subnetworks_.add(Networks.readNetwork(net.toString()));
        }

        Set<String> leafnames = new HashSet<>();
        for(Network net : subnetworks_) {
            for(Object leafObj : net.getLeaves()) {
                NetNode leaf = (NetNode) leafObj;
                leafnames.add(leaf.getName());
            }
        }
        leafnames_ = new ArrayList<>();
        leafnames_.addAll(leafnames);
        Collections.sort(leafnames_);

        n_ = leafnames_.size();
        MCRAs_ = new double[n_][n_];
        for(int i = 0 ; i < n_ ; i++) {
            for(int j = 0 ; j < n_ ; j++) {
                MCRAs_[i][j] = Double.MAX_VALUE;
            }
        }

        for(Network net : subnetworks_) {
            for (Object nodeObj : Networks.postTraversal(net)) {
                NetNode node = (NetNode) nodeObj;
                if (node.isLeaf()) {
                    heightMap_.put(node, 0.0);
                }
                // store height of parents
                for (Object parentObj : node.getParents()) {
                    NetNode parent = (NetNode) parentObj;
                    double parentHeight = node.getParentDistance(parent) + heightMap_.get(node);
                    // check if network is ultrametric
                    if (heightMap_.containsKey(parent)) {
                        if (Math.abs(heightMap_.get(parent) - parentHeight) > 1e-6) {
                            System.out.println("Input network is not ultrametric!");
                            System.out.println(net.toString());
                            System.exit(1);
                        }
                    } else {
                        heightMap_.put(parent, parentHeight);
                    }
                }
            }
        }

        for(Network net : subnetworks_) {
            for (Object nodeObj : Networks.postTraversal(net)) {
                NetNode node = (NetNode) nodeObj;
                resolved_.put(node, null);
            }
        }
    }

    static public Map<NetNode, Double> getNodeHeights(Network network) {
        Map<NetNode, Double> heightMap = new HashMap<>();
        for(Object nodeObj : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeObj;
            if(node.isLeaf()) {
                heightMap.put(node, 0.0);
            }
            // store height of parents
            for(Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;
                double parentHeight = node.getParentDistance(parent) + heightMap.get(node);
                // check if network is ultrametric
                if(heightMap.containsKey(parent)) {
                    if(Math.abs(heightMap.get(parent) - parentHeight) > 1e-6) {
                        System.out.println("Input network is not ultrametric!");
                        System.out.println(network.toString());
                        System.exit(1);
                    }
                } else {
                    heightMap.put(parent, parentHeight);
                }
            }
        }
        return heightMap;
    }

    // from parent to child
    static boolean pathExist(NetNode n1, NetNode n2) {
        if(n1 == n2) {
            return true;
        }
        Queue<NetNode> queue = new LinkedList<>();
        queue.add(n1);
        while(!queue.isEmpty()) {
            NetNode node = queue.poll();
            for(Object childObj : node.getChildren()) {
                NetNode child = (NetNode) childObj;
                if(child == n2) {
                    return true;
                }
                queue.add(child);
            }
        }
        return false;
    }

    static List<NetNode> getPath(NetNode n1, NetNode n2) {
        List<NetNode> result = new ArrayList<>();
        if(n1 == n2) {
            return result;
        }
        Queue<NetNode> queue = new LinkedList<>();
        Map<NetNode, NetNode> precursor = new HashMap<>();
        queue.add(n1);
        while(!queue.isEmpty()) {
            NetNode node = queue.poll();
            for(Object childObj : node.getChildren()) {
                NetNode child = (NetNode) childObj;
                precursor.put(child, node);
                if(child == n2) {
                    node = precursor.get(n2);
                    while(node != n1) {
                        result.add(0, node);
                        node = precursor.get(node);
                    }
                    return result;
                }
                queue.add(child);
            }
        }
        return null;
    }

    static public List<String> getLeafList(Network net) {
        List<String> result = new ArrayList<>();
        for(Object leafObj : net.getLeaves()) {
            NetNode leaf = (NetNode) leafObj;
            result.add(leaf.getName());
        }
        Collections.sort(result);
        return result;
    }

    static public List<String> getLeafListUnder(NetNode node) {
        Set<String> leaves = new HashSet<>();
        Queue<NetNode> q = new ArrayDeque<>();
        q.add(node);
        while(!q.isEmpty()) {
            NetNode n = q.poll();
            for(Object childObj : n.getChildren()) {
                NetNode child = (NetNode) childObj;
                if(child.isLeaf()) {
                    leaves.add(child.getName());
                } else {
                    q.add(child);
                }
            }
        }
        List<String> result = new ArrayList<>();
        result.addAll(leaves);
        Collections.sort(result);
        return result;
    }


    public Tuple<Network, Map<NetNode, NetNode>> getSubNetwork(Network trueNetwork, List<String> selectedLeaves, boolean removeBinaryNodes) {
        int subsize = selectedLeaves.size();
        Network currentNetwork = trueNetwork;

        // bottom-up build subnetwork
        Map<NetNode, BniNetNode> old2new = new HashMap<>();
        Map<NetNode, Integer> hitCount = new HashMap<>();
        for(String leaf : selectedLeaves) {
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
        if(removeBinaryNodes) {
            Networks.removeBinaryNodes(newsubnetwork);
        }

        Map<NetNode, NetNode> new2old = new HashMap<>();
        for(NetNode oldNode : old2new.keySet()) {
            new2old.put(old2new.get(oldNode), oldNode);
        }



        return new Tuple<>(newsubnetwork, new2old);
    }

    static Set<NetNode> getAllAncestors(NetNode node) {
        Queue<NetNode> q = new ArrayDeque<>();
        Set<NetNode> result = new HashSet<>();
        q.add(node);
        while(!q.isEmpty()) {
            NetNode curnode = q.poll();
            for(Object parentObj : curnode.getParents()) {
                NetNode parent = (NetNode) parentObj;
                result.add(parent);
                q.add(parent);
            }
        }
        return result;
    }

    NetNode getMostRecentRetiAncestor(NetNode node) {
        Set<NetNode> allAncestors = getAllAncestors(node);
        Map<NetNode, Double> bHeights = getNodeHeights(backbone_);
        double bestHeight = 1e99;
        NetNode bestNode = null;
        if(bHeights.get(node) == null) {
            bHeights = heightMap_;
        }
        for(NetNode anc : allAncestors) {
            if(anc.isNetworkNode() && bHeights.get(anc) < bestHeight) {
                bestHeight = bHeights.get(anc);
                bestNode = anc;
            }
        }
        return bestNode;
    }

    double minMCRA(Network net, NetNode node1, NetNode node2) {
        Set<NetNode> ac1 = getAllAncestors(node1);
        Set<NetNode> ac2 = getAllAncestors(node2);
        Set<NetNode> ca = new HashSet<>(ac1);
        ca.retainAll(ac2); // intersection
        double time = Double.MAX_VALUE;
        for(NetNode node : ca) {
            time = Math.min(time, heightMap_.get(node));
        }
        return time;
    }

    void computeMCRAs() {
        for(Network net : subnetworks_) {
            for(Object leafObj1 : net.getLeaves()) {
                for(Object leafObj2 : net.getLeaves()) {
                    if(leafObj1 == leafObj2) {
                        continue;
                    }
                    NetNode leaf1 = (NetNode) leafObj1;
                    NetNode leaf2 = (NetNode) leafObj2;
                    int i1 = leafnames_.indexOf(leaf1.getName());
                    int i2 = leafnames_.indexOf(leaf2.getName());
                    MCRAs_[i1][i2] = Math.min(minMCRA(net, leaf1, leaf2), MCRAs_[i1][i2]);

                }
            }
        }
    }

    static List<Tuple<NetNode, NetNode>> getIntersections(Network net, double height) {
        Map<NetNode, Double> heights = getNodeHeights(net);
        List<Tuple<NetNode, NetNode>> results = new ArrayList<>();
        for(Object nodeObj : Networks.postTraversal(net)) {
            NetNode node = (NetNode) nodeObj;
            for(Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;
                if(heights.get(parent) > height && heights.get(node) < height) {
                    results.add(new Tuple<>(node, parent));
                }
            }
        }
        return results;
    }

    public static List<NetNode> getParentsList(Network net, NetNode node) {
        List<NetNode> list = new ArrayList<>();
        NetNode tmp = node;
        while(tmp != net.getRoot()) {
            list.add(tmp);
            tmp = (NetNode) tmp.getParents().iterator().next();
        }
        return list;
    }

    public static int getDistance(NetNode node, NetNode mrca)  {
        int dist = 0;
        NetNode tmp = node;
        while(!tmp.equals(mrca)) {
            NetNode par = (NetNode) tmp.getParents().iterator().next();
            //dist += tmp.getParentDistance(par);
            dist += 1.0;
            tmp = par;
        }
        return dist;
    }

    public static int getReticulationNodeDiameter(Network net, Tuple<NetNode, NetNode> parents) {
        List<NetNode> list1 = getParentsList(net, parents.Item1);
        List<NetNode> list2 = getParentsList(net, parents.Item2);
        NetNode mrca = net.getRoot();
        for(NetNode n1: list1) {
            for(NetNode n2 : list2) {
                if(n1.equals(n2)) {
                    mrca = n1;
                    break;
                }
            }
        }
        return getDistance(parents.Item1, mrca) + getDistance(parents.Item2, mrca);
    }

    public static Map<NetNode, Integer> getDiameterMap(Network network) {
        Network net = network.clone();
        Map<NetNode, Integer> distMap = new HashMap<NetNode, Integer>();
        Map<NetNode, Tuple<NetNode, NetNode>> parentsMap = new HashMap<>();
        for(Object o : net.getNetworkNodes()) {
            NetNode node = (NetNode) o;
            List<NetNode> parents = IterableHelp.toList(node.getParents());
            if(parents.size() != 2) throw new IllegalArgumentException("Invalid network " + net.toString());
            distMap.put(node, 0);
            parentsMap.put(node, new Tuple<>(parents.get(0), parents.get(1)));
            for(NetNode par : parents) {
                //distMap.put(node, distMap.get(node) + node.getParentDistance(par));
                distMap.put(node, distMap.get(node) + 1);
                if(node.getParentProbability(par) < 0.50) {
                    par.removeChild(node);
                }
            }
        }
        for(NetNode key: distMap.keySet()) {
            distMap.put(key, distMap.get(key) + getReticulationNodeDiameter(net, parentsMap.get(key)));
        }
        return distMap;
    }

    double getClusterHeight(List<String> leaves1, List<String> leaves2) {
        double height = 0.0;
        for(String leaf1 : leaves1) {
            for(String leaf2 : leaves2) {
                int i1 = leafnames_.indexOf(leaf1);
                int i2 = leafnames_.indexOf(leaf2);
                height = Math.max(height, MCRAs_[i1][i2]);
            }
        }
        return height;
    }

    Network buildTree() {
        Set<Tuple3<BniNetNode, List<String>, Double>> nodes = new HashSet<>();
        for(String leafname : leafnames_) {
            BniNetNode node = new BniNetNode();
            node.setName(leafname);
            List<String> list = new ArrayList<>();
            list.add(leafname);
            nodes.add(new Tuple3<>(node, list, 0.0));
        }

        while(nodes.size() > 1) {
            double bestTime = 0.0;
            Tuple3<BniNetNode, List<String>, Double> bestTuple1 = null;
            Tuple3<BniNetNode, List<String>, Double> bestTuple2 = null;

            for (Tuple3<BniNetNode, List<String>, Double> tuple1 : nodes) {
                for (Tuple3<BniNetNode, List<String>, Double> tuple2 : nodes) {
                    if(tuple1 == tuple2) {
                        continue;
                    }
                    double time = getClusterHeight(tuple1.Item2, tuple2.Item2);
                    if(bestTuple1 == null || bestTime > time) {
                        bestTuple1 = tuple1;
                        bestTuple2 = tuple2;
                        bestTime = time;
                    }
                }
            }

            BniNetNode node = new BniNetNode();
            node.adoptChild(bestTuple1.Item1, bestTime - bestTuple1.Item3);
            node.adoptChild(bestTuple2.Item1, bestTime - bestTuple2.Item3);
            List<String> list = new ArrayList<>();
            list.addAll(bestTuple1.Item2);
            list.addAll(bestTuple2.Item2);
            nodes.remove(bestTuple1);
            nodes.remove(bestTuple2);
            nodes.add(new Tuple3<>(node, list, bestTime));

        }

        return new BniNetwork(nodes.iterator().next().Item1);
    }

    // Map nodes of trinets to nodes on backbone
    void resolveTreeNodes() {
        int sum = 0;

        for(Network net : subnetworks_) {
            for (Object nodeObj : Networks.postTraversal(net)) {
                NetNode node = (NetNode) nodeObj;
                if(node.isLeaf()) {
                    resolved_.put(node, backbone_.findNode(node.getName()));
                }
            }
        }

        for(NetNode node : resolved_.keySet()) {
            if(resolved_.get(node) != null && !getNodeHeights(backbone_).containsKey(resolved_.get(node))) {
                throw new RuntimeException("!!!");
            }
        }

        for(Network net : subnetworks_) {
            System.out.println(Arrays.toString(getLeafList(net).toArray()));

            Tuple<Network, Map<NetNode, NetNode>> rTuple = getSubNetwork(backbone_, getLeafList(net), true);
            Network restrictedBackbone = rTuple.Item1;
            Map<NetNode, NetNode> rNew2Old = rTuple.Item2;
            Map<NetNode, Double> restrictedHeights = getNodeHeights(restrictedBackbone);
            for(NetNode rNode : restrictedHeights.keySet()) {
                if(restrictedHeights.get(rNode) > 0)
                for (Object nodeObj : Networks.postTraversal(net)) {
                    NetNode node = (NetNode) nodeObj;
                    if(heightMap_.get(node) > restrictedHeights.get(rNode) - 1e-3) {
                        if(heightMap_.get(node) < restrictedHeights.get(rNode) + 1e-3) {
                            List<String> l1 = getLeafListUnder(node);
                            List<String> l2 = getLeafListUnder(rNode);
                            if(l2.containsAll(l1) && resolved_.get(node) == null) {
                                resolved_.put(node, rNew2Old.get(rNode));
                                if(!getNodeHeights(backbone_).containsKey(resolved_.get(node))) {
                                    throw new RuntimeException("!!!");
                                }
                                System.out.println("Resolved: " + heightMap_.get(node) + " -> " + restrictedHeights.get(rNode));
                                sum++;
                                break;
                            }
                        }
                    }
                }
            }

        }

        for(NetNode node : resolved_.keySet()) {
            if(resolved_.get(node) != null && !getNodeHeights(backbone_).containsKey(resolved_.get(node))) {
                throw new RuntimeException("!!!");
            }
        }

        System.out.println(sum);
    }

    // Map unresolved tree nodes of trinets
    void resolvedTreeNodes2() {
        while(true) {
            Map<NetNode, Double> bHeights = getNodeHeights(backbone_);

            NetNode pNode = null;
            NetNode ppNode = null;
            double pHeight = 1e99;
            double pDist = 1e99;
            for(Network net : subnetworks_) {
                for (Object nodeObj : Networks.postTraversal(net)) {
                    NetNode node = (NetNode) nodeObj;
                    if(resolved_.get(node) != null && (node.isTreeNode() || node.isLeaf()) && !node.isRoot()) {
                        NetNode parent = (NetNode) node.getParents().iterator().next();

                        if(parent != null
                                && parent.isTreeNode()
                                && resolved_.get(parent) == null
                                && pDist > heightMap_.get(parent) - heightMap_.get(node)) {
                            pHeight = heightMap_.get(parent);
                            pDist = heightMap_.get(parent) - heightMap_.get(node);
                            pNode = node;
                            ppNode = parent;
                        }
                    }
                }
            }

            if(pNode == null) {
                break;
            }

            NetNode bNode = resolved_.get(pNode);
            NetNode bpNode = (NetNode) bNode.getParents().iterator().next();
            if(heightMap_.get(ppNode) + 1e-3 < bHeights.get(bpNode)) {
                BniNetNode newnode = new BniNetNode();
                bpNode.removeChild(bNode);
                newnode.adoptChild(bNode, heightMap_.get(ppNode) - bHeights.get(bNode));
                bpNode.adoptChild(newnode, bHeights.get(bpNode) - heightMap_.get(ppNode));
                resolved_.put(ppNode, newnode);
            } else {
                Set<NetNode> a = getAllAncestors(bNode);
                pDist = 1e99;
                bpNode = null;
                for(NetNode node : a) {
                    if(pDist > Math.abs( bHeights.get(node) - heightMap_.get(ppNode))) {
                        pDist = Math.abs( bHeights.get(node) - heightMap_.get(ppNode));
                        bpNode = node;
                    }
                }
                resolved_.put(ppNode, bpNode);
            }

        }
    }

    boolean resolveNetNode() {
        List<NetNode> netNodes = new ArrayList<>();
        Map<NetNode, Network> node2network = new HashMap<>();
        for(Network net : subnetworks_) {
            Map<NetNode, Integer> distMap = getDiameterMap(net);
            boolean good = true;
            for(NetNode node : distMap.keySet()) {
                if(distMap.get(node) - 1 == 2 || distMap.get(node) - 1 == 3) {
                    good = false;
                }
            }
            if(!good) {
                continue;
            }

            for (Object nodeObj : Networks.postTraversal(net)) {
                NetNode node = (NetNode) nodeObj;
                if(node.isNetworkNode() && resolved_.get(node)==null) {
                    netNodes.add(node);
                    node2network.put(node, net);
                }
            }
        }
        Collections.sort(netNodes, (NetNode n1, NetNode n2)->heightMap_.get(n2).compareTo(heightMap_.get(n1)));
        for(NetNode netNode : netNodes) {
            System.out.println(heightMap_.get(netNode));
            Iterator it = netNode.getParents().iterator();
            double nheight = heightMap_.get(netNode);
            NetNode parent1 = (NetNode) it.next();
            NetNode parent2 = (NetNode) it.next();
            NetNode child = (NetNode) netNode.getChildren().iterator().next();
            NetNode rp1 = resolved_.get(parent1);
            NetNode rp2 = resolved_.get(parent2);
            NetNode rc = resolved_.get(child);
            Map<NetNode, Double> bHeights = getNodeHeights(backbone_);
            List<String> bNames = Networks.getTaxaNamesUnderReticulation(backbone_);
            List<String> sNames = Networks.getTaxaNamesUnderReticulation(node2network.get(netNode));
            boolean good = true;
            for(String bName : bNames) {
                if(node2network.get(netNode).findNode(bName)!=null && !sNames.contains(bName)) {
                    good = false;
                }
            }
            if(!good) {
                continue;
            }

            if(rp1!=null && rp2!=null && rc!=null) {


                if(pathExist(rp1, rc) && pathExist(rp2, rc)) {
                    NetNode r = getMostRecentRetiAncestor(rc);
                    if(r == null || !pathExist(rp1, r) || !pathExist(rp2, r)) {
                        continue;
                    }
                    resolved_.put(netNode, r);
                    System.out.println(heightMap_.get(netNode) + " -> "+bHeights.get(r));
                    System.out.println(backbone_);
                    System.out.println(node2network.get(netNode));
                    continue;
                } else if(pathExist(rp1, rc) || pathExist(rp2, rc)) {
                    if(pathExist(rp2, rc)) {
                        NetNode tt = rp1;
                        rp1 = rp2;
                        rp2 = tt;
                    }
                    List<NetNode> path = getPath(rp1, rc);
                    path.add(0, rp1);
                    if(nheight < bHeights.get(rc) || nheight > bHeights.get(rp1)) {
                        rp1 = path.get(0);
                        nheight = (bHeights.get(rc) + bHeights.get(rp1)) / 2.0;
                    } else {
                        path.add(rc);
                        for (int i = 0 ; i < path.size() ; i++) {
                            if(bHeights.get(path.get(i)) > nheight && bHeights.get(path.get(i + 1)) < nheight) {
                                rp1 = path.get(i);
                                break;
                            }
                        }
                    }
                    if(rp2.getChildCount() == 2 || !rc.hasParent(rp1)) {
                        continue;
                    }
                    BniNetNode newnode = new BniNetNode();
                    rp1.removeChild(rc);
                    rp1.adoptChild(newnode, bHeights.get(rp1) - nheight);
                    newnode.adoptChild(rc, nheight - bHeights.get(rc));
                    resolved_.put(netNode, newnode);
                    System.out.println(heightMap_.get(netNode) + " -> "+nheight);
                    rp2.adoptChild(newnode, bHeights.get(rp2) - nheight);
                    System.out.println(backbone_);
                    System.out.println(node2network.get(netNode));
                } else {
                    continue;
                }
            } else if(rp1!=null && rp2!=null && rc==null) {
                List<Tuple<NetNode, NetNode>> intersections = getIntersections(backbone_, nheight);
                for(Tuple<NetNode, NetNode> edge : intersections) {
                    NetNode bc = edge.Item1;
                    NetNode bp = edge.Item2;
                    Set<NetNode> anc = getAllAncestors(bp);
                    anc.add(bp);
                    for(NetNode node : anc) {
                        if(node.isNetworkNode()) {
                            Set<NetNode> anc1 = getAllAncestors(node);
                            if(anc1.contains(rp1) && anc1.contains(rp2)){
                                resolved_.put(netNode, node);
                                System.out.println(heightMap_.get(netNode) + " -> "+bHeights.get(node));
                                System.out.println(backbone_);
                                System.out.println(node2network.get(netNode));
                                break;
                            }
                        }
                    }
                    if(resolved_.get(netNode) != null) {
                        continue;
                    }

                }
            } else{
                continue;
            }
        }


        return true;
    }

    void compute() {
        computeMCRAs();
        backbone_ = buildTree();
        System.out.println(backbone_);
        resolveTreeNodes();
        resolvedTreeNodes2();
        System.out.println(backbone_);
        resolveNetNode();
        System.out.println(backbone_);
    }
}
