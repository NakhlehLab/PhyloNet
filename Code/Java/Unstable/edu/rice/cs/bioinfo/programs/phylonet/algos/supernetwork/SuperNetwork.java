package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.clustering.DataGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.counting.CoalescenceHistoriesCounting;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetworkMetricNakhleh;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import sun.nio.ch.Net;
import org.apache.commons.math3.util.Combinations;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 10/15/17
 * Time: 11:18 AM
 * To change this template use File | Settings | File Templates.
 */
public class SuperNetwork {
    private List<Network> _subnetworks = null;
    private Network _supernetwork = null;
    private Network _trueNetwork = null;
    private List<String> _nodeLabels = null;
    private Map<String, List<String>> _connections = null;
    private double eps = 1e-6;
    private Random _random = new Random(12345);

    public SuperNetwork(List<Network> subnetworks) {
        Long seed = new Random().nextLong();
        //System.out.println("Seed: " + seed);
        _random.setSeed(seed);
        for(Network net : subnetworks) {
            _subnetworks.add(Networks.readNetwork(net.toString()));
        }
    }

    public List getSubNetworks() {
        return _subnetworks;
    }

    public void sortLabels(List<String> labels, boolean reverse) {
        Collections.sort(labels, (s1, s2)-> {
            double t1, t2;
            try {
                t1 = Double.parseDouble(s1);
            } catch (NumberFormatException e) {
                t1 = 0.0;
            }

            try {
                t2 = Double.parseDouble(s2);
            } catch (NumberFormatException e) {
                t2 = 0.0;
            }

            return Double.compare(t1, t2);
        });
        if(reverse) {
            Collections.reverse(labels);
        }
    }

    public boolean correctHeights() {
        List<String> leafLabels = new ArrayList<>();
        for(Network net : _subnetworks) {
            for(Object leafObj : net.getLeaves()) {
                NetNode leaf = (NetNode) leafObj;
                if(!leafLabels.contains(leaf.getName())) {
                    leafLabels.add(leaf.getName());
                }
            }
        }
        Collections.sort(leafLabels);

        List<Set<NetNode>> corrected = new ArrayList<>();

        for(int i = 0 ; i < leafLabels.size() ; i++) {
            for(int j = i + 1 ; j < leafLabels.size() ; j++) {
                for(Network net : _subnetworks) {

                }
            }
        }

        return false;
    }

    public boolean findNodeLabels() {
        _nodeLabels = new ArrayList<>();
        for(Network net : _subnetworks) {
            Map<NetNode, Double> heightMap = new HashMap<>();
            for(Object nodeObj : Networks.postTraversal(net)) {
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
                        if(Math.abs(heightMap.get(parent) - parentHeight) > eps) {
                            System.out.println("Input network is not ultrametric!");
                            System.out.println(net.toString());
                            System.exit(1);
                        }
                    } else {
                        heightMap.put(parent, parentHeight);
                    }
                }
                // rename internal node according to height
                if(!node.isLeaf()) {
                    String newname = String.format("%.6f", heightMap.get(node));
                    node.setName(newname);
                }
                if(!_nodeLabels.contains(node.getName())) {
                    _nodeLabels.add(0, node.getName());
                }
            }
        }
        sortLabels(_nodeLabels, false);
        return true;
    }

    public boolean findConnections() {
        _connections = new HashMap<>();
        for(Network net : _subnetworks) {
            for(Object nodeObj : Networks.postTraversal(net)) {
                NetNode node = (NetNode) nodeObj;
                if(!_connections.containsKey(node.getName())) {
                    _connections.put(node.getName(), new ArrayList<>());
                }

                for(Object parentObj : node.getParents()) {
                    NetNode parent = (NetNode) parentObj;

                    if(!_connections.containsKey(parent.getName())) {
                        _connections.put(parent.getName(), new ArrayList<>());
                    }

                    if(!_connections.get(parent.getName()).contains(node.getName()))
                        _connections.get(parent.getName()).add(node.getName());
                }
            }
        }
        for(String s : _connections.keySet()) {
            sortLabels(_connections.get(s), true);
        }
        return true;
    }

    public boolean pathExist(NetNode n1, NetNode n2) {
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

    public boolean buildNetwork() {

        // construct connected component
        Map<String, BniNetNode> label2newnode = new HashMap<>();
        for(String s : _nodeLabels) {
            BniNetNode newnode = new BniNetNode();
            newnode.setName(s);
            label2newnode.put(s, newnode);

            for(String s2 : _connections.get(s)) {
                BniNetNode newchild = label2newnode.get(s2);
                if(!pathExist(newnode, newchild) && (newchild.getIndeg() != 1 || newchild.getOutdeg() > 0 && newchild.getOutdeg() < 2)) {
                    newnode.adoptChild(newchild, NetNode.NO_DISTANCE);
                }
            }


            List<NetNode> clusterRoots = new ArrayList<>();
            for(String s2 : _connections.get(s)) {
                Queue<NetNode> queue = new LinkedList<>();
                NetNode node = label2newnode.get(s2);
                if(pathExist(newnode, node)) {
                    continue;
                }
                queue.add(node);
                while(!queue.isEmpty()) {
                    node = queue.poll();
                    if(node.getIndeg() == 0 && !clusterRoots.contains(node)) {
                        clusterRoots.add(node);
                    }
                    for(Object parentObj : node.getParents()) {
                        NetNode parent = (NetNode) parentObj;
                        queue.add(parent);
                    }
                }
            }
            for(NetNode newchild : clusterRoots) {
                if(!pathExist(newnode, newchild)) {
                    newnode.adoptChild(newchild, NetNode.NO_DISTANCE);
                }
            }

        }

        // find root
        BniNetNode newroot = null;
        for(String s : label2newnode.keySet()) {
            BniNetNode node = label2newnode.get(s);
            if(node.getIndeg() == 0) {
                if(newroot == null) {
                    newroot = node;
                } else {
                    System.out.println("Multiple roots!");
                    return false;
                }
            }
        }

        Network result = new BniNetwork(newroot);
        System.out.println("C^3 network: " + result.toString());

        // rebuild relationship
        /*List<String> reversedNodeLabels = new ArrayList<>(_nodeLabels);
        Collections.reverse(reversedNodeLabels);
        for(String s : reversedNodeLabels) {
            BniNetNode newnode = label2newnode.get(s);
            if(newnode.getIndeg() <= 1 && newnode.getOutdeg() == 1) {
                for(String s2 : _connections.get(s)) {
                    BniNetNode newchildnode = label2newnode.get(s2);
                    if(newchildnode.hasParent(newnode)) {
                        continue;
                    }
                    if(newchildnode.getIndeg() == 1 && newchildnode.getOutdeg() == 1) {
                        newnode.adoptChild(newchildnode, NetNode.NO_DISTANCE);
                        break;
                    }
                }
            }
        }*/

        List<String> reversedNodeLabels = new ArrayList<>(_nodeLabels);
        Collections.reverse(reversedNodeLabels);
        List<NetNode> binaryNodes = new ArrayList<>();
        for(String s : reversedNodeLabels) {
            BniNetNode newnode = label2newnode.get(s);
            if(newnode.getIndeg() <= 1 && newnode.getOutdeg() == 1) {
                binaryNodes.add(newnode);
            }
        }

        for(Network net : _subnetworks) {
            List<String> leafLabels = new ArrayList<>();
            for(Object leafObj : net.getLeaves()) {
                NetNode leaf = (NetNode) leafObj;
                leafLabels.add(leaf.getName());
            }
            Network newsubnet = getSubNetwork(result, leafLabels, false, false);

            for(Object nodeObj : net.bfs()) {
                NetNode node = (NetNode) nodeObj;
                String label = node.getName();
                NetNode newdummynode = newsubnet.findNode(label);
                NetNode newnode = label2newnode.get(label);

                if(binaryNodes.contains(newnode)) {

                    if(newnode.getIndeg() <= 1 && newnode.getOutdeg() == 1) {
                        for (Object childObj : node.getChildren()) {
                            NetNode child = (NetNode) childObj;
                            String childLabel = child.getName();

                            NetNode newchild = label2newnode.get(childLabel);
                            if(newchild.hasParent(newnode)) {
                                continue;
                            }
                            if(newchild.getIndeg() == 1 && newchild.getOutdeg() == 1) {
                                // Lian le bai lian
                                Network newsubnetcopy = Networks.readNetwork(newsubnet.toString());
                                newsubnetcopy.findNode(label).adoptChild(newsubnetcopy.findNode(childLabel), NetNode.NO_DISTANCE);
                                Networks.removeBinaryNodes(newsubnetcopy);
                                if(newsubnetcopy.findNode(label) == null || newsubnetcopy.findNode(childLabel) == null) {
                                    continue;
                                }

                                newnode.adoptChild(newchild, NetNode.NO_DISTANCE);
                                break;
                            }
                        }
                    }
                }
                /*
                if(newdummynode != null) {
                    for (Object childObj : node.getChildren()) {
                        NetNode child = (NetNode) childObj;
                        String childLabel = child.getName();
                        boolean newchildExist = false;
                        for (Object newchildObj : newdummynode.getChildren()) {
                            NetNode newchild = (NetNode) newchildObj;
                            if (newchild.getName().equals(childLabel)) {
                                newchildExist = true;
                            }
                        }
                        if (!newchildExist) {
                            NetNode newnode = label2newnode.get(label);
                            NetNode newchild = label2newnode.get(childLabel);
                            newnode.adoptChild(newchild, NetNode.NO_DISTANCE);
                        }
                    }
                } else {
                    newdummynode = label2newnode.get(label);
                    for (Object childObj : node.getChildren()) {
                        NetNode child = (NetNode) childObj;
                        String childLabel = child.getName();
                        boolean newchildExist = false;
                        for (Object newchildObj : newdummynode.getChildren()) {
                            NetNode newchild = (NetNode) newchildObj;
                            if (newchild.getName().equals(childLabel)) {
                                newchildExist = true;
                            }
                        }
                        if (!newchildExist) {
                            NetNode newnode = label2newnode.get(label);
                            NetNode newchild = label2newnode.get(childLabel);
                            if(newchild.getOutdeg() > 0 && newchild.getOutdeg() < 2) {
                                if(newchild.getIndeg() == 2) {
                                    double value = Double.parseDouble(label);
                                    for(Object parentObj : newchild.getParents()) {
                                        NetNode parent = (NetNode) parentObj;
                                        if(value < Double.parseDouble(parent.getName())) {
                                            parent.removeChild(newchild);
                                            newnode.adoptChild(newchild, NetNode.NO_DISTANCE);
                                            break;
                                        }
                                    }
                                } else {
                                    newnode.adoptChild(newchild, NetNode.NO_DISTANCE);
                                }
                            }
                        }
                    }
                }*/
            }
        }

        System.out.println("Super-Network: " + result.toString());
        _supernetwork = result;
        return true;

        /*
        //compute in degrees
        Map<String, Integer> indegree = new HashMap<>();
        for(String s1 : _connections.keySet()) {
            if(!indegree.containsKey(s1)) {
                indegree.put(s1, 0);
            }
            for(String s2 : _connections.get(s1)) {
                if(!indegree.containsKey(s2)) {
                    indegree.put(s2, 0);
                }
                indegree.put(s2, indegree.get(s2) + 1);
            }
        }

        Map<String, List<String>> newReversedConnections = new HashMap<>();

        //topological sorting
        String rootName = null;
        while(indegree.size() > 0) {
            String s = null;
            for(String s1 : indegree.keySet()) {
                if(indegree.get(s1) == 0) {
                    s = s1;
                    if(rootName == null) {
                        rootName = s1;
                    }
                }
            }
            if(s == null) {
                System.out.println("Input contains cycles!");
                System.out.println();
                System.exit(1);
            }
            for(String s2 : _connections.get(s)) {
                if(!newReversedConnections.containsKey(s2)) {
                    newReversedConnections.put(s2, new ArrayList<>());
                }
                if(newReversedConnections.get(s2).contains(s)) {
                    newReversedConnections.get(s2).remove(s);
                }
                newReversedConnections.get(s2).add(s);

                indegree.put(s2, indegree.get(s2) - 1);
            }
            indegree.remove(s);
        }

        // filter reversed connections
        for(String s2 : newReversedConnections.keySet()) {
            while(newReversedConnections.get(s2).size() > 2) {
                newReversedConnections.get(s2).remove(0);
            }
        }

        // build network
        Map<String, BniNetNode> label2node = new HashMap<>();
        for(String s : _nodeLabels) {
            label2node.put(s, new BniNetNode());
            label2node.get(s).setName(s);
        }

        for(String s2 : newReversedConnections.keySet()) {
            for(String s1 : newReversedConnections.get(s2)) {
                label2node.get(s1).adoptChild(label2node.get(s2), NetNode.NO_DISTANCE);
            }
        }

        Network result = new BniNetwork(label2node.get(rootName));
        */

    }

    public boolean checkResult() {
        for(Network net : _subnetworks) {
            List<String> leafLabels = new ArrayList<>();
            for(Object leafObj : net.getLeaves()) {
                NetNode leaf = (NetNode) leafObj;
                leafLabels.add(leaf.getName());
            }
            Network newsubnet = getSubNetwork(_supernetwork, leafLabels, true, false);

            for(Object nodeObj : Networks.postTraversal(net)) {
                NetNode node = (NetNode) nodeObj;
                String label = node.getName();
                NetNode newdummynode = newsubnet.findNode(label);
                for(Object childObj : node.getChildren()) {
                    NetNode child = (NetNode) childObj;
                    String childLabel = child.getName();
                    boolean newchildExist = false;
                    for(Object newchildObj : newdummynode.getChildren()) {
                        NetNode newchild = (NetNode) newchildObj;
                        if(newchild.getName().equals(childLabel)) {
                            newchildExist = true;
                        }
                    }
                    if(!newchildExist) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    public boolean checkResult2() {
        NetworkMetricNakhleh metric = new NetworkMetricNakhleh();
        double dist = metric.computeDistanceBetweenTwoNetworks(_supernetwork, _trueNetwork);
        return Math.abs(dist) < 1e-6;
    }

    public Network getSubNetwork(Network trueNetwork, List<String> selectedLeaves, boolean removeBinaryNodes, boolean alterHeights) {
        int subsize = selectedLeaves.size();
        Network currentNetwork = trueNetwork.clone();
        if(alterHeights) {
            changeHeight(currentNetwork);
        }

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
        return newsubnetwork.clone();
    }

    public void genRandomSubNetworks(Network trueNetwork, int num, int subsize) {
        _subnetworks = new ArrayList<>();

        for(int i = 0 ; i < num ; i++) {
            List<String> leavesToSelect = new ArrayList<>();
            List<String> selectedLeaves = new ArrayList<>();
            for(Object leafNodeObj : trueNetwork.getLeaves()) {
                NetNode leafNode = (NetNode) leafNodeObj;
                leavesToSelect.add(leafNode.getName());
            }
            // pick leaves randomly
            for(int j = 0 ; j < subsize ; j++) {
                int k = Math.abs(_random.nextInt()) % leavesToSelect.size();
                selectedLeaves.add(leavesToSelect.get(k));
                leavesToSelect.remove(k);
            }


            _subnetworks.add(getSubNetwork(trueNetwork, selectedLeaves, true, false));
        }
    }

    public List<Network> genAllSubNetworks(Network trueNetwork, int subsize) {
        _subnetworks = new ArrayList<>();

        List<String> allLeaves = new ArrayList<>();
        for(Object leafObj : trueNetwork.getLeaves()) {
            NetNode leaf = (NetNode) leafObj;
            allLeaves.add(leaf.getName());
        }
        Collections.sort(allLeaves);

        Combinations combinations = new Combinations(trueNetwork.getLeafCount(), subsize);
        Iterator<int[]> combination = combinations.iterator();
        while(combination.hasNext()) {
            int[] indices = combination.next();
            List<String> selectedLeaves = new ArrayList<>();
            for(int i = 0 ; i < indices.length ; i++) {
                selectedLeaves.add(allLeaves.get(indices[i]));
            }
            _subnetworks.add(getSubNetwork(trueNetwork, selectedLeaves, true, true));
        }

        return _subnetworks;

    }

    public double[] getLowerAndUpperBoundOfHeight(NetNode node, Map<NetNode, Double> heights) {
        double[] bounds = new double[] {Double.MIN_VALUE, Double.MAX_VALUE};
        for(Object childObj: node.getChildren()) {
            NetNode child = (NetNode) childObj;
            bounds[0] = Math.max(bounds[0], heights.get(child));
        }
        for(Object parentObj: node.getParents()) {
            NetNode parent = (NetNode) parentObj;
            bounds[1] = Math.min(bounds[1], heights.get(parent));
        }
        return bounds;
    }

    public void setNodeHeight(NetNode node, double newHeight, Map<NetNode, Double> heights) {
        heights.put(node, newHeight);
        for(Object childObj: node.getChildren()) {
            NetNode child = (NetNode) childObj;
            child.setParentDistance(node, newHeight - heights.get(child));
        }
        for(Object parentObj: node.getParents()) {
            NetNode parent = (NetNode) parentObj;
            node.setParentDistance(parent, heights.get(parent) - newHeight);
        }
    }

    public Map<NetNode, Double> getNodeHeights(Network network) {
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
                    if(Math.abs(heightMap.get(parent) - parentHeight) > eps) {
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

    public void changeHeight(Network network) {
        Map<NetNode, Double> heights = getNodeHeights(network);

        for(Object nodeObj : network.bfs()) {
            NetNode node = (NetNode) nodeObj;
            if(node.isLeaf()) {
                continue;
            }
            double[] bounds = getLowerAndUpperBoundOfHeight(node, heights);
            double oldHeight = heights.get(node);
            double bound = oldHeight * 0.1;
            bound = Math.min(bound, oldHeight - bounds[0]);
            bound = Math.min(bound, bounds[1] - oldHeight);
            double newHeight = oldHeight + (_random.nextDouble() - 0.5) * 2 * bound;
            setNodeHeight(node, newHeight, heights);
        }
    }

    protected void adopt(NetNode par, NetNode child, double[] params, Map<NetNode, Double> height) {
        par.adoptChild(child, height.get(par) - height.get(child));
        child.setParentProbability(par, params[0]);
        child.setParentSupport(par, params[1]);
    }

    protected double[] getParameters(NetNode par, NetNode child) {
        return new double[] {child.getParentProbability(par), child.getParentSupport(par)};
    }

    public void addReticulation(Network network, Map<NetNode, Double> height) {
        List<Tuple<NetNode, NetNode>> edges = Networks.getAllEdges(network);
        int numEdges = edges.size();
        int numRetiNodes = network.getReticulationCount();

        boolean success = false;

        while(!success) {
            Tuple<NetNode, NetNode> edge1, edge2;
            edge1 = edges.get(_random.nextInt(numEdges));
            do {
                edge2 = edges.get(_random.nextInt(numEdges));
            } while (edge1 == edge2);

            NetNode v1, v2, v3, v4, v5, v6;
            v3 = edge1.Item2;
            v4 = edge1.Item1;
            v5 = edge2.Item2;
            v6 = edge2.Item1;

            double t3 = height.get(v3);
            double t4 = height.get(v4);
            double l1 = t3 - t4;
            double t1 = t4 + 0.25 * l1 + 0.5 * l1 * _random.nextDouble();

            double t5 = height.get(v5);
            double t6 = height.get(v6);
            double l2 = t5 - t6;
            double t2 = t6 + 0.25 * l2 + 0.25 * l2 * _random.nextDouble();

            boolean validT = true;

            for(NetNode node : height.keySet()) {
                if(Math.abs(t1 - height.get(node)) < 0.005) validT = false;
                if(Math.abs(t2 - height.get(node)) < 0.005) validT = false;
            }

            if(!validT) {
                continue;
            }

            if(t1 < t2 + 0.0001) {
                continue;
            }

            v1 = new BniNetNode<>();
            v2 = new BniNetNode<>();

            height.put(v1, t1);
            height.put(v2, t2);

            double[] paramV3V4 = getParameters(v3, v4);
            double[] paramV5V6 = getParameters(v5, v6);

            double gamma = _random.nextDouble() * 0.5 + 0.25;

            v3.removeChild(v4);
            v5.removeChild(v6);

            adopt(v3, v1, new double[] {NetNode.NO_PROBABILITY, 0.036} , height);
            adopt(v1, v4, paramV3V4, height);
            adopt(v5, v2, new double[] {1.0-gamma, 0.036}, height );
            adopt(v2, v6, paramV5V6, height);

            adopt(v1, v2, new double[] {gamma, 0.036}, height );

            success = true;
        }
    }

    public Network genRandomNetwork(int numTaxa, int numReti) {
        String taxaNames[] = new String[numTaxa];
        for(int i = 0 ; i < numTaxa ; i++) {
            taxaNames[i] = Character.toString((char)('A' + i));
        }
        Tree tree = Trees.generateRandomTree(taxaNames);
        Network network = Networks.readNetwork(tree.toNewick());
        Map<NetNode, Double> height = new HashMap<>();
        double currentHeight = 0;
        // assign heights
        for(Object nodeObj : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeObj;
            if(node.isLeaf()) {
                height.put(node, 0.0);
            } else {
                currentHeight += 1.0 + _random.nextDouble() * 3.0;
                height.put(node, currentHeight);
                for(Object childObj : node.getChildren()) {
                    NetNode child = (NetNode) childObj;
                    child.setParentDistance(node, height.get(node) - height.get(child));
                    child.setParentSupport(node, 0.036);
                }
            }
        }

        network.getRoot().setRootPopSize(0.036);

        if(numReti == 0) {
            return network;
        }

        for(int i = 0 ; i < numReti ; i++) {
            addReticulation(network, height);
            if(Networks.hasCycle(network)) {
                System.err.println("has cycle!");
                System.exit(1);
            }
            if(!Networks.isDisconnectedNetwork(network, null)) {
                System.err.println("disconnected!");
                System.exit(1);
            }
        }

        return network;
    }

    public Network addRandomReticulations(Network net, int numReti) {
        Network network = net.clone();
        Map<NetNode, Double> height = getNodeHeights(network);

        if(numReti == 0) {
            return network;
        }

        for(int i = 0 ; i < numReti ; i++) {
            addReticulation(network, height);
            if(Networks.hasCycle(network)) {
                System.err.println("has cycle!");
                System.exit(1);
            }
            if(!Networks.isDisconnectedNetwork(network, null)) {
                System.err.println("disconnected!");
                System.exit(1);
            }
        }

        return network;
    }

    public static void main(String[] args) {
        final SuperNetwork superNetwork0 = new SuperNetwork(new ArrayList<>());
        for(int i = 0 ; i < 1000 ; ) {
            Network network = Networks.readNetworkWithRootPop("[0.01](A:0.19138876485031905:0.01,(B:0.12682431493516344:0.01,(C:0.09614402921753337:0.01,(((D:0.06864873769305778:0.01,((E:0.04784096272944867:0.01,(F:0.03053172961628988:0.01,G:0.03053172961628988:0.01):0.017309233113158788:0.01):0.007107765870675464:0.01,(H:0.04381981472730136:0.01,I:0.04381981472730136:0.01):0.011128913872822777:0.01):0.013700009092933646:0.01):0.003126652986989553:0.01,(J:0.06689143414599366:0.01,K:0.06689143414599366:0.01):0.004883956534053671:0.01):0.0016928336785998616:0.01,((L:0.044882393305892485:0.01,M:0.044882393305892485:0.01):0.022226283989763562:0.01,(N:0.06296395354420453:0.01,(O:0.03361286516308759:0.01,P:0.03361286516308759:0.01):0.02935108838111694:0.01):0.004144723751451515:0.01):0.006359547062991147:0.01):0.022675804858886178:0.01):0.030680285717630068:0.01):0.06456444991515561:0.01);");
            network = superNetwork0.addRandomReticulations(network, 3);
            if(CoalescenceHistoriesCounting.getNumTaxaUnderReticulation(network) != 5) continue;
            System.out.println(CoalescenceHistoriesCounting.getNumTaxaUnderReticulation(network));
            System.out.println(Networks.getFullString(network));
            System.out.println(Networks.getDendroscopeCompatibleString(network));
        }

        for(int i = 0 ; i < 20000 ; ) {
            SuperNetwork superNetwork = new SuperNetwork(new ArrayList<>());
            Network trueNetwork = superNetwork.genRandomNetwork(20, superNetwork._random.nextInt(20));//Networks.readNetwork("(((((((C:1.0,D:1.0)I1:1.0)I2#H1:2.0,B:4.0)I4:3.0,(I2#H1:3.0,E:5.0)I5:2.0)I7:1.0)I8#H2:1.0,A:9.0)I9:2.0,(((F:3.0,G:3.0)I3:3.0,H:6.0)I6:4.0,I8#H2:2.0)I10:1.0)I11;");
            //trueNetwork = Networks.readNetwork("(((((((C:1.0,D:1.0)I1:1.0)I2#H1:2.0,B:4.0)I4:3.0,(I2#H1:3.0,E:5.0)I5:2.0)I7:1.0)I8#H2:1.0,A:9.0)I9:2.0,(((F:3.0,G:3.0)I3:3.0,H:6.0)I6:4.0,I8#H2:2.0)I10:1.0)I11;");
            System.out.println("True network: " + trueNetwork.toString());

            superNetwork._subnetworks = new ArrayList<>();
            superNetwork._subnetworks.add(trueNetwork);
            superNetwork.findNodeLabels();
            List<Network> subNetworks = superNetwork.getSubNetworks();
            System.out.println("Converted: " + subNetworks.get(0));
            superNetwork._trueNetwork = subNetworks.get(0);

            superNetwork.genAllSubNetworks(trueNetwork, 3);
            //List<Network> subNetworks = superNetwork.getSubNetworks();
            for (Network net : superNetwork._subnetworks) {
                System.out.println(net.toString());
            }

            superNetwork.findNodeLabels();
            superNetwork.findConnections();
            if(!superNetwork.buildNetwork()) {
                continue;
            }
            boolean correctness = superNetwork.checkResult2();
            System.out.println("Correctness: " + correctness);
            if(!correctness) {
                break;
            }
        }
    }
}
