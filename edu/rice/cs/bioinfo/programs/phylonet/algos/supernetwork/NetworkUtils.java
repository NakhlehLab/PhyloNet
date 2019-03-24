package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import org.apache.commons.math3.util.Combinations;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 11/24/18
 * Time: 11:30 AM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkUtils {

    public static Network getSubNetwork(Network trueNetwork, List<String> selectedLeaves, boolean removeBinaryNodes) {
        int subsize = selectedLeaves.size();
        Network currentNetwork = trueNetwork.clone();

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

    public static List<Network> genAllSubNetworks(Network trueNetwork, int subsize) {
        List<Network> _subnetworks = new ArrayList<>();

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
            _subnetworks.add(getSubNetwork(trueNetwork, selectedLeaves, true));
        }

        return _subnetworks;

    }

    public static Map<NetNode, Double> getNodeHeightMapping(Network network) {
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
                    if(Math.abs(heightMap.get(parent) - parentHeight) > 1e-5) {
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

    public static void setNodeHeightMapping(NetNode node, double newHeight, Map<NetNode, Double> heights) {
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

    public static double[] getLowerAndUpperBoundOfHeight(NetNode node, Map<NetNode, Double> heights) {
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

    public static void alterHeights(Network network, Random random) {
        alterHeights(network, random, 0.001);
    }

    public static void alterHeights(Network network, Random random, double range) {
        Map<NetNode, Double> heights = getNodeHeightMapping(network);

        for(Object nodeObj : network.bfs()) {
            NetNode node = (NetNode) nodeObj;
            if(node.isLeaf()) {
                continue;
            }
            double[] bounds = getLowerAndUpperBoundOfHeight(node, heights);
            double oldHeight = heights.get(node);
            double bound = oldHeight * range;
            bound = Math.min(bound, oldHeight - bounds[0]);
            bound = Math.min(bound, bounds[1] - oldHeight);
            double newHeight = oldHeight + (random.nextDouble() - 0.5) * 2 * bound;
            setNodeHeightMapping(node, newHeight, heights);
        }
    }

    public static Set<NetNode> GetAllRetiAboveLeaf(Network network, String leafname) {
        Set<NetNode> set1 = SuperNetwork3.getAllAncestors(network.findNode(leafname));
        Set<NetNode> set2 = new HashSet<>();
        for(NetNode node : set1) {
            if(node.isNetworkNode()) {
                set2.add(node);
            }
        }
        return set2;
    }

    public static int GetNumAllRetiAboveLeaf(Network network, String leafname) {
        return GetAllRetiAboveLeaf(network, leafname).size();
    }

    public static double ComputeScore(Network netToTry, SuperNetwork3 sn) {
        double score = 0.0;
        //List<Network> backbonesToTry = Pipeline.getAllBackboneNets(netToTry);
        List<Network> backbonesToTry = new ArrayList<>();
        backbonesToTry.add(netToTry.clone());
        sn.Prepare();

        for (SuperNetwork3.NetworkWithInfo netinfo : sn.subnetworks_) {
            if(!netinfo.trustTime) {
                continue;
            }

            if (netinfo.dirty) {
                netinfo.network = Networks.readNetwork(netinfo.backup);
                SuperNetwork3.initNetHeights(netinfo.network);
                netinfo.dirty = false;
            }
            List<String> leavesIntersection = new ArrayList<>();
            for(Object leafObj : netToTry.getLeaves()) {
                NetNode leaf = (NetNode) leafObj;
                leavesIntersection.add(leaf.getName());
            }
            leavesIntersection.retainAll(netinfo.taxa);
            if(leavesIntersection.size() < 2) continue;

            double minScore = Double.MAX_VALUE;
            for(Network backbone : backbonesToTry) {
                double curScore = SuperNetwork3.compareTwoSubNetwork(netinfo.network, backbone, leavesIntersection);
                minScore = Math.min(minScore, curScore);
                if(minScore == 0) break;
            }

            if(minScore == 0) minScore = 1;
            else minScore = 0;
            score += minScore;

            //score += minScore * netinfo.percentage;
        }

        //score += netToTry.getReticulationCount() * 8.0;

        return score;
    }
}
