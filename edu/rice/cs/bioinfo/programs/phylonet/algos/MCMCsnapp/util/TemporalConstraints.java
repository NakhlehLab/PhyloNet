package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util;

import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;

import java.util.*;

/**
 * Temporal constraints between gene trees and network
 * Created by wendingqiao on 2/25/16.
 */
public class TemporalConstraints {

    public static Map<String, Double> getTemporalConstraints(UltrametricTree gt, Map<String, String> alleles2species) {
        Map<String, Double> constraints = new HashMap<>();
        Map<TNode, Set<String>> clusters = new HashMap<>();
        for(TNode node : gt.getTree().postTraverse()) {
            Set<String> cluster = new HashSet<>();
            if(node.isLeaf()) {
                if(alleles2species == null) {
                    cluster.add(node.getName());
                } else {
                    cluster.add(alleles2species.get(node.getName()));
                }
            } else {
                List<Set<String>> childClusters = new ArrayList<>();
                for(TNode child : node.getChildren()) {
                    childClusters.add(clusters.get(child));
                    cluster.addAll(clusters.get(child));
                }
                if(childClusters.size() != 2) {
                    throw new IllegalArgumentException("TreeNode " + node.getName() + " has " + childClusters.size() + " children.");
                }
                for(String taxa1 : childClusters.get(0)) {
                    for(String taxa2 : childClusters.get(1)) {
                        String key = getTaxaPair(taxa1, taxa2);
                        if(key == null) continue;
                        if(constraints.containsKey(key)) {
                            constraints.put(key, Math.min(constraints.get(key), gt.getNodeHeight(node)));
                        } else {
                            constraints.put(key, gt.getNodeHeight(node));
                        }
                    }
                }
            }
            clusters.put(node, cluster);
        }
        return constraints;
    }

    public static Map<String, Double> getTemporalConstraints(List<UltrametricTree> gts,
                                                             Map<String, List<String>> species2alleles,
                                                             Map<String, String> alleles2species) {
        if(gts == null) {
            return null;
        }
        if(gts.size() == 0) {
            throw new IllegalArgumentException("List of gene trees cannot be null");
        }
        Map<String, Double> constraints = new HashMap<>();
        Set<String> keys = (species2alleles == null) ? getTaxaPairs(gts.get(0).getTree().getLeaves()) :
                getTaxaPairs(species2alleles.keySet());
        for(String key : keys) {
            constraints.put(key, Double.MAX_VALUE);
        }
        for(int i = 0; i < gts.size(); i++) {
            Map<String, Double> temp = getTemporalConstraints(gts.get(i), alleles2species);
            for(String key : temp.keySet()) {
                if(!constraints.keySet().contains(key)) {
                    throw new IllegalArgumentException("Cannot locate key " + key + " in constraint map " + (alleles2species == null));
                }
                constraints.put(key, Math.min(constraints.get(key), temp.get(key)));
            }
        }
        return constraints;
    }

    public static Map<String, Double> getLowerBounds(Network<NetNodeInfo> net) {
        Map<String, Double> bounds = new HashMap<>();
        Map<NetNode<NetNodeInfo>, Set<String>> clusters = new HashMap<>();
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(net)) {
            Set<String> cluster = new HashSet<>();
            if(node.isLeaf()) {
                cluster.add(node.getName());
            } else if(node.isNetworkNode()) {
                NetNode<NetNodeInfo> child = node.getChildren().iterator().next();
                if(node.getChildCount() != 1 || !clusters.containsKey(child)) {
                    throw new RuntimeException("Wrong network/ post traversal: " + net.toString());
                }
                cluster.addAll(clusters.get(child));
            } else {
                if(node.getChildCount() != 2) {
                    throw new RuntimeException("Wrong network: " + net.toString());
                }
                List<NetNode<NetNodeInfo>> children = IterableHelp.toList(node.getChildren());
                Set<String> cluster0 = clusters.get(children.get(0));
                Set<String> cluster1 = clusters.get(children.get(1));
                for(String taxa1 : cluster0) {
                    for(String taxa2 : cluster1) {
                        String key = getTaxaPair(taxa1, taxa2);
                        if(key == null) continue;
                        if(bounds.containsKey(key)) {
                            bounds.put(key, Math.min(bounds.get(key), node.getData().getHeight()));
                        } else {
                            bounds.put(key, node.getData().getHeight());
                        }
                    }
                }
                cluster.addAll(cluster0);
                cluster.addAll(cluster1);
            }
            clusters.put(node, cluster);
        }
        return bounds;
    }

    public static <T> Map<String, NetNode<T>> getRestrictedNodes(Network<T> net) {
        Map<String, NetNode<T>> restrictedNodes = new HashMap<>();
        Map<NetNode<T>, Set<String>> clusters = new HashMap<>();
        for(NetNode<T> node : Networks.postTraversal(net)) {
            Set<String> cluster = new HashSet<>();
            if(node.isLeaf()) {
                cluster.add(node.getName());
            } else if(node.isNetworkNode()) {
                NetNode<T> child = node.getChildren().iterator().next();
                if(node.getChildCount() != 1 || !clusters.containsKey(child)) {
                    throw new RuntimeException("Wrong network/ post traversal: " + net.toString());
                }
                cluster.addAll(clusters.get(child));
            } else {
                if(node.getChildCount() != 2) {
                    throw new RuntimeException("Wrong network: " + net.toString());
                }
                List<NetNode<T>> children = IterableHelp.toList(node.getChildren());
                Set<String> cluster0 = clusters.get(children.get(0));
                Set<String> cluster1 = clusters.get(children.get(1));
                for(String taxa1 : cluster0) {
                    for(String taxa2 : cluster1) {
                        String key = getTaxaPair(taxa1, taxa2);
                        if(key == null || restrictedNodes.containsKey(key)) continue;
                        restrictedNodes.put(key, node);
                    }
                }
                cluster.addAll(cluster0);
                cluster.addAll(cluster1);
            }
            clusters.put(node, cluster);
        }
        return restrictedNodes;
    }

    public static <T> Map<NetNode<T>, Set<String>> getNodeRestriction(Network<T> net) {
        Map<NetNode<T>, Set<String>> restrictedNodes = new HashMap<>();
        Map<NetNode<T>, Set<String>> clusters = new HashMap<>();
        Set<String> usedKey = new HashSet<>();
        for(NetNode<T> node : Networks.postTraversal(net)) {
            Set<String> cluster = new HashSet<>();
            if(node.isLeaf()) {
                cluster.add(node.getName());
            } else if(node.isNetworkNode()) {
                NetNode<T> child = node.getChildren().iterator().next();
                if(node.getChildCount() != 1 || !clusters.containsKey(child)) {
                    throw new RuntimeException("Wrong network/ post traversal: " + net.toString());
                }
                cluster.addAll(clusters.get(child));
            } else {
                if(node.getChildCount() != 2) {
                    throw new RuntimeException("Wrong network: " + net.toString());
                }
                List<NetNode<T>> children = IterableHelp.toList(node.getChildren());
                Set<String> cluster0 = clusters.get(children.get(0));
                Set<String> cluster1 = clusters.get(children.get(1));
                for(String taxa1 : cluster0) {
                    for(String taxa2 : cluster1) {
                        String key = getTaxaPair(taxa1, taxa2);
                        if(key == null || usedKey.contains(key)) continue;
                        if(!restrictedNodes.containsKey(node)) {
                            restrictedNodes.put(node, new HashSet<String>());
                        }
                        restrictedNodes.get(node).add(key);
                        usedKey.add(key);
                    }
                }
                cluster.addAll(cluster0);
                cluster.addAll(cluster1);
            }
            clusters.put(node, cluster);
        }
        return restrictedNodes;
    }

    public static <T> Set<String> getTaxaPairs(Set<String> taxon) {
        Set<String> pairs = new HashSet<>();
        List<String> keys = new ArrayList<>(taxon);
        Collections.sort(keys);
        for(int i = 0; i < keys.size(); i++) {
            for(int j = i+1; j < keys.size(); j++) {
                pairs.add(keys.get(i) + ";" + keys.get(j));
            }
        }
        return pairs;
    }

    public static <T> Set<String> getTaxaPairs(String[] taxon) {
        Set<String> pairs = new HashSet<>();
        Arrays.sort(taxon);
        for(int i = 0; i < taxon.length; i++) {
            for(int j = i+1; j < taxon.length; j++) {
                pairs.add(taxon[i] + ";" + taxon[j]);
            }
        }
        return pairs;
    }

    public static String getTaxaPair(String taxa1, String taxa2) {
        if(taxa1.compareTo(taxa2) == 0) return null;
        return taxa1.compareTo(taxa2) < 0 ? taxa1+";"+taxa2 : taxa2+";"+taxa1;
    }

}
