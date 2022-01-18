package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

public class MatrixBasedDissimilarity<T> {
    private Network<T> network1;
    private Network<T> network2;

    /**
     * Constructor for Matrix-Based Dissimilarity metric.
     * @param network1 First network.
     * @param network2 Second network.
     */
    public MatrixBasedDissimilarity(Network<T> network1, Network<T> network2){
        this.network1 = network1;
        this.network2 = network2;
    }

    public double compute() {
        if (!Networks.leafSetsAgree(this.network1, this.network2)) {
            throw new RuntimeException("Networks must have identical leaf sets");
        }

        computeDistanceMatrix(this.network1);
        computeDistanceMatrix(this.network2);

        return Double.NaN;
    }

    /**
     * Builds and returns a distance matrix for all leaves of the given network.
     *
     * @param network Network to build the distance matrix from.
     * @param <T> Indicates the type of additional data this node will store.
     * @return A distance matrix for all leaves of the given network.
     */
    private static <T> HashMap<NetNode<T>, HashMap<NetNode<T>, List<Double>>> computeDistanceMatrix(Network<T> network) {
        // Create a distance matrix and set diagonals to 0.
        HashMap<NetNode<T>, HashMap<NetNode<T>, List<Double>>> distanceMatrix = new HashMap<>();
        for (NetNode<T> leaf1 : network.getLeaves()) {
            distanceMatrix.put(leaf1, new HashMap<>());

            for (NetNode<T> leaf2 : network.getLeaves()) {
                distanceMatrix.get(leaf1).put(leaf2, new ArrayList<>());

                if (leaf1.equals(leaf2)) {
                    distanceMatrix.get(leaf1).get(leaf2).add(0.0);
                    distanceMatrix.get(leaf2).get(leaf1).add(0.0);
                }
            }
        }

        // Create a hashmap to store how many branches are explored that go out of the node. This acts as a semaphore,
        // so that we won't look at the node unless all of its branches are already explored.
        HashMap<NetNode<T>, Integer> numBranchesExplored = new HashMap<>();
        Networks.getInternalNodes(network).forEach(i -> numBranchesExplored.put(i, 0));

        // Create a very nested hashmap:
        // Current Node:
        //   Leaf Node:
        //     List<Distance from Current Node to Leaf Node>
        HashMap<NetNode<T>, HashMap<NetNode<T>, List<Double>>> distances = new HashMap<>();

        // Create a queue and add all leaf nodes to the queue.
        List<NetNode<T>> queue = new ArrayList<>();
        for (NetNode<T> node : network.getLeaves()) {
            queue.add(node);
        }

        while (queue.size() > 0) {
            NetNode<T> curNode = queue.remove(0);

            distances.put(curNode, new HashMap<>());

            if (curNode.isLeaf()) {
                // If it is a leaf node, then add itself as with distance of 0.0.
                distances.get(curNode).putIfAbsent(curNode, new ArrayList<>());
                distances.get(curNode).get(curNode).add(0.0);
            } else {
                // For any internal node, obtain distances from its children, and add distance between current node and
                // the distances retrieved from children. If there are multiple paths leading to a leaf node from
                // multiple branches, then we add a new entry to the list mapped to the leaf node.
                for (NetNode<T> childNode : curNode.getChildren()) {
                    double dist = childNode.getParentDistance(curNode);
                    for (Map.Entry<NetNode<T>, List<Double>> entry : distances.get(childNode).entrySet()) {
                        NetNode<T> leafNode = entry.getKey();
                        distances.get(curNode).putIfAbsent(leafNode, new ArrayList<>());

                        for (Double oldDist : entry.getValue()) {
                            distances.get(curNode).get(leafNode).add(oldDist + dist);
                        }
                    }
                }

                if (curNode.getChildCount() > 1) {
                    // For any network node, where the there are more than 2 branches, we look at if this node connects
                    // two leaves through a new path, and if so, we add the distance between two leaves in our distance
                    // matrix.
                    List<Set<NetNode<T>>> considered = new ArrayList<>();

                    for (NetNode<T> childNode1 : curNode.getChildren()) {
                        for (NetNode<T> childNode2 : curNode.getChildren()) {
                            if (childNode1.equals(childNode2))
                                continue;

                            if (considered.contains(Set.of(childNode1, childNode2)))
                                continue;

                            considered.add(Set.of(childNode1, childNode2));

                            double dist1 = childNode1.getParentDistance(curNode);
                            double dist2 = childNode2.getParentDistance(curNode);

                            for (Map.Entry<NetNode<T>, List<Double>> entry1 : distances.get(childNode1).entrySet()) {
                                for (Map.Entry<NetNode<T>, List<Double>> entry2 : distances.get(childNode2).entrySet()) {
                                    NetNode<T> leafNode1 = entry1.getKey();
                                    NetNode<T> leafNode2 = entry2.getKey();

                                    if (leafNode1.equals(leafNode2))
                                        continue;

                                    for (double oldDist1 : entry1.getValue()) {
                                        for (double oldDist2 : entry2.getValue()) {
                                            distanceMatrix.get(leafNode1).get(leafNode2).add(dist1 + dist2 + oldDist1 + oldDist2);
                                            distanceMatrix.get(leafNode2).get(leafNode1).add(dist1 + dist2 + oldDist1 + oldDist2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // If the current node's parent(s) have all of their branches explored, add that parent to queue.
            curNode.getParents().forEach(parent -> {
                numBranchesExplored.put(parent, numBranchesExplored.get(parent) + 1);
                if (numBranchesExplored.get(parent) == parent.getChildCount()) {
                    queue.add(parent);
                }
            });
        }

        // DEBUG PRINT FOR DISTANCES VARIABLE
//        for (Map.Entry<NetNode<T>, HashMap<NetNode<T>, List<Double>>> entry1 : distanceMatrix.entrySet()) {
//            for (Map.Entry<NetNode<T>, List<Double>> entry2 : entry1.getValue().entrySet()) {
//                System.out.print(entry1.getKey().getName());
//                System.out.print(" -> ");
//                System.out.print(entry2.getKey().getName());
//                System.out.print(" : ");
//                System.out.println(entry2.getValue());
//            }
//        }

        // DEBUG PRINT FOR DISTANCES MATRIX
//        for (Map.Entry<NetNode<T>, HashMap<NetNode<T>, List<Double>>> entry1 : distanceMatrix.entrySet()) {
//            for (Map.Entry<NetNode<T>, List<Double>> entry2 : entry1.getValue().entrySet()) {
//                System.out.print(entry1.getKey().getName());
//                System.out.print(" -> ");
//                System.out.print(entry2.getKey().getName());
//                System.out.print(" : ");
//                System.out.println(entry2.getValue());
//            }
//        }

        return distanceMatrix;
    }
}
