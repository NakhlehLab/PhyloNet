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

    public double computeMeanDistance() {
        if (!Networks.leafSetsAgree(this.network1, this.network2)) {
            throw new RuntimeException("Networks must have identical leaf sets");
        }

        double[][] matrix1 = computeMeanDistanceMatrix(this.network1);
        double[][] matrix2 = computeMeanDistanceMatrix(this.network2);
        double[][] differenceMatrix = matrixDifference(matrix1, matrix2);

        return frobeniusNorm(differenceMatrix);
    }

    /**
     * Given network, we compute mean distance matrix, which is the average of all distances in between leaves i and j.
     *
     * @param network Network to get list of taxa from.
     * @param <T> Indicates the type of additional data this node will store.
     * @return Mean distance matrix.
     */
    private static <T> double[][] computeMeanDistanceMatrix(Network<T> network) {
        // Get distances in between each pair of leaves.
        HashMap<NetNode<T>, HashMap<NetNode<T>, List<Double>>> leafDistances = computeLeafDistances(network);

        // Get taxa to index mapping to build the matrix.
        Map<String, Integer> taxaToIndexMapping = computeTaxaToIndexMapping(network);

        // Get number of taxa.
        int numTaxa = taxaToIndexMapping.keySet().size();

        // Initialize our matrix.
        double[][] meanDistanceMatrix = new double[numTaxa][numTaxa];

        for (Map.Entry<NetNode<T>, HashMap<NetNode<T>, List<Double>>> entry : leafDistances.entrySet()) {
            NetNode<T> leaf1 = entry.getKey();

            for (Map.Entry<NetNode<T>, List<Double>> entry2 : entry.getValue().entrySet()) {
                NetNode<T> leaf2 = entry2.getKey();

                // Compute average distance by taking mean of all distances and place it in the matrix.
                meanDistanceMatrix[taxaToIndexMapping.get(leaf1.getName())]
                        [taxaToIndexMapping.get(leaf2.getName())] =
                        entry2.getValue().stream().mapToDouble(i -> i).average().orElse(0);
            }
        }

        return meanDistanceMatrix;
    }

    /**
     * Compute and return taxa name to their index in the matrix mapping.
     *
     * @param network Network to get list of taxa from.
     * @param <T> Indicates the type of additional data this node will store.
     * @return Taxa name to their index in the matrix mapping.
     */
    private static <T> Map<String, Integer> computeTaxaToIndexMapping(Network<T> network) {
        Map<String, Integer> mapping = new HashMap<>();

        List<String> taxa = new ArrayList<>();
        network.getLeaves().forEach(leaf -> taxa.add(leaf.getName()));

        // Sort taxa in alphabetical order.
        Collections.sort(taxa);

        int i = 0;
        for (String taxon : taxa) {
            mapping.put(taxon, i);
            i += 1;
        }

        return mapping;
    }

    /**
     * Builds and returns a distances for all leaves of the given network.
     *
     * @param network Network to build the distances from.
     * @param <T> Indicates the type of additional data this node will store.
     * @return A distance mapping for all leaves of the given network. For each pair of leaves, it encodes all possible
     * distances inside a nested hashmap.
     */
    private static <T> HashMap<NetNode<T>, HashMap<NetNode<T>, List<Double>>> computeLeafDistances(Network<T> network) {
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

                    if (Double.isInfinite(dist) || Double.isNaN(dist))
                        throw new RuntimeException("Branch length cannot be infinite or undefined");

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

                            Set<NetNode<T>> set = new HashSet<>();
                            set.add(childNode1);
                            set.add(childNode2);

                            if (considered.contains(set))
                                continue;

                            considered.add(set);

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

    private static float frobeniusNorm(double[][] matrix)
    {
        // To store the sum of squares of the
        // elements of the given matrix
        int sumSq = 0;
        for (int i = 0; i < matrix.length; i++)
        {
            for (int j = 0; j < matrix[i].length; j++)
            {
                sumSq += (int)Math.pow(matrix[i][j], 2);
            }
        }

        // Return the square root of
        // the sum of squares
        float res = (float)Math.sqrt(sumSq);
        return res;
    }

    private static double[][] matrixDifference(double[][] matrix1, double[][] matrix2) {
        double[][] newMatrix = new double[matrix1.length][matrix1[0].length];

        if (matrix1.length != matrix2.length)
            throw new RuntimeException("Both matrices should have the same dimension");

        for (int i = 0; i < matrix1.length; i++) {
            if (matrix1[i].length != matrix2[i].length)
                throw new RuntimeException("Both matrices should have the same dimension");

            for (int j = 0; j < matrix1[i].length; j++) {
                newMatrix[i][j] = matrix1[i][j] - matrix2[i][j];
            }
        }

        return newMatrix;
    }
}
