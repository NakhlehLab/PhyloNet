package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;
import java.util.stream.Collectors;

public class AveragePathDistance {
    /**
     * Compute APD.
     *
     * @param <T> Type of the network.
     * @param network1 First network.
     * @param network2 Second network.
     * @return APD of two networks.
     */
    public static <T> double compute(Network<T> network1, Network<T> network2) {
        if (!Networks.leafSetsAgree(network1, network2)) {
            throw new RuntimeException("Networks must have identical leaf sets");
        }

        double[][] matrix1 = computeMatrix(network1);
        double[][] matrix2 = computeMatrix(network2);
        double[][] differenceMatrix = matrixDifference(matrix1, matrix2);

        return frobeniusNorm(differenceMatrix);
    }

    /**
     * Compute APD, where matrices are pre-computed.
     *
     * @param matrix1 M1 representing first network.
     * @param matrix2 M2 representing second network.
     * @return APD of two networks.
     */
    public static double compute(double[][] matrix1, double[][] matrix2) {
        double[][] differenceMatrix = matrixDifference(matrix1, matrix2);

        return frobeniusNorm(differenceMatrix);
    }

    /**
     * Compute NormAPD.
     *
     * @param <T> Type of the network.
     * @param network1 First network.
     * @param network2 Second network.
     * @return NormAPD of two networks.
     */
    public static <T> double computeNormalized(Network<T> network1, Network<T> network2) {
        if (!Networks.leafSetsAgree(network1, network2)) {
            throw new RuntimeException("Networks must have identical leaf sets");
        }

        double[][] matrix1 = computeMatrix(network1);
        double[][] matrix2 = computeMatrix(network2);
        double[][] differenceMatrix = matrixDifference(matrix1, matrix2);

        return frobeniusNorm(differenceMatrix) / Math.min(frobeniusNorm(matrix1), frobeniusNorm(matrix2));
    }

    /**
     * Compute NormAPD, where matrices are pre-computed.
     *
     * @param matrix1 M1 representing first network.
     * @param matrix2 M2 representing second network.
     * @return NormAPD of two networks.
     */
    public static double computeNormalized(double[][] matrix1, double[][] matrix2) {
        double[][] differenceMatrix = matrixDifference(matrix1, matrix2);

        return frobeniusNorm(differenceMatrix) / Math.min(frobeniusNorm(matrix1), frobeniusNorm(matrix2));
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
     * Given network, we compute average path distance.
     *
     * @param network Network to get list of taxa from.
     * @param <T> Indicates the type of additional data this node will store.
     * @return Mean distance matrix.
     */
    public static <T> double[][] computeMatrix(Network<T> network) {
        // Create a distance matrix and set diagonals to 0.
        HashMap<NetNode<T>, HashMap<NetNode<T>, List<Double>>> pathDistanceMatrix = new HashMap<>();
        for (NetNode<T> leaf1 : network.getLeaves()) {
            pathDistanceMatrix.put(leaf1, new HashMap<>());

            for (NetNode<T> leaf2 : network.getLeaves()) {
                pathDistanceMatrix.get(leaf1).put(leaf2, new ArrayList<>());

                if (leaf1.equals(leaf2)) {
                    pathDistanceMatrix.get(leaf1).get(leaf2).add(0.0);
                }
            }
        }

        // Create a probability matrix and set diagonals to 1.
        HashMap<NetNode<T>, HashMap<NetNode<T>, List<Double>>> probabilityMatrix = new HashMap<>();
        for (NetNode<T> leaf1 : network.getLeaves()) {
            probabilityMatrix.put(leaf1, new HashMap<>());

            for (NetNode<T> leaf2 : network.getLeaves()) {
                probabilityMatrix.get(leaf1).put(leaf2, new ArrayList<>());

                if (leaf1.equals(leaf2)) {
                    probabilityMatrix.get(leaf1).get(leaf2).add(1.0);
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

        // Create a very nested hashmap:
        // Current Node:
        //   Leaf Node:
        //     List<Probability from Current Node to Leaf Node>
        HashMap<NetNode<T>, HashMap<NetNode<T>, List<Double>>> probabilities = new HashMap<>();

        // Create a very nested hashmap:
        // Current Node:
        //   Leaf Node:
        //     List<List<Path from Current Node to Leaf Node>>
        HashMap<NetNode<T>, HashMap<NetNode<T>, List<List<NetNode<T>>>>> paths = new HashMap<>();

        // Create a queue and add all leaf nodes to the queue.
        List<NetNode<T>> queue = new ArrayList<>();
        for (NetNode<T> node : network.getLeaves()) {
            queue.add(node);
        }

        while (queue.size() > 0) {
            NetNode<T> curNode = queue.remove(0);

            distances.put(curNode, new HashMap<>());
            probabilities.put(curNode, new HashMap<>());
            paths.put(curNode, new HashMap<>());

            if (curNode.isLeaf()) {
                // If it is a leaf node, then add itself as with distance of 0.0, and set probability to 1.0.
                distances.get(curNode).putIfAbsent(curNode, new ArrayList<>());
                distances.get(curNode).get(curNode).add(0.0);

                probabilities.get(curNode).putIfAbsent(curNode, new ArrayList<>());
                probabilities.get(curNode).get(curNode).add(1.0);

                paths.get(curNode).putIfAbsent(curNode, new ArrayList<>());
                paths.get(curNode).get(curNode).add(new ArrayList<>());
                paths.get(curNode).get(curNode).get(0).add(curNode);
            } else {
                // For any internal node, obtain distances from its children, and add distance between current node and
                // the distances retrieved from children. If there are multiple paths leading to a leaf node from
                // multiple branches, then we add a new entry to the list mapped to the leaf node.
                for (NetNode<T> childNode : curNode.getChildren()) {
                    double distFromChildToCur = childNode.getParentDistance(curNode);
                    if (Double.isInfinite(distFromChildToCur) || Double.isNaN(distFromChildToCur))
                        throw new RuntimeException("Branch length cannot be infinite or undefined");

                    double probFromChildToCur = 1.0;

                    for (Map.Entry<NetNode<T>, List<Double>> entry : distances.get(childNode).entrySet()) {
                        NetNode<T> leafNode = entry.getKey();
                        distances.get(curNode).putIfAbsent(leafNode, new ArrayList<>());

                        for (Double distFromChildToLeaf : entry.getValue()) {
                            distances.get(curNode).get(leafNode).add(distFromChildToLeaf + distFromChildToCur);
                        }
                    }

                    for (Map.Entry<NetNode<T>, List<Double>> entry : probabilities.get(childNode).entrySet()) {
                        NetNode<T> leafNode = entry.getKey();
                        probabilities.get(curNode).putIfAbsent(leafNode, new ArrayList<>());

                        for (Double probFromChildToLeaf : entry.getValue()) {
                            probabilities.get(curNode).get(leafNode).add(probFromChildToLeaf * probFromChildToCur);
                        }
                    }

                    for (Map.Entry<NetNode<T>, List<List<NetNode<T>>>> entry : paths.get(childNode).entrySet()) {
                        NetNode<T> leafNode = entry.getKey();
                        paths.get(curNode).putIfAbsent(leafNode, new ArrayList<>());

                        for (List<NetNode<T>> pathFromChildToLeaf : entry.getValue()) {
                            List<NetNode<T>> pathFromCurToLeaf = new ArrayList<>();
                            pathFromCurToLeaf.add(curNode);
                            pathFromCurToLeaf.addAll(pathFromChildToLeaf);
                            paths.get(curNode).get(leafNode).add(pathFromCurToLeaf);
                        }
                    }
                }

                if (curNode.getChildCount() > 1) {
                    HashMap<NetNode<T>, List<Double>> curDistances = distances.get(curNode);
                    HashMap<NetNode<T>, List<Double>> curProbabilities = probabilities.get(curNode);
                    HashMap<NetNode<T>, List<List<NetNode<T>>>> curPaths = paths.get(curNode);

                    if (!curDistances.keySet().containsAll(curProbabilities.keySet()))
                        throw new RuntimeException("Assertion failure. There's a bug in the code.");

                    if (!curProbabilities.keySet().containsAll(curDistances.keySet()))
                        throw new RuntimeException("Assertion failure. There's a bug in the code.");

                    if (!curDistances.keySet().containsAll(curPaths.keySet()))
                        throw new RuntimeException("Assertion failure. There's a bug in the code.");

                    if (!curPaths.keySet().containsAll(curDistances.keySet()))
                        throw new RuntimeException("Assertion failure. There's a bug in the code.");

                    Set<NetNode<T>> accessibleLeaves = curDistances.keySet();

                    for (NetNode<T> leaf1 : accessibleLeaves) {
                        for (NetNode<T> leaf2 : accessibleLeaves) {
                            if (leaf1.equals(leaf2))
                                continue;

                            if (curDistances.get(leaf1).size() != curProbabilities.get(leaf1).size())
                                throw new RuntimeException("Assertion failure. There's a bug in the code.");

                            if (curDistances.get(leaf1).size() != curDistances.get(leaf1).size())
                                throw new RuntimeException("Assertion failure. There's a bug in the code.");

                            if (curDistances.get(leaf2).size() != curProbabilities.get(leaf2).size())
                                throw new RuntimeException("Assertion failure. There's a bug in the code.");

                            if (curDistances.get(leaf2).size() != curDistances.get(leaf2).size())
                                throw new RuntimeException("Assertion failure. There's a bug in the code.");

                            for (int i = 0; i < curDistances.get(leaf1).size(); i++) {
                                for (int j = 0; j < curDistances.get(leaf2).size(); j++) {
                                    Set<NetNode<T>> commonNodes = curPaths.get(leaf1).get(i).stream()
                                    .distinct()
                                    .filter(curPaths.get(leaf2).get(j)::contains)
                                    .collect(Collectors.toSet());

                                    if (commonNodes.size() == 1 && commonNodes.contains(curNode)) {
                                        pathDistanceMatrix.get(leaf1).get(leaf2).add(
                                                curDistances.get(leaf1).get(i) * curProbabilities.get(leaf1).get(i) +
                                                curDistances.get(leaf2).get(j) * curProbabilities.get(leaf2).get(j));
                                        probabilityMatrix.get(leaf1).get(leaf2).add(curProbabilities.get(leaf1).get(i) * curProbabilities.get(leaf2).get(j));
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

        // Get taxa to index mapping to build the matrix.
        Map<String, Integer> taxaToIndexMapping = computeTaxaToIndexMapping(network);

        // Get number of taxa.
        int numTaxa = taxaToIndexMapping.keySet().size();

        // Initialize our matrix.
        double[][] averagePathDistanceMatrix = new double[numTaxa][numTaxa];

        for (Map.Entry<NetNode<T>, HashMap<NetNode<T>, List<Double>>> entry : pathDistanceMatrix.entrySet()) {
            NetNode<T> leaf1 = entry.getKey();

            for (Map.Entry<NetNode<T>, List<Double>> entry2 : entry.getValue().entrySet()) {
                NetNode<T> leaf2 = entry2.getKey();

                // Assert each distance has a probability.
                if (probabilityMatrix.get(leaf1).get(leaf2).size() != pathDistanceMatrix.get(leaf1).get(leaf2).size())
                    throw new RuntimeException("Assertion failure. There's a bug in the code.");

                averagePathDistanceMatrix[taxaToIndexMapping.get(leaf1.getName())]
                        [taxaToIndexMapping.get(leaf2.getName())] = pathDistanceMatrix.get(leaf1).get(leaf2).stream().mapToDouble(i -> i).average().getAsDouble();
            }
        }

        return averagePathDistanceMatrix;
    }

    /**
     * Compute Frobenius norm.
     * @param matrix Matrix.
     * @return Frobenius norm of matrix.
     */
    private static double frobeniusNorm(double[][] matrix) {
        double sumSq = 0;
        for (int i = 0; i < matrix.length; i++)
        {
            for (int j = 0; j < matrix[i].length; j++)
            {
                sumSq += Math.pow(matrix[i][j], 2);
            }
        }

        return Math.sqrt(sumSq);
    }

    /**
     * Compute matrix difference.
     * @param matrix1 First matrix.
     * @param matrix2 Second matrix.
     * @return Difference matrix.
     */
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
