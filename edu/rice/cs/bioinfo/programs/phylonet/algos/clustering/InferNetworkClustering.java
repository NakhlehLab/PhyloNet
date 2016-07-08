package edu.rice.cs.bioinfo.programs.phylonet.algos.clustering;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 7/7/16
 * Time: 2:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class InferNetworkClustering {
    private String _clusteringMethod = "MDS+KMeans";
    private String _inferringMethod = "MDC";
    private Random _random = new Random();
    private int _MDSdim = 3;
    private InferTreeWrapper _inferTreeWrapper = new InferTreeWrapper();

    public void getMDScoordinates(int n, int d, double[][] distMatrix, double [][] coordinates) {
        MultidimensionalScaling multidimensionalScaling = new MultidimensionalScaling();
        multidimensionalScaling.getMDScoordinates(n, d, distMatrix, coordinates);
    }

    private double getRootedRobinsonFouldsDistance(Tree tree1, Tree tree2) {
        SymmetricDifference symmetricDifference = new SymmetricDifference();
        symmetricDifference.computeDifference(tree1, tree2, true);
        double diff = symmetricDifference.getWeightedAverage();

        return diff;
    }

    public double getEuclideanDistance(int dim, double [] x1, double [] x2) {
        double d = 0;
        for(int i = 0 ; i < dim ; i++)
            d += (x1[i] - x2[i]) * (x1[i] - x2[i]);
        d = Math.sqrt(d);
        return d;
    }

    public double getSilhouetteIndex(List<List<Integer>> cluster, double [][] distMatrix, int n, int m) {
        double s = 0;
        for(int k = 0 ; k < m ; k++) {
            for(int i : cluster.get(k)) {
                double ai = 0;
                for(int j : cluster.get(k)) {
                    ai += distMatrix[i][j];
                }
                ai /= cluster.get(k).size();

                double bi = Double.MAX_VALUE;
                for(int k0 = 0 ; k0 < m ; k0++) {
                    if(k0 == k) continue;
                    double b = 0;
                    for(int j : cluster.get(k0)) {
                        b += distMatrix[i][j];
                    }
                    bi = Math.min(bi, b / cluster.get(k0).size());
                }

                double si = 0;
                if(cluster.get(k).size() > 1)
                    si = (bi - ai) / Math.max(ai, bi);
                s += si;
            }
        }
        s /= n;
        return s;
    }

    private void KMeansClustering_MDC(List<List<Integer>> cluster, List<List<MutableTuple<Tree, Double>>> trees, int n, int m) {
        int maxIter = 1000;
        List<Tree> centroids = new ArrayList<>();
        //double centroids[][] = new double[m][dim];
        cluster.clear();
        for(int i = 0 ; i < m ; i++) {
            cluster.add(new ArrayList<Integer>());
            centroids.add(trees.get(i).get(0).Item1);
        }

        int assignment[] = new int[n];
        boolean changed = true;
        int iter = 0;

        while(changed && iter < maxIter) {
            //System.out.println("Iteration: " + iter);
            changed = false;
            iter++;

            for (int i = 0; i < n; i++) {
                int tmpCluster = 0;
                double minD = getRootedRobinsonFouldsDistance(trees.get(i).get(0).Item1, centroids.get(0));
                //double minD = getEuclideanDistance(dim, coordinates[i], centroids[0]);
                for (int j = 1; j < m; j++) {
                    //double D = getEuclideanDistance(dim, coordinates[i], centroids[j]);
                    double D = getRootedRobinsonFouldsDistance(trees.get(i).get(0).Item1, centroids.get(j));
                    if (D < minD) {
                        minD = D;
                        tmpCluster = j;
                    }
                }
                assignment[i] = tmpCluster;
            }

            List<Tree> newcentroids = new ArrayList<>();
            List<List<Tree>> gts = new ArrayList<>();
            //double newcentroids[][] = new double[m][dim];
            int count[] = new int[m];

            for (int i = 0; i < m; i++) {
                count[i] = 0;
                newcentroids.add(null);
                gts.add(new ArrayList<>());
            }


            for (int i = 0; i < n; i++) {
                count[assignment[i]]++;
                gts.get(assignment[i]).add(trees.get(i).get(0).Item1);
            }

            for (int i = 0; i < m; i++) {
                if (count[i] > 0) {
                    newcentroids.set(i, Trees.readTree(_inferTreeWrapper.inferTreeByMethod(gts.get(i), null, "MDC")));

                    if (getRootedRobinsonFouldsDistance(newcentroids.get(i), centroids.get(i)) > 1e-6) {
                        centroids.set(i, newcentroids.get(i));
                        changed = true;
                    }
                } else {
                    changed = true;
                    centroids.set(i, trees.get(_random.nextInt(n)).get(0).Item1);
                }
            }
        }

        for(int i = 0 ; i < n ; i++) {
            cluster.get(assignment[i]).add(i);
        }
        //System.out.println("Iteration: " + iter);

    }

    private void KMeansClustering(List<List<Integer>> cluster, double [][] coordinates, int n, int m, int dim) {
        int maxIter = 1000;
        double centroids[][] = new double[m][dim];
        cluster.clear();
        for(int i = 0 ; i < m ; i++) {
            cluster.add(new ArrayList<Integer>());
            centroids[i] = coordinates[i];
        }

        int assignment[] = new int[n];
        boolean changed = true;
        int iter = 0;

        while(changed && iter < maxIter) {
            //System.out.println("Iteration: " + iter);
            changed = false;
            iter++;

            for (int i = 0; i < n; i++) {
                int tmpCluster = 0;
                double minD = getEuclideanDistance(dim, coordinates[i], centroids[0]);
                for (int j = 1; j < m; j++) {
                    double D = getEuclideanDistance(dim, coordinates[i], centroids[j]);
                    if (D < minD) {
                        minD = D;
                        tmpCluster = j;
                    }
                }
                assignment[i] = tmpCluster;
            }

            double newcentroids[][] = new double[m][dim];
            int count[] = new int[m];

            for (int i = 0; i < m; i++) {
                count[i] = 0;
                for (int j = 0; j < dim; j++)
                    newcentroids[i][j] = 0;
            }


            for (int i = 0; i < n; i++) {
                count[assignment[i]]++;
                for (int j = 0; j < dim; j++)
                    newcentroids[assignment[i]][j] += coordinates[i][j];
            }

            for (int i = 0; i < m; i++) {
                if (count[i] > 0) {
                    for (int j = 0; j < dim; j++) {
                        newcentroids[i][j] /= count[i];
                    }

                    if (getEuclideanDistance(dim, newcentroids[i], centroids[i]) > 1e-6) {
                        centroids[i] = newcentroids[i];
                        changed = true;
                    }
                } else {
                    changed = true;
                    centroids[i] = coordinates[_random.nextInt(n)];
                }
            }
        }

        for(int i = 0 ; i < n ; i++) {
            cluster.get(assignment[i]).add(i);
        }
        //System.out.println("Iteration: " + iter);

    }

    //n points, m clusters
    private void KMedoidsClustering(List<List<Integer>> cluster, double [][] distMatrix, int n, int m) {
        int maxIter = 1000;
        int medoids[] = new int[m];
        cluster.clear();
        for(int i = 0 ; i < m ; i++) {
            medoids[i] = _random.nextInt(n);
            cluster.add(new ArrayList<Integer>());
        }
        int assignment[] = new int[n];

        boolean changed = true;
        int count = 0;
        while(changed && count < maxIter) {
            System.out.println("Iteration: " + count);
            changed = false;
            count++;
            for(int i = 0 ; i < n ; i++) {
                double bestDist = distMatrix[i][medoids[0]];
                int bestIndex = 0;
                for(int j = 1 ; j < m ; j++) {
                    double tmpDist = distMatrix[i][medoids[j]];
                    if(tmpDist < bestDist) {
                        bestDist = tmpDist;
                        bestIndex = j;
                    }
                }
                if(assignment[i] != bestIndex) {
                    assignment[i] = bestIndex;
                    changed = true;
                }
            }

            for(int k = 0 ; k < m ; k++) {
                int medoid = medoids[k];
                for(int i = 0 ; i < n ; i++) {
                    int bestMedoid = medoid;
                    double lowestCostDelta = 0;
                    if(i != medoid && assignment[i] == k) {
                        double costDelta = 0;
                        for(int j = 0 ; j < n ; j++) {
                            if(assignment[j] == k) {
                                double oldDist = distMatrix[medoid][j];
                                double newDist = distMatrix[i][j];
                                costDelta += newDist - oldDist;
                            }
                        }
                        if(costDelta < lowestCostDelta) {
                            bestMedoid = i;
                            lowestCostDelta = costDelta;
                        }

                        if(bestMedoid != medoid) {
                            medoids[k] = bestMedoid;
                            changed = true;
                        }
                    }
                }
            }
        }

        for(int i = 0 ; i < n ; i++) {
            cluster.get(assignment[i]).add(i);
        }
    }


    public String inferNetwork(List<List<MutableTuple<Tree, Double>>> gts) {

        int size = gts.size();

        double distMatrix[][] = new double[size][size];
        for (int i = 0; i < size; i++)
            for (int j = i + 1; j < size; j++)
                distMatrix[i][j] = distMatrix[j][i] = getRootedRobinsonFouldsDistance(gts.get(i).get(0).Item1, gts.get(j).get(0).Item1);

        int dim = _MDSdim;
        double coordinates[][] = new double[size][dim];
        getMDScoordinates(size, dim, distMatrix, coordinates);


        double maxSilhouetteIndex = Double.MIN_VALUE;
        int bestNumOfClusters = 2;
        int numOfClusters = 2;
        //determine number of clusters
        while (numOfClusters < 100) {
            List<List<Integer>> cluster = new ArrayList<>();


            if(_clusteringMethod.equals("MDS+KMeans")) {
                KMeansClustering(cluster, coordinates, size, numOfClusters, dim);
            }
            else if(_clusteringMethod.equals("MDC+KMeans")) {
                KMeansClustering_MDC(cluster, gts, size, numOfClusters);
            }
            else if(_clusteringMethod.equals("MDS+KMedoids")) {
                double newdistMatrix[][] = new double[size][size];
                for (int i = 0; i < size; i++)
                    for (int j = i + 1; j < size; j++)
                        newdistMatrix[i][j] = newdistMatrix[j][i] = getEuclideanDistance(dim, coordinates[i], coordinates[j]);

                KMedoidsClustering(cluster, newdistMatrix, size, numOfClusters);
            }
            else
                throw new RuntimeException("Unknown clustering method: " + _clusteringMethod);

            double score = getSilhouetteIndex(cluster, distMatrix, size, numOfClusters);
            //System.out.println("numOfCluster = " + numOfClusters);
            //System.out.println("SilhouetteIndex = " + score);

            if (maxSilhouetteIndex < score) {
                maxSilhouetteIndex = score;
                bestNumOfClusters = numOfClusters;
            } else {
                bestNumOfClusters = numOfClusters - 1;
                break;
            }

            numOfClusters++;


        }

        numOfClusters = bestNumOfClusters;
        List<List<Integer>> cluster = new ArrayList<>();

        List<List<MutableTuple<Tree, Double>>> inferredTrees = new ArrayList<>();

        if(_clusteringMethod.equals("MDS+KMeans")) {
            KMeansClustering(cluster, coordinates, size, numOfClusters, dim);
        }
        else if(_clusteringMethod.equals("MDC+KMeans")) {
            KMeansClustering_MDC(cluster, gts, size, numOfClusters);
        }
        else if(_clusteringMethod.equals("MDS+KMedoids")) {
            double newdistMatrix[][] = new double[size][size];
            for (int i = 0; i < size; i++)
                for (int j = i + 1; j < size; j++)
                    newdistMatrix[i][j] = newdistMatrix[j][i] = getEuclideanDistance(dim, coordinates[i], coordinates[j]);

            KMedoidsClustering(cluster, newdistMatrix, size, numOfClusters);
        }
        else
            throw new RuntimeException("Unknown clustering method: " + _clusteringMethod);

        for (int k = 0; k < numOfClusters; k++) {
            List<Tree> currentGTs = new ArrayList<>();
            for (Integer i : cluster.get(k))
                currentGTs.add(gts.get(i).get(0).Item1);
            System.out.println("Cluster size = " + currentGTs.size());
            if (currentGTs.size() == 0)
                continue;
            String speciesTree = _inferTreeWrapper.inferTreeByMethod(currentGTs, null, new String(_inferringMethod));
            System.out.println(speciesTree);
            inferredTrees.add(Arrays.asList(new MutableTuple<Tree, Double>(Trees.readTree(speciesTree), new Double(currentGTs.size()))));

        }

        InferNetworkFromParentalTrees inferNetworkFromParentalTrees = new InferNetworkFromParentalTrees();
        return inferNetworkFromParentalTrees.inferNetwork(inferredTrees).toString();
    }
}
