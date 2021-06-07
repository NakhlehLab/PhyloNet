package edu.rice.cs.bioinfo.programs.phylonet.algos.clustering;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 7/28/16
 * Time: 11:56 AM
 * To change this template use File | Settings | File Templates.
 */
public class HungarianMatching {
    private double[][] _costMatrix;
    private int _rows, _cols, _dim;


    public void setCostMatrix(double[][] costMatrix) {
        _dim = Math.max(costMatrix.length, costMatrix[0].length);
        _rows = costMatrix.length;
        _cols = costMatrix[0].length;
        _costMatrix = new double[_dim][_dim];
        for (int w = 0; w < _dim; w++) {
            if (w < costMatrix.length) {
                if (costMatrix[w].length != _cols) {
                    throw new IllegalArgumentException("Illegal dimension of cost matrix");
                }
                _costMatrix[w] = Arrays.copyOf(costMatrix[w], _dim);
            } else {
                _costMatrix[w] = new double[_dim];
            }
        }

    }

    public HungarianMatching() {

    }

    public HungarianMatching(double[][] costMatrix) {
        setCostMatrix(costMatrix);
    }

    //execute algorithm
    public int[] compute() {
        //initialize arrays

        int[] workerToJob = new int[_dim];
        int[] jobToWorker = new int[_dim];
        int[] jobToParentWorker = new int[_dim];
        int[] jobToNextWorker = new int[_dim];
        double[] jobToNextWorkerValue = new double[_dim];
        boolean[] workerUsed = new boolean[_dim];
        double[] workerLabel = new double[_dim];
        double[] jobLabel = new double[_dim];

        //initialize first solution

        Arrays.fill(workerToJob, -1);
        Arrays.fill(jobToWorker, -1);

        for (int w = 0; w < _dim; w++) {
            double min = Double.POSITIVE_INFINITY;
            for (int j = 0; j < _dim; j++) {
                if (_costMatrix[w][j] < min) {
                    min = _costMatrix[w][j];
                }
            }
            for (int j = 0; j < _dim; j++) {
                _costMatrix[w][j] -= min;
            }
        }
        double[] min = new double[_dim];
        for (int j = 0; j < _dim; j++) {
            min[j] = Double.POSITIVE_INFINITY;
        }
        for (int w = 0; w < _dim; w++) {
            for (int j = 0; j < _dim; j++) {
                if (_costMatrix[w][j] < min[j]) {
                    min[j] = _costMatrix[w][j];
                }
            }
        }
        for (int w = 0; w < _dim; w++) {
            for (int j = 0; j < _dim; j++) {
                _costMatrix[w][j] -= min[j];
            }
        }

        for (int j = 0; j < _dim; j++) {
            jobLabel[j] = Double.POSITIVE_INFINITY;
        }
        for (int w = 0; w < _dim; w++) {
            for (int j = 0; j < _dim; j++) {
                if (_costMatrix[w][j] < jobLabel[j]) {
                    jobLabel[j] = _costMatrix[w][j];
                }
            }
        }

        for (int w = 0; w < _dim; w++) {
            for (int j = 0; j < _dim; j++) {
                if (workerToJob[w] == -1 && jobToWorker[j] == -1 && _costMatrix[w][j] - workerLabel[w] - jobLabel[j] == 0) {
                    workerToJob[w] = j;
                    jobToWorker[j] = w;
                }
            }
        }

        //find a unmatched worker
        int w;
        for (w = 0; w < _dim; w++) {
            if (workerToJob[w] == -1) {
                break;
            }
        }

        //begin iteration
        while (w < _dim) {
            Arrays.fill(workerUsed, false);
            Arrays.fill(jobToParentWorker, -1);
            workerUsed[w] = true;
            for (int j = 0; j < _dim; j++) {
                jobToNextWorkerValue[j] = _costMatrix[w][j] - workerLabel[w] - jobLabel[j];
                jobToNextWorker[j] = w;
            }

            while (true) {
                int bestWorker = -1, bestJob = -1;
                double minv = Double.POSITIVE_INFINITY;
                for (int j = 0; j < _dim; j++) {
                    if (jobToParentWorker[j] == -1) {
                        if (jobToNextWorkerValue[j] < minv) {
                            minv = jobToNextWorkerValue[j];
                            bestWorker = jobToNextWorker[j];
                            bestJob = j;
                        }
                    }
                }
                if (minv > 0) {
                    for (int i = 0; i < _dim; i++) {
                        if (workerUsed[i]) {
                            workerLabel[i] += minv;
                        }
                    }
                    for (int j = 0; j < _dim; j++) {
                        if (jobToParentWorker[j] != -1) {
                            jobLabel[j] -= minv;
                        } else {
                            jobToNextWorkerValue[j] -= minv;
                        }
                    }
                }
                jobToParentWorker[bestJob] = bestWorker;
                if (jobToWorker[bestJob] == -1) {
                    int committedJob = bestJob;
                    int parentWorker = jobToParentWorker[committedJob];
                    while (true) {
                        int temp = workerToJob[parentWorker];
                        workerToJob[parentWorker] = committedJob;
                        jobToWorker[committedJob] = parentWorker;
                        committedJob = temp;
                        if (committedJob == -1) {
                            break;
                        }
                        parentWorker = jobToParentWorker[committedJob];
                    }
                    break;
                } else {

                    int worker = jobToWorker[bestJob];
                    workerUsed[worker] = true;
                    for (int j = 0; j < _dim; j++) {
                        if (jobToParentWorker[j] == -1) {
                            double slack = _costMatrix[worker][j] - workerLabel[worker] - jobLabel[j];
                            if (jobToNextWorkerValue[j] > slack) {
                                jobToNextWorkerValue[j] = slack;
                                jobToNextWorker[j] = worker;
                            }
                        }
                    }
                }
            }

            //find next unassigned worker
            for (w = 0; w < _dim; w++) {
                if (workerToJob[w] == -1) {
                    break;
                }
            }
        }
        int[] result = Arrays.copyOf(workerToJob, _rows);
        for (w = 0; w < result.length; w++) {
            if (result[w] >= _cols) {
                result[w] = -1;
            }
        }
        return result;
    }

}