package edu.rice.cs.bioinfo.programs.phylonet.algos.clustering;

import java.util.Random;
import java.util.concurrent.RunnableFuture;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 7/4/16
 * Time: 3:24 PM
 * To change this template use File | Settings | File Templates.
 */
public class MultidimensionalScaling {
    private double _convergeThreshold;
    private int _iterations;

    public double prod(double[] x, double[] y) {
        double result = 0.0D;
        int length = Math.min(x.length, y.length);

        for(int i = 0; i < length; ++i) {
            result += x[i] * y[i];
        }

        return result;
    }

    public void multiply(double[][] matrix, double factor) {
        for(int i = 0; i < matrix.length; ++i) {
            for(int j = 0; j < matrix[0].length; ++j) {
                matrix[i][j] *= factor;
            }
        }

    }

    public double normalize(double[] x) {
        double norm = Math.sqrt(prod(x, x));

        for(int i = 0; i < x.length; ++i) {
            x[i] /= norm;
        }

        return norm;
    }


    public void randomize(double[][] matrix) {
        Random random = new Random(1L);

        for(int i = 0; i < matrix.length; ++i) {
            for(int j = 0; j < matrix[0].length; ++j) {
                matrix[i][j] = 0.5D - random.nextDouble();
            }
        }

    }

    public void squareEntries(double[][] matrix) {
        int n = matrix[0].length;
        int k = matrix.length;

        for(int i = 0; i < k; ++i) {
            for(int j = 0; j < n; ++j) {
                matrix[i][j] = Math.pow(matrix[i][j], 2.0);
            }
        }

    }

    public int[] landmarkIndices(double[][] matrix) {
        int k = matrix.length;
        int n = matrix.length;
        int[] result = new int[k];
        double eps = 1e-6;

        for(int i = 0; i < k; ++i) {
            for(int j = 0; j < n; ++j) {
                if(Math.abs(matrix[i][j]) < eps) {
                    result[i] = j;
                }
            }
        }

        return result;
    }

    public double[][] getLandmarkMatrix(double[][] matrix) {
        int n = matrix.length;
        double[][] result = new double[n][n];
        int[] index = landmarkIndices(matrix);

        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                result[i][j] = matrix[i][index[j]];
            }
        }

        return result;
    }

    public void doubleCenter(double[][] matrix) {
        int n = matrix[0].length;
        int k = matrix.length;

        for(int i = 0; i < k; ++i) {
            double avg = 0.0;

            for(int j = 0; j < n; ++j) {
                avg += matrix[i][j];
            }

            avg /= (double)n;

            for(int j = 0; j < n; ++j) {
                matrix[i][j] -= avg;
            }
        }

        for(int i = 0; i < n; ++i) {
            double avg = 0.0;

            for(int j = 0; j < k; ++j) {
                avg += matrix[j][i];
            }

            avg /= (double)matrix.length;

            for(int j = 0; j < k; ++j) {
                matrix[j][i] -= avg;
            }
        }

    }

    public void eigen(double[][] matrix, double[][] eigenvectors, double[] eigenvalues) {
        int d = eigenvalues.length;
        int n = matrix.length;
        double eps = 1e-6;

        for(int m = 0; m < d; ++m) {
            if(m > 0) {
                for(int i = 0; i < n; ++i) {
                    for(int j = 0; j < n; ++j) {
                        matrix[i][j] -= eigenvalues[m - 1] * eigenvectors[m - 1][i] * eigenvectors[m - 1][j];
                    }
                }
            }

            for(int i = 0; i < n; ++i) {
                eigenvectors[m][i] = Math.random();
            }

            normalize(eigenvectors[m]);
            double r = 0.0;

            for(int iter = 0; iter < 100; ++iter) {
                double[] q = new double[n];

                for(int i = 0; i < n; ++i) {
                    for(int j1 = 0; j1 < n; ++j1) {
                        q[i] += matrix[i][j1] * eigenvectors[m][j1];
                    }
                }

                eigenvalues[m] = prod(eigenvectors[m], q);
                normalize(q);
                r = Math.abs(prod(eigenvectors[m], q));
                eigenvectors[m] = q;

                if(Math.abs(1.0 - r) > eps)
                    break;
            }
        }

    }

    public void landmarkMDS(double[][] distMatrix, double[][] coordinates) {
        double[][] distances = new double[distMatrix.length][distMatrix[0].length];

        for(int i = 0; i < distances.length; ++i) {
            for(int j = 0; j < distances[0].length; ++j) {
                distances[i][j] = distMatrix[i][j];
            }
        }

        squareEntries(distances);
        int n = distances.length;
        int d = coordinates[0].length;
        double[] mean = new double[n];

        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                mean[i] += distances[j][i];
            }
        }

        for(int i = 0; i < n; ++i) {
            mean[i] /= (double)n;
        }

        double[] eigenvalues = new double[d];
        double[][] eigenvectors = new double[d][n];
        randomize(eigenvectors);
        double[][] landmarkMatrix = getLandmarkMatrix(distMatrix);
        squareEntries(landmarkMatrix);
        doubleCenter(landmarkMatrix);
        multiply(landmarkMatrix, -0.5);
        eigen(landmarkMatrix, eigenvectors, eigenvalues);

        for(int m = 0; m < eigenvectors.length; ++m) {
            for(int i = 0; i < eigenvectors[0].length; ++i) {
                eigenvectors[m][i] *= Math.sqrt(eigenvalues[m]);
            }
        }

        for(int m = 0; m < d; ++m) {
            for(int i = 0; i < n; ++i) {
                coordinates[i][m] = 0.0;

                for(int j = 0; j < n; ++j) {
                    coordinates[i][m] -= 0.5 * (distances[j][i] - mean[i]) * eigenvectors[m][j] / eigenvalues[m];
                }
            }
        }

    }

    public void majorize(double[][] distMatrix, double[][] coordinates) {
        double eps = 1e-6;

        if(distMatrix.length != distMatrix[0].length)
            throw new RuntimeException("Distance Matrix Error: wrong dimension " + distMatrix.length + ',' + distMatrix[0].length);

        if(distMatrix.length != coordinates.length) {
            throw new RuntimeException("Invaild coordinates array");
        }

        int n = coordinates.length;
        int dim = coordinates[0].length;
        int[] index = landmarkIndices(distMatrix);

        double[][] weightMatrix = new double[n][n];
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                if(distMatrix[i][j] > eps) {
                    weightMatrix[i][j] = Math.pow(distMatrix[i][j], 0.0);
                }
            }
        }

        double[] weightSum = new double[n];

        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                weightSum[i] += weightMatrix[j][i];
            }
        }

        for(int round = 0; round < _iterations; ++round) {
            double delta = 0.0;
            double magnitude = 0.0;

            for(int i = 0; i < n; ++i) {
                double[] xnew = new double[dim];

                for(int m = 0; m < n; ++m) {
                    double inv = 0.0;

                    for(int m1 = 0; m1 < dim; ++m1) {
                        inv += Math.pow(coordinates[i][m1] - coordinates[index[m]][m1], 2.0);
                    }

                    if(Math.abs(inv) > eps) {
                        inv = Math.pow(inv, -0.5);
                    }

                    for(int m1 = 0; m1 < dim; ++m1) {
                        xnew[m1] += weightMatrix[m][i] * (coordinates[index[m]][m1] + distMatrix[m][i] * (coordinates[i][m1] - coordinates[index[m]][m1]) * inv);
                    }
                }

                if(Math.abs(weightSum[i]) > eps) {
                    for(int m = 0; m < dim; ++m) {
                        delta += Math.pow(xnew[m] / weightSum[i] - coordinates[i][m], 2.0);
                        magnitude += Math.pow(coordinates[i][m], 2.0);
                        coordinates[i][m] = xnew[m] / weightSum[i];
                    }
                }
            }

            delta = Math.sqrt(delta / magnitude);
            if(delta < _convergeThreshold) {
                return ;
            }

        }

    }

    public void getMDScoordinates(int n, int d, double[][] distMatrix, double [][] coordinates) {
        _iterations = 1000;
        _convergeThreshold = Math.pow(10.0, -3);

        randomize(coordinates);
        landmarkMDS(distMatrix, coordinates);
        majorize(distMatrix, coordinates);


    }
}
