package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model;

import jeigen.DenseMatrix;

public class JCModel extends SubstitutionModel {
    DenseMatrix rateMatrix;

    public JCModel(double[] equilibriumFrequencies, double[] transitionFrequencies) {
        this(transitionFrequencies[0]);
    }

    public JCModel(double rate)
    {
        rateMatrix = createMatrix(rate);
    }

    DenseMatrix createMatrix(double transition)
    {
        double[][] temp = new double[4][4];

        for (int y = 0; y < temp.length; y++) {
            double[] row = temp[y];
            for (int x = 0; x < row.length; x++) {
                if (x == y)
                    temp[y][x] = -3.0/4 * transition;
                else
                    temp[y][x] = 1.0/4 * transition;
            }
        }
        return new DenseMatrix(temp);
    }


    @Override
    public DenseMatrix getEquilibriumVector() {
        return new DenseMatrix(new double[][] {{.25},{.25},{.25},{.25}});
    }

    @Override
    public DenseMatrix getRateMatrix() {
        return rateMatrix;
    }


    public static void main(String[] args)
    {
        JCModel foo = new JCModel(new double[] {.25,.25,.25,.25},new double[] {.1,0,0,0,0,0});

        System.out.println(foo.getProbabilityMatrix(2));
    }
}
