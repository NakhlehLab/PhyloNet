package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model;

import jeigen.DenseMatrix;

public class JCModel extends SubstitutionModel {
    DenseMatrix rateMatrix;
    double m;	//rate or mu.
    DenseMatrix rateMatrixIntegrated;


    public JCModel(double[] equilibriumFrequencies, double[] transitionFrequencies) {
        //Pass it to next constructor.
        this(transitionFrequencies[0]);
    }

    public JCModel(double rate)
    {
        //System.out.println("Rate: " + rate);
        m = rate;
        rateMatrix = createMatrix(rate);
        rateMatrixIntegrated = setProbabilityMatrixIntegrated();

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
        //System.out.println("Temp Matrix:\n" + temp);
        //System.out.println("Rate Matrix:\n" + new DenseMatrix(temp));

        return new DenseMatrix(temp);

    }


    //@Override
    // It is a column vector.
    public DenseMatrix getEquilibriumVector() {
        return new DenseMatrix(new double[][] {{.25},{.25},{.25},{.25}});
    }

    //@Override
    public DenseMatrix getRateMatrix() {
        return rateMatrix;
    }



    public DenseMatrix setProbabilityMatrixIntegrated(){

        double[][] temp = new double[4][4];

        for (int y = 0; y < temp.length; y++) {
            double[] row = temp[y];
            for (int x = 0; x < row.length; x++) {
                if (x == y)
                    temp[y][x] = (m+4)/ (4 * (m+1));
                else
                    temp[y][x] = (m)/ (4 * (m+1));
            }
        }

        return new DenseMatrix(temp);
    }
    public DenseMatrix getProbabilityMatrixIntegrated(){
        return rateMatrixIntegrated;
    }

    public static void main(String[] args)
    {
        JCModel foo = new JCModel(new double[] {.25,.25,.25,.25},new double[] {.1,0,0,0,0,0});

        System.out.println("Rate Matrix:");
        System.out.println(foo.getRateMatrix());
        System.out.println("Probability Matrix:");
        System.out.println(foo.getProbabilityMatrix(2));
        System.out.println("Integrated Matrix:");
        System.out.println(foo.getProbabilityMatrixIntegrated());

    }
}
