package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model;

import jeigen.DenseMatrix;

public class BiAllelicGTR implements RateModel {

    DenseMatrix rateMatrix;
    DenseMatrix equilibriumMatrix;

    private static double sum(double[] arr)
    {
        double result= 0;
        for (double v : arr) {
            result += v;
        }
        return result;
    }

    public BiAllelicGTR(double[] equilibriumFrequencies, double[] transitionFrequencies) {
        //Check if the parameters are the correct length.
        if (equilibriumFrequencies.length != 2) {
            throw new IllegalArgumentException("Error: There should be 2 equilibrium frequencies.");
        }

        if (Math.abs(sum(equilibriumFrequencies)-1) >.005) {
            throw new IllegalArgumentException("The equilibrium frequencies need to add to 1.");
        }

        if (transitionFrequencies.length != 1) {
            throw new IllegalArgumentException("Error: There should be 1 transition frequencies.");
        }

        createRateMatrix(equilibriumFrequencies,transitionFrequencies);
        createEquilibriumMatrix(equilibriumFrequencies);
    }

    private void createEquilibriumMatrix(double[] e) {
        double [][] temporary = {
                {e[0]},
                {e[1]}
        };
        equilibriumMatrix = new DenseMatrix(temporary);
    }

    private void createRateMatrix(double[] e, double[] t)
    {
        double[][] temporary = {{-t[0],  (e[0]/e[1])*t[0] },
                                { t[0], -(e[0]/e[1])*t[0] }};

        //Create the matrix Q, as described by WikiPedia's Substitution Model page.
        rateMatrix = new DenseMatrix(temporary);
    }

    @Override
    public DenseMatrix getEquilibriumVector(){
        return equilibriumMatrix;
    }

    @Override
    public DenseMatrix getRateMatrix() {
        return rateMatrix;    }


    public static void main(String[] args)
    {
        BiAllelicGTR f = new BiAllelicGTR(new double[]{.1, .9}, new double[]{.1});

        System.out.println(f.getEquilibriumVector());
        System.out.println(f.getRateMatrix());

        System.out.println(f.getRateMatrix().mmul(f.getEquilibriumVector()));
    }

}
