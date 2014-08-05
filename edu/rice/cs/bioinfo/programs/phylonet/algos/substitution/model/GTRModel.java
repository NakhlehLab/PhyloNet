package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model;

//Created 5-22-14, our own implementation different from GTRSubstitutionModel.

import jeigen.DenseMatrix;

import java.util.Arrays;
import java.util.Random;

/**
 * This class represents a GTR substition model.
 */
public class GTRModel extends SubstitutionModel
{

    private final DenseMatrix rateMatrix;
    private final DenseMatrix equilibriumVector;

    /**
     * Creates a GTR model with scaled transition frequencies.
     * Both equilibrium and pre-scaled transition frequencies need to sum to one.
     * @param scale The scale constant.
     * @param equilibriumFrequencies The equilibrium frequencies.
     * @param transitionFrequencies The pre-scaled transition frequencies.
     */
    public GTRModel(double scale, double[] equilibriumFrequencies, double[] transitionFrequencies)
    {
        this(equilibriumFrequencies, scale(scale, transitionFrequencies));

        if (Math.abs(sum(transitionFrequencies) - 1) > .005)
        {
            throw new IllegalArgumentException("The transition frequencies need to add to 1.");
        }
    }

    /**
     * Creates a GTR model with equilibrium and transition frequencies.
     * Equilibrium frequencies need to sum to one.
     * @param equilibriumFrequencies The equilibrium frequencies.
     * @param transitionFrequencies The transition frequencies.
     */
    public GTRModel(double[] equilibriumFrequencies, double[] transitionFrequencies)
    {
        if (equilibriumFrequencies.length != 4)
        {
            throw new IllegalArgumentException("Error: There should be 4 equilibrium frequencies.");
        }

        if (Math.abs(sum(equilibriumFrequencies) - 1) > .005)
        {
            throw new IllegalArgumentException("The equilibrium frequencies need to add to 1.");
        }

        if (transitionFrequencies.length != 6)
        {
            throw new IllegalArgumentException("Error: There should be 6 transition frequencies.");
        }


        rateMatrix = createRateMatrix(equilibriumFrequencies, transitionFrequencies);
        equilibriumVector = createEquilibriumVector(equilibriumFrequencies);
    }

    /**
     * Calculates the sum of an array.
     * @param arr An array.
     * @return The resulting sum.
     */
    private static double sum(double[] arr)
    {
        double result = 0;
        for (double v : arr)
        {
            result += v;
        }

        return result;
    }

    /**
     * Scales an array by a constant.
     * Does not mutate the passed in array.
     * @param scale The constant to scale by.
     * @param array The array to scale.
     * @return The scaled array.
     */
    private static double[] scale(double scale, double[] array)
    {
        double[] result = array.clone();
        for (int i = 0; i < result.length; i++)
        {
            result[i] *= scale;
        }

        return result;
    }


    /**
     * Creates an equilibrium vector given equilibrium frequencies.
     * @param e The equilibrium frequencies.
     * @return The equilibrium vector.
     */
    private static DenseMatrix createEquilibriumVector(double[] e)
    {
        double[][] temporary = {
                {e[0]},
                {e[1]},
                {e[2]},
                {e[3]}
        };
        return new DenseMatrix(temporary);
    }

    /**
     * Creates a rate matrix from the given equilibrium and transition frequencies.
     * @param e The equilibrium frequencies.
     * @param t The transition frequencies.
     * @return The rate matrix.
     */
    private static DenseMatrix createRateMatrix(double[] e, double[] t)
    {
        double[][] temporary = {{-(t[0] + t[1] + t[2]), e[0] * t[0] / e[1], e[0] * t[1] / e[2], e[0] * t[2] / e[3]},
                {t[0], -((e[0] * t[0] / e[1]) + t[3] + t[4]), e[1] * t[3] / e[2], e[1] * t[4] / e[3]},
                {t[1], t[3], -((e[0] * t[1] / e[2]) + (e[1] * t[3] / e[2]) + t[5]), e[2] * t[5] / e[3]},
                {t[2], t[4], t[5], -((e[0] * t[2] / e[3]) + (e[1] * t[4] / e[3]) + (e[2] * t[5] / e[3]))}};

        //Create the matrix Q, as described by WikiPedia's Substitution Model page.
        return new DenseMatrix(temporary);
    }

    @Override
    public DenseMatrix getEquilibriumVector()
    {
        return equilibriumVector;
    }

    @Override
    public DenseMatrix getRateMatrix()
    {
        return rateMatrix;
    }

    public static void main(String[] Args)
    {

        Random rand = new Random();


        double[] equalibriums = new double[4];

        for (int i = 0; i < 4; i++)
        {
            equalibriums[i] = rand.nextDouble();
        }

        double sum = sum(equalibriums);
        for (int i = 0; i < equalibriums.length; i++)
        {
            equalibriums[i] /= sum;
        }

        double[] transitions = new double[6];

        for (int i = 0; i < transitions.length; i++)
        {
            transitions[i] = rand.nextDouble() * 10;
        }

        System.out.println(Arrays.toString(equalibriums));
        GTRModel model = new GTRModel(equalibriums, transitions);

        System.out.println(model.rateMatrix);


        System.out.println(model.rateMatrix.mmul(model.getEquilibriumVector()));

        DenseMatrix testColumns = new DenseMatrix(new double[][]{{1, 1, 1, 1}});
        System.out.println(testColumns.mmul(model.rateMatrix));


        System.out.println(model.getProbabilityMatrix(1));

    }


}
