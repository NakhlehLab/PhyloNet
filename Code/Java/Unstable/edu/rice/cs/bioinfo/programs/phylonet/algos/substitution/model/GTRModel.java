package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model;

//Created 5-22-14, our own implementation different from GTRSubstitutionModel.

import jeigen.DenseMatrix;
//import jeigen.DenseMatrix.PseudoEigenResult;
import Jama.Matrix;


import java.util.Arrays;


/**
 * This class represents a GTR substition model.
 */
public class GTRModel extends SubstitutionModel
{

    private final DenseMatrix rateMatrix;
    private final DenseMatrix equilibriumVector;
    DenseMatrix rateMatrixIntegrated;


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
        rateMatrixIntegrated = setProbabilityMatrixIntegrated();
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

    //@Override
    public DenseMatrix getEquilibriumVector()
    {
        return equilibriumVector;
    }

    //@Override
    public DenseMatrix getRateMatrix()
    {
        return rateMatrix;
    }

    public DenseMatrix getProbabilityMatrixIntegrated(){
        return rateMatrixIntegrated;
    }


    public DenseMatrix setProbabilityMatrixIntegrated(){

        int rows = rateMatrix.rows;
        int cols = rateMatrix.cols;
        //todo
        /*
        PseudoEigenResult eig = rateMatrix.peig();
        DenseMatrix D = eig.values;
        DenseMatrix V = eig.vectors;
        DenseMatrix Vinv = inverse(eig.vectors);
        */
    	
    	
    	/*
    	System.out.println("values");
    	System.out.println(D);
    	System.out.println("vectors");
    	System.out.println(V);
    	System.out.println("vectors inverse");
    	System.out.println(Vinv);    	
    	
    	
    	for (int i = 0; i < rows; i++) {
    		D.set(i, i, Math.exp(2 * D.get(i, i)));
    	}
    	System.out.println("New values");
    	System.out.println(D);  	
    	
    	
    	System.out.println(V.mmul(D).mmul(Vinv));
    	System.out.println("Inside");
    	*/


        double [][] temp = new double[rows][cols];
        double pivot = rows;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                for (int k = 0; k < pivot; k++) {
                    //todo
                    //temp[i][j] += V.get(i, k)* Vinv.get(k, j)* (1/(-D.get(k,k)+1)) ;
                }
            }
        }

        return new DenseMatrix(temp);
    }



//  //todo
    /*
    public DenseMatrix inverse (DenseMatrix M) {
        // Use JAMA instead of Jeigen, as it doesn't have support for inverse.
        int rows = M.rows;
        int cols = M.cols;
        double [] values = M.getValues();


        //Convert "Jeigen" to "JAMA".
        double [][] array = new double[rows][cols];
        int count = 0;
        for(int i = 0; i < cols; i++) {
            for (int j = 0; j < rows; j++) {
                array[j][i] = values[count++];
            }
        }
        Matrix A = new Matrix(array);
        Matrix Ainv = A.inverse();

        //Convert "JAMA" to "Jeigen".
        double [][] temp = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                temp[i][j] = Ainv.get(i, j);
            }
        }

        return new DenseMatrix(temp);
    }
    */




    public static void main(String[] Args)
    {
        double[] equilibriums = {.1,.2,.3,.4};
        double[] transitions = {.4,.9,.3,.2,.5,.2};

        //double[] equilibriums = {.25,.25,.25,.25};
        //double[] transitions = {.025,.025,.025,.025,.025,.025};    	

        System.out.println(Arrays.toString(equilibriums));
        GTRModel model = new GTRModel(equilibriums, transitions);

        //System.out.println(model.rateMatrix);


        //System.out.println(model.rateMatrix.mmul(model.getEquilibriumVector()));

        //DenseMatrix testColumns = new DenseMatrix(new double[][]{{1, 1, 1, 1}});
        //System.out.println(testColumns.mmul(model.rateMatrix));

        System.out.println(model.getRateMatrix());
        System.out.println(model.getProbabilityMatrixIntegrated());
    }


}
