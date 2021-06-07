package edu.rice.cs.bioinfo.programs.phylonet.algos.matrixexponentiation;

import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.GTRModel;
import jeigen.DenseMatrix;

import java.util.Arrays;

/**
 * This class implements the Krylov method of matrix exponentiation.
 * Calculates: y = e^A * V
 * Please check the estimatedError after using to make sure the results are reasonable.
 *
 * Source:
 * http://www-users.cs.umn.edu/~saad/PDF/RIACS-90-ExpTh.pdf
 */
public class GenericKrylovMethod
{
    /**
     * This is the interface for a matrix type data structure which this method can work on.
     */
    public interface MatrixLike
    {
        /**
         * Multiplies this matrix by a column vector, like A * vector = result
         *
         * @param vector A column vector.
         * @return A row vector for the result.
         */
        double[] multiplyByVector(double[] vector);
    }

    /**
     * Creates a MatrixLike wrapper for a DenseMatrix.
     *
     * @param matrix The matrix to wrap.
     * @return A MatrixLike object representing this matrix.
     */
    public static MatrixLike onADense(final DenseMatrix matrix)
    {
        return new MatrixLike()
        {
            @Override
            public double[] multiplyByVector(double[] vector)
            {
                DenseMatrix vectorCopy = new DenseMatrix(new double[][]{vector}).t();
                DenseMatrix result1 = matrix.mmul(vectorCopy);

                double[] actualResutl = new double[result1.rows];
                for (int i = 0; i < result1.rows; i++)
                {
                    actualResutl[i] = result1.get(i, 0);
                }
                return actualResutl;
            }
        };
    }


    private MatrixLike matrixA;
    private double[] vectorV;
    private int M;

    private DenseMatrix Hm;
    private double[][] Vm;


    private double errorFactor;
    private double vectorVnorm;

    private DenseMatrix result;
    private double estimatedError;

    /**
     * Creates and runs the Krylov method on a matrix and vector with the given number of iterations.
     * Calculates: y = e^A * V
     * @param matrixA The matrix to do the exponentiation on.
     * @param vectorV The vector to multiply by.
     * @param iterations The number of iterations.
     */
    public GenericKrylovMethod(MatrixLike matrixA, double[] vectorV, int iterations)
    {
        this.matrixA = matrixA;
        this.vectorV = vectorV.clone();
        M = iterations;
        ArnoldiAlgorithm();

        DenseMatrix y = Hm.mexp().col(0).mul(vectorVnorm);

        DenseMatrix VmDense = new DenseMatrix(Vm).t();
        result = VmDense.mmul(y);

        estimatedError = errorFactor * y.get(y.rows - 1, 0);
    }


    private static double getNorm(double[] vec)
    {
        double sum = 0;
        for (double elem : vec)
        {
            sum += elem * elem;
        }

        return Math.sqrt(sum);
    }

    private static void scale(double[] vec, double scale)
    {
        for (int i = 0; i < vec.length; i++)
        {
            vec[i] *= scale;
        }
    }

    private static double dotProduct(double[] a, double[] b)
    {
        double sum = 0;
        for (int i = 0; i < a.length; i++)
        {
            sum += a[i] * b[i];
        }
        return sum;
    }

    private static void subScale(double[] a, double[] b, double scale)
    {
        for (int i = 0; i < a.length; i++)
        {
            a[i] -= b[i] * scale;
        }
    }

    private void ArnoldiAlgorithm()
    {
        //Initialization.
        Hm = new DenseMatrix(M, M);
        Vm = new double[M][];

        vectorVnorm = getNorm(vectorV);
        scale(vectorV, 1 / vectorVnorm);

        Vm[0] = vectorV;

        //Iterate.
        for (int j = 0; j < M; j++)
        {
            double[] currentVector = matrixA.multiplyByVector(Vm[j]);
            for (int i = 0; i <= j; i++)
            {
                Hm.set(i, j, dotProduct(Vm[i], currentVector));
                subScale(currentVector, Vm[i], Hm.get(i, j));
            }

            //Check if done.
            double currentVectorNorm = getNorm(currentVector);
            if (currentVectorNorm == 0)
            {
                Hm = Hm.cols(0, j + 1).rows(0, j + 1);
                Vm = Arrays.copyOf(Vm, j + 1);
                return;
            }

            if (j != M - 1)
            {
                Hm.set(j + 1, j, currentVectorNorm);

                scale(currentVector, 1 / currentVectorNorm);
                Vm[j + 1] = currentVector;
            } else
            {
                errorFactor = currentVectorNorm;
            }
        }
    }

    /**
     * Gets the result of the matrix exponentiation and multiplication.
     * @return The resulting column vector.
     */
    public DenseMatrix getResult()
    {
        return result;
    }

    /**
     * Gets the a number which approximates the error (to within roughly an order of magnitude).
     * @return The estimated error.
     */
    public double getEstimatedError()
    {
        return estimatedError;
    }


    public static void main(String[] args)
    {
        testKrylov();
    }

    private static void testKrylov()
    {

        GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4}, new double[]{1, 1.5, 2, 2.5, 3, 3.5});
        DenseMatrix rateMatrix = gtrModel.getRateMatrix();

        double[] specialVector2 = new double[]{5, 8, 2, 7};
        DenseMatrix specialVector = new DenseMatrix(4, 1);
        specialVector.set(0, 0, 5);
        specialVector.set(1, 0, 8);
        specialVector.set(2, 0, 2);
        specialVector.set(3, 0, 7);
        System.out.println("Rate Matrix: " + rateMatrix);
        System.out.println("Old MX = e^M * vectorV = " + rateMatrix.mexp().mmul(specialVector));

        long startTime = System.currentTimeMillis();
        GenericKrylovMethod specialK = new GenericKrylovMethod(onADense(rateMatrix), specialVector2, 1000);
        System.out.println("New MX = Krylov Method = " + specialK.getResult());
        System.out.println(System.currentTimeMillis() - startTime);

    }


}
