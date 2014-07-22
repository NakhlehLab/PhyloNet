package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import edu.rice.cs.bioinfo.programs.phylonet.algos.matrixexponentiation.GenericKrylovMethod;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.RateModel;
import jeigen.DenseMatrix;
import jeigen.Shortcuts;

import java.util.Arrays;

public class MatrixQ{

    final int M;

    final RateModel rateModel;
    final DenseMatrix rateMatrix;


    final DenseMatrix transposedQ;
    final DenseMatrix equilibrium;


    public MatrixQ(RateModel rModel, int M){
        this.M = M;
        rateMatrix = rModel.getRateMatrix();
        rateModel = rModel;
        transposedQ = initializeMatrixQ(rateMatrix,M);
        equilibrium = computeEquilibrium();

    }



    private static DenseMatrix initializeMatrixQ(DenseMatrix rateMatrix, int M){
        DenseMatrix transposedQ = Shortcuts.zeros(R.getMatrixSize(M), R.getMatrixSize(M));

        for(int n=1; n<=M; n++){
            for(R r : R.loopOver(n)){

                for (int fromType = 0; fromType < r.getNumberOfTypes(); fromType++) {
                    if (r.getNum(fromType) != 0)
                    {
                        for (int toType = 0; toType < r.getNumberOfTypes(); toType++) {
                            if (toType == fromType)
                                continue;

                            transposedQ.set(r.getIndex(),r.transition(fromType, toType).getIndex(), r.getNum(fromType) * rateMatrix.get(toType, fromType));
                        }
                    }
                }

                if (n != M)
                {
                    for (int splitType = 0; splitType < r.getNumberOfTypes(); splitType++) {
                        transposedQ.set(r.getIndex(),r.split(splitType).getIndex(), r.getNum(splitType)*(n+1)/2.0);
                    }
                }

                double sum = -n*(n-1)/2.0;

                for (int stayType = 0; stayType < r.getNumberOfTypes(); stayType++) {
                    sum += r.getNum(stayType) * rateMatrix.get(stayType,stayType);
                }

                transposedQ.set(r.getIndex(), r.getIndex(),sum);
            }
        }

        return transposedQ;


    }


    public double[][] getProbabilityForMatrix(double time, double[][] matrix)
    {

        double[][] resultMatrix = new double[matrix.length][matrix[0].length];

        for(int i=0; i<matrix[0].length; i++){
            double[] column = new double[matrix.length];
            for(int j=0; j<matrix.length; j++){
                column[j] = matrix[j][i];
            }
            DenseMatrix resultColumn = getProbabilityForColumn(time, column);
            for (int j = 0; j < resultColumn.rows; j++)
                resultMatrix[j][i] = resultColumn.get(j,0);
        }

        return resultMatrix;
    }


    public GenericKrylovMethod.MatrixLike getScaled(final double timeScale)
    {
        return new GenericKrylovMethod.MatrixLike()
        {
            @Override
            public double[] multiplyByVector(double[] vector)
            {
                return multiplyByVectorScaled(vector,timeScale);
            }
        };
    }


    public DenseMatrix getProbabilityForColumn(double time, double[] actualColumn)
    {
        //KrylovResult krylovResult = getProbabilityForColumnKrylov(time,actualColumn);
        //boolean isGood = isGood(krylovResult);
        KrylovResult krylovResult = null;
        boolean isGood = false;
        if (isGood)
            return krylovResult.result;
        else
        {
            //System.out.println(Arrays.toString(actualColumn) + " Krylov is returning bad approximations. You might want to consider upping the iteration count in getProbabilityForColumnKrylov.");
            DenseMatrix mexpResult = getProbabiltyForColumnMexp(time, actualColumn);
            //checkAccuracyReal(mexpResult,krylovResult);
            return mexpResult;
        }

    }

    private boolean isGood(KrylovResult result) {
        return result.errorGuess < .01;
    }



    private void checkAccuracyReal(DenseMatrix mexpResult, KrylovResult krylovResult) {
        DenseMatrix diff = krylovResult.result.sub(mexpResult);

        double diffMagnitude = diff.pow(2).sum().sqrt().s();
        //double originalMagnitude = mexpResult.pow(2).sum().sqrt().s();
        //System.out.println(diff + " " + diffMagnitude);

        if (diffMagnitude > .001)
        {
            System.out.print(' ');
        }

        if (diffMagnitude > (krylovResult.errorGuess + 1e-10))
            System.out.println("Actually bad");
    }



    private static class KrylovResult
    {
        final double errorGuess;
        final DenseMatrix result;

        public KrylovResult(DenseMatrix result, double errorGuess) {
            this.result = result;
            this.errorGuess = errorGuess;
        }
    }

    private KrylovResult getProbabilityForColumnKrylov(double time, double[] actualColumn) {

        GenericKrylovMethod method = new GenericKrylovMethod(getScaled(time),actualColumn,5);

        return new KrylovResult(method.getResult(),method.getEstimatedError());
    }

    private DenseMatrix getProbabiltyForColumnMexp(double time, double[] column) {
        int size = column.length;
        DenseMatrix columnCopy = new DenseMatrix(new double[][]{column}).t();
        DenseMatrix first = transposedQ.mul(time).cols(0,size).rows(0,size);
        return first.mexp().mmul(columnCopy);
    }


    public DenseMatrix computeEquilibrium() {

        DenseMatrix result = Shortcuts.ones(R.getMatrixSize(M),1);
        result = result.div(0);


        for (R r: R.loopOver(1))
        {
            result.set(r.getIndex(),0,rateModel.getEquilibriumVector().get(r.getType(),0));
        }


        for (int n = 2; n <= M; n++)
        {
            DenseMatrix coefficients = Shortcuts.zeros(R.getMatrixSizeWithoutN(n), R.getMatrixSizeWithoutN(n));
            DenseMatrix constants = Shortcuts.zeros(R.getMatrixSizeWithoutN(n),1);
            for (R r : R.loopOver(n))
            {
                double myCoeff = -n * (n-1)/2.0;

                for (int stayType = 0; stayType < r.getNumberOfTypes(); stayType++) {
                    myCoeff += r.getNum(stayType) * rateMatrix.get(stayType,stayType);
                }
                coefficients.set(r.getIndexWithoutN(),r.getIndexWithoutN(),myCoeff);


                for (int fromType = 0; fromType < r.getNumberOfTypes(); fromType++) {
                    if (r.getNum(fromType) != 0)
                    {
                        for (int toType = 0; toType < r.getNumberOfTypes(); toType++) {
                            if (toType == fromType)
                                continue;

                            double prev = (r.getNum(toType)+1)*rateMatrix.get(fromType,toType);
                            coefficients.set(r.getIndexWithoutN(),r.transition(fromType,toType).getIndexWithoutN(),prev);
                        }
                    }
                }

                if (n != 1)
                {
                    double constant = 0;
                    for (int coalesceType = 0; coalesceType < r.getNumberOfTypes(); coalesceType++) {
                        if (r.getNum(coalesceType) != 0)
                        {
                            constant += -(r.getNum(coalesceType)-1)/2.0 * n * result.get(r.coalesce(coalesceType).getIndex(),0);
                        }
                    }
                    constants.set(r.getIndexWithoutN(), 0, constant);
                }
            }

            DenseMatrix myResult = coefficients.fullPivHouseholderQRSolve(constants);

            for (R r: R.loopOver(n))
            {
                result.set(r.getIndex(),0,myResult.get(r.getIndexWithoutN(),0));
            }

        }
        return result;
    }

    public DenseMatrix getEquilibrium()
    {
        return equilibrium;
    }

    public static void main(String[] args){
        //R.dims = 1;
        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {5/6.0, 1/6.0}, new double[] {0.1});
        System.out.println("The BiAllelicGTR Model is: " + BAGTRModel.getRateMatrix());
        MatrixQ tester = new MatrixQ(BAGTRModel , 5);

        System.out.println(tester.transposedQ);
    }


    private double[] multiplyByVectorScaled(double[] vector, double timeScale) {

        final int numberOfTypes = R.getNumberOfTypes();

        int actualMaxSize = 1;

        while (R.getMatrixSize(actualMaxSize) != vector.length)
            actualMaxSize++;


        double[] result = new double[R.getMatrixSize(actualMaxSize)];

        for(int n=1; n<=actualMaxSize; n++){
            for(R r : R.loopOver(n)){

                double sum = 0;
                for (int fromType = 0; fromType < numberOfTypes; fromType++) {
                    if (r.getNum(fromType) != 0)
                    {
                        for (int toType = 0; toType < numberOfTypes; toType++) {
                            if (toType == fromType)
                                continue;

                            sum += vector[r.transitionIndex(fromType,toType)] * r.getNum(fromType) * rateMatrix.get(toType, fromType);
                        }
                    }
                }

                if (n != actualMaxSize)
                {
                    for (int splitType = 0; splitType < numberOfTypes; splitType++) {
                        sum += vector[r.splitIndex(splitType)] * r.getNum(splitType) * (n + 1) / 2.0;
                    }
                }

                double selfSum = -n*(n-1)/2.0;

                for (int stayType = 0; stayType < numberOfTypes; stayType++) {
                    selfSum += r.getNum(stayType) * rateMatrix.get(stayType,stayType);
                }

                sum += vector[r.getIndex()] * selfSum;

                result[r.getIndex()] = sum*timeScale;
            }
        }
        return result;

    }
}
