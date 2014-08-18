package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPP;

import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.GTRModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.RateModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.matrixexponentiation.GenericKrylovMethod;
import jeigen.DenseMatrix;
import jeigen.Shortcuts;

public class MatrixQ{

    final int M;

    final RateModel rateModel;
    final DenseMatrix rateMatrix;


    final DenseMatrix transposedQ;
    final DenseMatrix equilibrium;

    final double theta;

    final int krylovIterations;

    public MatrixQ(RateModel rModel, int M, double theta){
        this.M = M;
        this.theta = theta;

        rateMatrix = rModel.getRateMatrix();
        rateModel = rModel;
        transposedQ = initializeMatrixQ(rateMatrix,M,theta);
        equilibrium = computeEquilibrium();

        krylovIterations = 20;
    }

    public MatrixQ(RateModel rModel, int M)
    {
        this(rModel, M, 0.02);
    }


    private static DenseMatrix initializeMatrixQ(DenseMatrix rateMatrix, int M, double theta){
        DenseMatrix transposedQ = Shortcuts.zeros(R.getMatrixSize(M), R.getMatrixSize(M));

        for(int n=1; n<=M; n++){
            for(R r : R.loopOver(n)){

                for (int fromType = 0; fromType < R.getNumberOfTypes(); fromType++) {
                    if (r.getNum(fromType) != 0)
                    {
                        for (int toType = 0; toType < R.getNumberOfTypes(); toType++) {
                            if (toType == fromType)
                                continue;

                            transposedQ.set(r.getIndex(),r.transition(fromType, toType).getIndex(), r.getNum(fromType) * rateMatrix.get(toType, fromType));
                        }
                    }
                }

                if (n != M)
                {
                    for (int splitType = 0; splitType < R.getNumberOfTypes(); splitType++) {
                        transposedQ.set(r.getIndex(),r.split(splitType).getIndex(), r.getNum(splitType)*(n+1)/theta);
                    }
                }

                double sum = -n*(n-1)/theta;

                for (int stayType = 0; stayType < R.getNumberOfTypes(); stayType++) {
                    sum += r.getNum(stayType) * rateMatrix.get(stayType,stayType);
                }

                transposedQ.set(r.getIndex(), r.getIndex(),sum);
            }
        }

        return transposedQ;


    }


    private GenericKrylovMethod.MatrixLike getScaled(final double timeScale)
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

        KrylovResult krylovResult = getProbabilityForColumnKrylov(time,actualColumn,krylovIterations);

        if (isGood(krylovResult))
            return krylovResult.result;
        else
        {
            System.out.println("Krylov is returning bad approximations. You might want to consider upping the iteration count in MatrixQ.getProbabilityForColumnKrylov.");
            DenseMatrix mexpResult = getProbabiltyForColumnMexp(time, actualColumn);
            //checkAccuracyReal(mexpResult,krylovResult);
            return mexpResult;
        }

    }

    private boolean isGood(KrylovResult result) {
        return result.errorGuess < .01;
    }

    private boolean isActuallyGood(DenseMatrix mexpResult, KrylovResult krylovResult)
    {
        DenseMatrix diff = krylovResult.result.sub(mexpResult);

        double diffMagnitude = diff.pow(2).sum().sqrt().s();

        return diffMagnitude < .01;

    }

    private void checkAccuracyReal(DenseMatrix mexpResult, KrylovResult krylovResult) {
        DenseMatrix diff = krylovResult.result.sub(mexpResult);

        double diffMagnitude = diff.pow(2).sum().sqrt().s();
        //double originalMagnitude = mexpResult.pow(2).sum().sqrt().s();
        //System.out.println(diff + " " + diffMagnitude);

        if (diffMagnitude > .01)
        {
            System.out.println("Really BAD :(");
        }

        if (diffMagnitude > (krylovResult.errorGuess + 1e-10))
            System.out.println("Error guess is  bad " + diffMagnitude + " vs " + krylovResult.errorGuess);
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

    private KrylovResult getProbabilityForColumnKrylov(double time, double[] actualColumn, int iterations) {

        GenericKrylovMethod method = new GenericKrylovMethod(getScaled(time),actualColumn,iterations);

        return new KrylovResult(method.getResult(),method.getEstimatedError());
    }

    private DenseMatrix getProbabiltyForColumnMexp(double time, double[] column) {
        int size = column.length;
        DenseMatrix columnCopy = new DenseMatrix(new double[][]{column}).t();
        DenseMatrix first = transposedQ.mul(time).cols(0,size).rows(0,size);
        return first.mexp().mmul(columnCopy);
    }


    private DenseMatrix computeEquilibrium() {

        DenseMatrix result = Shortcuts.ones(R.getMatrixSize(M),1);
        result = result.div(0);


        for (R r: R.loopOver(1))
        {
            result.set(r.getIndex(),0,rateModel.getEquilibriumVector().get(r.getType(),0));
        }


        for (int n = 2; n <= M; n++)
        {
            DenseMatrix coefficients = Shortcuts.zeros(R.getMatrixSizeWithoutN(n),R.getMatrixSizeWithoutN(n));
            DenseMatrix constants = Shortcuts.zeros(R.getMatrixSizeWithoutN(n),1);
            for (R r : R.loopOver(n))
            {
                double myCoeff = -n * (n-1)/theta;

                for (int stayType = 0; stayType < R.getNumberOfTypes(); stayType++) {
                    myCoeff += r.getNum(stayType) * rateMatrix.get(stayType,stayType);
                }
                coefficients.set(r.getIndexWithoutN(),r.getIndexWithoutN(),myCoeff);


                for (int fromType = 0; fromType < R.getNumberOfTypes(); fromType++) {
                    if (r.getNum(fromType) != 0)
                    {
                        for (int toType = 0; toType < R.getNumberOfTypes(); toType++) {
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
                    for (int coalesceType = 0; coalesceType < R.getNumberOfTypes(); coalesceType++) {
                        if (r.getNum(coalesceType) != 0)
                        {
                            constant += -(r.getNum(coalesceType)-1)/theta * n * result.get(r.coalesce(coalesceType).getIndex(),0);
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
//        //R.dims = 1;
//        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {5/6.0, 1/6.0}, new double[] {0.1});
//        System.out.println("The BiAllelicGTR Model is: " + BAGTRModel.getRateMatrix());
//        MatrixQ tester = new MatrixQ(BAGTRModel , 5);
//
//        System.out.println(tester.transposedQ);




        RateModel jc = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});
        for (int n = 1; n <= 100 ;n++)
        {
            FMatrix test = new FMatrix(n);
            test.set(R.loopOver(n).iterator().next(),1);
            MatrixQ q = new MatrixQ(jc,n);

            for (int a = 0; a< 100; a++)
            {
                q.transposedQ.mexp();//(1,test.arr);
            }

            long startTime = System.nanoTime();
            for (int a = 0; a< 100; a++)
            {
                q.transposedQ.mexp();//getProbabilityForColumn(1,test.arr);
            }
            long totalTime = System.nanoTime() - startTime;

            System.out.println(n + "," + totalTime/100.0);

        }
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
                        sum += vector[r.splitIndex(splitType)] * r.getNum(splitType) * (n + 1) / theta;
                    }
                }

                double selfSum = -n*(n-1)/theta;

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
