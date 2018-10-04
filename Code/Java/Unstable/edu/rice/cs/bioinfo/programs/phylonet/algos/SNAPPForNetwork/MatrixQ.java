package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import Jama.Matrix;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.RateModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.matrixexponentiation.GenericKrylovMethod;
import jeigen.DenseMatrix;
import jeigen.Shortcuts;
import org.apache.commons.math3.complex.Complex;

public class MatrixQ{

    final int M;

    final RateModel rateModel;
    final DenseMatrix rateMatrix;


    final DenseMatrix transposedQ;
    final DenseMatrix equilibrium;

    final double theta;

    double _time;
    DenseMatrix _expQ;

    double _u;
    double _v;

    public MatrixQ(RateModel rModel, int M, double theta){
        this.M = M;
        this.theta = theta;

        rateMatrix = rModel.getRateMatrix();
        rateModel = rModel;
        transposedQ = initializeMatrixQ(rateMatrix,M,theta);
        equilibrium = computeEquilibrium();

        if(rateModel instanceof BiAllelicGTR) {
            BiAllelicGTR ba = (BiAllelicGTR) rateModel;
            _u = ba.getU();
            _v = ba.getV();
        }
    }

    public MatrixQ(RateModel rModel, int M, double theta, boolean computeMatrix){
        this.M = M;
        this.theta = theta;

        rateMatrix = rModel.getRateMatrix();
        rateModel = rModel;

        if(computeMatrix) {
            transposedQ = initializeMatrixQ(rateMatrix, M, theta);
            equilibrium = computeEquilibrium();
        } else {
            transposedQ = null;
            equilibrium = null;
        }

        if(rateModel instanceof BiAllelicGTR) {
            BiAllelicGTR ba = (BiAllelicGTR) rateModel;
            _u = ba.getU();
            _v = ba.getV();
        }
    }

    public MatrixQ(RateModel rModel, int M)
    {
        this(rModel,M,1);
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

    public void setTime(double time, int maxColumnSize) {
        _time = time;
        //maxColumnSize = transposedQ.cols;
        _expQ = transposedQ.mul(time).cols(0,maxColumnSize).rows(0,maxColumnSize).mexp();
        //_expQ = transposedQ.mul(time).mexp();
    }

    public DenseMatrix getProbabilityForColumn(double[] actualColumn)
    {
        int size = actualColumn.length;
        DenseMatrix columnCopy = new DenseMatrix(new double[][]{actualColumn}).t();
        DenseMatrix ret = _expQ.cols(0,size).rows(0,size).mmul(columnCopy);
        return ret;

    }

    public DenseMatrix getProbabilityForColumn(double time, double[] actualColumn)
    {
        DenseMatrix mexpResult = getProbabiltyForColumnMexp(time, actualColumn);
        //checkAccuracyReal(mexpResult,krylovResult);
        return mexpResult;
        /*
        KrylovResult krylovResult = getProbabilityForColumnKrylov(time,actualColumn);
        //KrylovResult krylovResult = null;
        if (isGood(krylovResult))
            return krylovResult.result;
        else
        {
            //System.out.println("Krylov is returning bad approximations. You might want to consider upping the iteration count in MatrixQ.getProbabilityForColumnKrylov.");
            DenseMatrix mexpResult = getProbabiltyForColumnMexp(time, actualColumn);
            //checkAccuracyReal(mexpResult,krylovResult);
            return mexpResult;
        }*/

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

        GenericKrylovMethod method = new GenericKrylovMethod(getScaled(time),actualColumn,10);

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
            DenseMatrix coefficients = Shortcuts.zeros(R.getMatrixSizeWithoutN(n),R.getMatrixSizeWithoutN(n));
            DenseMatrix constants = Shortcuts.zeros(R.getMatrixSizeWithoutN(n),1);
            for (R r : R.loopOver(n))
            {
                double myCoeff = -n * (n-1)/ theta;

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

    private void solveOneBlock(ComplexCP[] y, ComplexCP offset, int n, ComplexCP[] x ) {
        ComplexCP t = new ComplexCP();
        ComplexCP d[] = new ComplexCP[n + 1];
        ComplexCP e[] = new ComplexCP[n + 1];

        for(int i = 0 ; i < n + 1 ; i++) {
            d[i] = new ComplexCP();
            e[i] = new ComplexCP();
        }
        d[0]._re = offset._re -1.0 * n * (n - 1) / theta - n * _v;
        d[0]._im = offset._im;
        e[0]._re = y[0]._re;
        e[0]._im = y[0]._im;
        ComplexCP m = new ComplexCP();
        for(int i = 1 ; i <= n ; i++) {
            m.setToQuotientOf(i * _u, 0.0, d[i - 1]._re, d[i - 1]._im);
            d[i]._re = d[0]._re + i * (_v - _u) - m._re * _v * (n - i + 1);
            d[i]._im = d[0]._im - m._im * _v * (n - i + 1);
            t.setToProductOf(m._re, m._im, e[i - 1]._re, e[i - 1]._im);
            e[i]._re = y[i]._re - t._re;
            e[i]._im = y[i]._im - t._im;
        }

        x[n].setToQuotientOf(e[n]._re, e[n]._im, d[n]._re, d[n]._im);
        for(int i = n - 1 ; i >= 0 ; i--) {
            x[i].setToQuotientOf(e[i]._re - x[i + 1]._re * (_v * (n - i)), e[i]._im - x[i + 1]._im * (_v * (n - i)), d[i]._re, d[i]._im);
        }

    }

    private void solve(ComplexCP[] y, ComplexCP offset, ComplexCP[] x, int mx) {

        for(int i = 0 ; i < x.length ; i++) {
            x[i]._re = 0.0;
            x[i]._im = 0.0;
        }
        ComplexCP x1[] = new ComplexCP[mx + 1];
        ComplexCP y1[] = new ComplexCP[mx + 1];
        for(int i = 0 ; i < mx + 1 ; i++){
            x1[i] = new ComplexCP();
            y1[i] = new ComplexCP();
        }

        int xIndex = x.length - 1 - mx;
        int yIndex = y.length - 1 - mx;

        for(int i = 0 ; i < y1.length ; i++) {
            y1[i]._re = y[yIndex + i]._re;
            y1[i]._im = y[yIndex + i]._im;
        }

        solveOneBlock(y1, offset, mx, x1);

        for(int i = 0 ; i < x1.length ; i++) {
            x[xIndex + i]._re = x1[i]._re;
            x[xIndex + i]._im = x1[i]._im;
        }

        for(int m = mx - 1 ; m >= 1 ; m--) {
            for(int i = 0 ; i <= m ; i++) {
                y1[i].setToSumOfScaled(x1[i], 1.0 * (m - i) * (m + 1)/ theta, x1[i + 1], 1.0 * i * (m + 1)/ theta);
            }
            xIndex -= m + 1;
            yIndex -= m + 1;
            for(int i = 0 ; i <= m ; i++) {
                y1[i]._re = y[yIndex + i]._re - y1[i]._re;
                y1[i]._im = y[yIndex + i]._im - y1[i]._im;
            }

            solveOneBlock(y1, offset, m, x1);

            for(int i = 0 ; i <= m ; i++) {
                x[xIndex + i]._re = x1[i]._re;
                x[xIndex + i]._im = x1[i]._im;
            }
        }


    }

    public double [] expQTtx(double time, double [] actualCol, int mx) {

        double xcol[] = new double[actualCol.length + 1];
        System.arraycopy(actualCol, 0, xcol, 1, actualCol.length);
        double y1[] = cfexp(time, xcol, mx);

        double y0[] = new double[y1.length - 1];
        System.arraycopy(y1, 1, y0, 0, y0.length);

        return y0;
    }

    // those numbers can be computed from Matlab package:
    // https://www.mathworks.com/matlabcentral/fileexchange/22055-carathéodory-fejér-approximation?s_tid=gn_loc_drop
    private static final int CF_DEG = 12;
    private static final double [] ci_real = {
            0.000818433612497,
            -0.068571505514864,
            1.319411815998137,
            -8.238258033274786,
            18.785982629476070,
            -11.799383335697918,
            -11.799383335697890,
            18.785982629476067,
            -8.238258033274763,
            1.319411815998138,
            -0.068571505514865,
            0.000818433612497};
    private static final double [] ci_imag = {
            0.000581353207069,
            -0.038419074245887,
            0.183523497750480,
            2.796192505614474,
            -20.237292093573895,
            46.411650777279597,
            -46.411650777279569,
            20.237292093573895,
            -2.796192505614448,
            -0.183523497750480,
            0.038419074245888,
            -0.000581353207069};
    private static final double [] zi_real = {
            -6.998688082445778,
            -2.235968223749446,
            0.851707264834878,
            2.917868800307170,
            4.206124506834328,
            4.827493775040721,
            4.827493775040721,
            4.206124506834328,
            2.917868800307170,
            0.851707264834878,
            -2.235968223749446,
            -6.998688082445778};
    private static final double [] zi_imag = {
            -13.995917029301355,
            -11.109296400461870,
            -8.503832905826961,
            -6.017345968518187,
            -3.590920783130140,
            -1.193987999180278,
            1.193987999180278,
            3.590920783130140,
            6.017345968518187,
            8.503832905826961,
            11.109296400461870,
            13.995917029301355};


    // Caratheodory-Fejer approximation
    private double [] cfexp(double time, double [] v, int mx) {
        int n = v.length;

        ComplexCP x[] = new ComplexCP[n];
        ComplexCP y[] = new ComplexCP[n];
        ComplexCP r[] = new ComplexCP[n];

        r[0] = new ComplexCP();
        for(int i = 1 ; i < n ; i++) {
            r[i] = new ComplexCP(v[i]);
        }

        for(int i = 0 ; i < n ; i++) {
            x[i] = new ComplexCP();
            y[i] = new ComplexCP();
        }

        int nsteps = 1;
        double stepsize = time / nsteps;
        for(int step = 0 ; step < nsteps ; step++) {
            for(int i = 1; i < n ; i++) {
                y[i]._re = r[i]._re / stepsize;
                y[i]._im = r[i]._im / stepsize;
                r[i]._re = 0.0;
                r[i]._im = 0.0;
            }

            ComplexCP offset = new ComplexCP();
            for(int deg = 0 ; deg < CF_DEG ; deg++) {
                offset._re = -zi_real[deg] / stepsize;
                offset._im = -zi_imag[deg] / stepsize;

                solve(y, offset, x, mx);
                for(int j = 1 ; j < n ; j++) {
                    r[j].addProductOf(ci_real[deg], ci_imag[deg], x[j]._re, x[j]._im);
                }
            }
        }

        double real[] = new double[n];
        for(int i = 1 ; i < n ; i++)
            real[i] = r[i]._re;
        return real;
    }

    public static void main(String []args) {
        double pi0 = 0.5;
        double pi1 = 1.0 - pi0;
        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[]{pi0, pi1}, new double[]{1.0 / (2.0 * pi0)});

        MatrixQ q = new MatrixQ(BAGTRModel, 1, 0.01);
        double f[] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}; // first element is reserved
        double g[] = q.cfexp(5, f, 2);
        System.out.println(g);
    }

}
