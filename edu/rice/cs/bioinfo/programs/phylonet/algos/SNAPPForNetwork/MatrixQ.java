package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

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

    public MatrixQ(RateModel rModel, int M, double theta){
        this.M = M;
        this.theta = theta;

        rateMatrix = rModel.getRateMatrix();
        rateModel = rModel;
        transposedQ = initializeMatrixQ(rateMatrix,M,theta);
        equilibrium = computeEquilibrium();
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

    public void solveCentralBlock(Complex[] y, Complex offset, int n, double u, double v, double theta, Complex[] x ) {
        Complex K = offset.add(-1.0 * n * (n - 1) / theta - n * v);
        Complex tmp = Complex.ZERO;
        Complex d[] = new Complex[n + 1];
        Complex e[] = new Complex[n + 1];
        d[0] = K;
        e[0] = y[0];
        for(int i = 1 ; i <= n ; i++) {
            tmp = new Complex(i * u);
            Complex m = tmp.divide(d[i - 1]);
            d[i] = K.subtract(m.multiply(v * (n - i + 1))).add(i * (v - u));
            tmp = m.multiply(e[i - 1]);
            e[i] = y[i].subtract(tmp);
        }

        x[n] = e[n].divide(d[n]);
        for(int i = n - 1 ; i >= 0 ; i--) {
            tmp =e[i].subtract(x[i + 1].multiply(v * (n - i)));
            x[i] = tmp.divide(d[i]);
        }

    }

    public void multiplyUpperBlock(Complex[] x, int n, double theta, Complex[] y) {
        for(int i = 0 ; i <= n ; i++) {
            y[i] = x[i].multiply(1.0 * (n - i) * (n + 1)/ theta).add(x[i + 1].multiply(1.0 * i * (n + 1)/ theta));
        }
    }

    public void solve(Complex[] y, Complex offset, Complex[] x, int mx) {
        double u = 0.0;
        double v = 0.0;

        if(rateModel instanceof BiAllelicGTR) {
            BiAllelicGTR ba = (BiAllelicGTR) rateModel;
            u = ba.getU();
            v = ba.getV();
        } else {
            throw new RuntimeException("only for bi-allele");
        }

        for(int i = 0 ; i < x.length ; i++) {
            x[i] = Complex.ZERO;
        }
        Complex xn[] = new Complex[mx + 1];
        Complex yn[] = new Complex[mx + 1];
        for(int i = 0 ; i < mx + 1 ; i++){
            xn[i] = Complex.ZERO;
            yn[i] = Complex.ZERO;
        }

        int xptr = x.length - 1 - mx;
        int yptr = y.length - 1 - mx;

        System.arraycopy(y, yptr, yn, 0, yn.length);

        solveCentralBlock(yn, offset, mx, u, v, theta, xn);

        System.arraycopy(xn, 0, x, xptr, xn.length);

        for(int m = mx - 1 ; m >= 1 ; m--) {
            xptr = xptr - (m + 1);
            multiplyUpperBlock(xn, m, theta, yn);

            yptr = yptr - (m + 1);
            for(int i = 0 ; i <= m ; i++) {
                yn[i] = y[yptr + i].subtract(yn[i]);
            }

            solveCentralBlock(yn, offset, m, u, v, theta, xn);

            System.arraycopy(xn, 0, x, xptr, m + 1);
        }


    }

    public DenseMatrix expQTtx(double time, double [] actualCol, int mx) {
        double u = 0.0;
        double v = 0.0;

        if(rateModel instanceof BiAllelicGTR) {
            BiAllelicGTR ba = (BiAllelicGTR) rateModel;
            u = ba.getU();
            v = ba.getV();
        } else {
            throw new RuntimeException("only for bi-allele");
        }

        double xcol[] = new double[actualCol.length + 1];
        System.arraycopy(actualCol, 0, xcol, 1, actualCol.length);
        double y1[] = cfexp(time, xcol, mx);

        double y0[] = new double[y1.length - 1];
        System.arraycopy(y1, 1, y0, 0, y0.length);

        return new DenseMatrix(new double[][]{y0}).t();
    }

    public static int CF_DEG = 12;
    public static double [] ci_real = {
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
    public static double [] ci_imag = {
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
    public static double [] zi_real = {
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
    public static double [] zi_imag = {
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



    public double [] cfexp(double time, double [] v, int mx) {
        int n = v.length;
        double w[] = new double[n];

        Complex xc[] = new Complex[n];
        Complex wc[] = new Complex[n];
        Complex vc[] = new Complex[n];

        wc[0] = Complex.ZERO;
        for(int i = 1 ; i < n ; i++) {
            wc[i] = new Complex(v[i]);
        }

        for(int i = 0 ; i < n ; i++) {
            xc[i] = Complex.ZERO;
            vc[i] = Complex.ZERO;
        }

        int steps = (int)(Math.ceil(Math.log(n) / Math.log(2)));
        steps = 1;
        double stepsize = time / steps;
        for(int k = 0 ; k < steps ; k++) {
            for(int i = 1; i < n ; i++) {
                vc[i] = wc[i].divide(stepsize);
                wc[i] = Complex.ZERO;
            }

            for(int i = 0 ; i < CF_DEG ; i++) {
                Complex offset = new Complex(-zi_real[i] / stepsize, -zi_imag[i] / stepsize);
                solve(vc, offset, xc, mx);
                Complex ci = new Complex(ci_real[i], ci_imag[i]);
                for(int j = 1 ; j < n ; j++) {
                    wc[j] = wc[j].add(ci.multiply(xc[j]));
                }
            }
        }

        for(int i = 1 ; i < n ; i++)
            w[i] = wc[i].getReal();
        return w;
    }

    public static void testExp() {
        long start = System.currentTimeMillis();

        FMatrix fb = new FMatrix(30, false);
        R r = new R(30, new int[]{5});
        fb.set(r, 1.0);
        //System.out.println(fb);

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {0.5, 0.5}, new double[] {1.0});
        MatrixQ Q = new MatrixQ(BAGTRModel , 30, 0.04);
        FMatrix ft = new FMatrix(30, false);
        for(int i = 0 ; i < 10000 ; i++)
        ft.setMatrix(Q.expQTtx(0.5, fb.getArr(), fb.mx));
        //System.out.println(ft);

        System.out.println("Seconds: " + (System.currentTimeMillis()-start)/1000.0);
    }

    public static void testComplex() {
        long start = System.currentTimeMillis();

        Complex c = new Complex(1.0, 1.0);
        Complex t1 = new Complex(1.0, 0.0);
        Complex t2 = new Complex(0.0, 1.0);
        for(int i = 0 ; i < 10000000 ; i++) {
            c = c.divide(1.0);
        }
        System.out.println(c);
        System.out.println((System.currentTimeMillis()-start)/1000.0);
    }

    public static void main(String[] args){
        testExp();
        //testComplex();
        //R.dims = 1;
//        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {5/6.0, 1/6.0}, new double[] {0.1});
//        System.out.println("The BiAllelicGTR Model is: " + BAGTRModel.getRateMatrix());
//        MatrixQ tester = new MatrixQ(BAGTRModel , 5, 1);
//
//        System.out.println(tester.transposedQ);
    }

}
