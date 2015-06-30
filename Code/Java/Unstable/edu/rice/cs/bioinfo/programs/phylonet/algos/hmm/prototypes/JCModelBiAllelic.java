package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.prototypes;

import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.SubstitutionModel;
import jeigen.DenseMatrix;

public class JCModelBiAllelic extends SubstitutionModel {
    DenseMatrix rateMatrix;
    double U;
    double V;
    DenseMatrix rateMatrixIntegrated;

    public JCModelBiAllelic(double u, double v)
    {
        U = u;
        V = v;
        rateMatrix = createMatrix(u, v);
    }

    static DenseMatrix createMatrix(double u, double v)
    {
        double[][] temp = new double[2][2];
        temp[0][0] = -u;
        temp[0][1] = v;
        temp[1][0] = u;
        temp[1][1] = -v;

        return new DenseMatrix(temp);
    }

    @Override
    public DenseMatrix getEquilibriumVector() {
        return new DenseMatrix(new double[][] {{V/(U+V)},{U/(U+V)}});
    }

    @Override
    public DenseMatrix getRateMatrix() {
        return rateMatrix;
    }


    public static void main(String[] args)
    {
        JCModelBiAllelic foo = new JCModelBiAllelic(1,5);
        System.out.println("Test of JCModelBiAlellic");


        System.out.println("Matrix Result: = " + foo.rateMatrix.mmul(foo.getEquilibriumVector()));
        System.out.println("get prb matrix = ");
        System.out.println(foo.getProbabilityMatrix(.5));

        System.out.println("EQUILIBRIUM MATRIX = " + foo.getEquilibriumVector());

    }

    public DenseMatrix getProbabilityMatrixIntegrated(){
        return rateMatrixIntegrated;
    }
}
