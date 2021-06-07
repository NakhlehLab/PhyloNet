package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.substitution;

/**
 * Created by wendingqiao on 5/3/16.
 */
public class EigenDecomposition {

    // Eigenvalues, eigenvectors, and inverse eigenvectors
    private final double[] Evec;
    private final double[] Ievc;
    private final double[] Eval;
    private final double[] Evali; // imaginary part of eigenvalues

    public EigenDecomposition(double[] evec, double[] ievc, double[] eval) {
        Evec = evec;
        Ievc = ievc;
        Eval = eval;
        Evali = null;
    }

    public EigenDecomposition(double[] evec, double[] ievc, double[] eval, double[] evali) {
        Evec = evec;
        Ievc = ievc;
        Eval = eval;
        Evali = evali;
    }

    public EigenDecomposition copy() {
        double[] evec = Evec.clone();
        double[] ievc = Ievc.clone();
        double[] eval = Eval.clone();
        return new EigenDecomposition(evec, ievc, eval);
    }

    public final double[] getEigenVectors() {
        return Evec;
    }

    public final double[] getInverseEigenVectors() {
        return Ievc;
    }

    public final double[] getEigenValues() {
        return Eval;
    }

    public final double[] getImEigenValues() {
        return Evali;
    }

    public boolean canReturnComplexDiagonalization() {
        return false;
    }

    /**
     * Rescales the eigenvalues.
     */
    public void normalizeEigenValues(double scale) {
        int dim = Eval.length;
        for (int i = 0; i < dim; i++) {
            Eval[i] /= scale;
        }
    }

    public Boolean hasImagEigenvectors() {
        if (Evali == null) return false;
        for (int i = 0; i < Evali.length; i++)
            if (Evali[i] != 0) {
                return true;
            }
        return false;
    }
}