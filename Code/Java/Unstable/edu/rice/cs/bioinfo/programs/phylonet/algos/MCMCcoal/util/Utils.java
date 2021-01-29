package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.alignment.Alignment;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.List;
import java.util.Set;

/**
 * Utils fields and methods for the whole MCMCcoal program
 * Created by Xinhao Liu on 11/4/19
 * Same usage as edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils
 */
public class Utils {
    public static final long _SEED = 12345678;

    // --- simulation ---
    public static int N0 = 10000;  // N0 for ms
    public static int sequenceLength = -1; // simulated sequence length for ms, be modified according to -r
    public static int CROSS_OVER_RATE = 1000; // -r parameter of ms
    public static final int hmmSmoothingParam = 3; // smoothing on transition probability of hmm hidden states 3

    // --- likelihood calculation ---
    public static double MUTATION_RATE = 0.001; // this is HCG mutation rate of 2.5e-8/site/generation (normal rate)
    //public static double MUTATION_RATE = 0.01; // this is HCG mutation rate of 2.5e-7 (one order of magnitude larger)
//    public static double MUTATION_RATE = 0.1; // this is HCG mutation rate of 2.5e-6 (two orders of magnitude larger)
//    public static double MUTATION_RATE = 0.02; // this is HCG mutation rate of 3.0e-7 (20 folds)
    //public static double MUTATION_RATE = 0.05; // this is HCG mutation rate of 1.25e-6 (50 folds)
    //public static double MUTATION_RATE = 0.00094; // this is HCG mutation rate of 0.094% change per million years (from 09 coalhmm paper)
//    public static double MUTATION_RATE = 0.0116; // this is butterfly mutation rate of 2.9e-9/site/generation, assuming pop size 1000000
//    public static final double MUTATION_RATE = 0.04; // this is butterfly mutation rate of 1e-8/site/generation, assuming pop size 1000000
//    public static final double MUTATION_RATE = 0.116; // this is butterfly mutation rate of 2.9e-8/site/generation, assuming pop size 1000000
    public static Alignment DATA = null;
    public static double ILLEGAL_LIKELIHOOD = 0; // likelihood for models with illegal parameter values
    public static int NUM_BIN = 2;

    // --- priors ---
    public static final boolean DISABLE_ALL_PRIOR = false;
    public static final double GAMMA_SHAPE = 2; // *BEAST, MCMC_SEQ

    // --- tree ---
    public static final int DEFAULT_TREE_ROOT_HEIGHT = 720000; // This is human-orangutan divergence time in unit of generations (18 Mya, 25 y generation time)
    public static final int DEFAULT_TREE_LEAF_HEIGHT = 0;
    public static final double TREE_INTI_SCALE = 0.5;
    public static int POP_SIZE_MEAN = 50000;

    // --- variational inference settings, learning settings---
    public static boolean ILLEGAL_SAMPLE_GENERATED = false;
//    public static final double NODE_HEIGHT_INIT_STDDEV = 10000.0;
//    public static final double POP_SIZE_INIT_STDDEV = 5000.0;
    public static double NODE_HEIGHT_INIT_STDDEV = 20000.0;
    public static double POP_SIZE_INIT_STDDEV = 10000.0;
    public static final double RECOMB_RATE_INIT_STDDEV = 0.8;
    public static int nSamples = 50;
    public static int nIterations = 200;
//    public static double NODE_HEIGHT_MEAN_LEARNING_RATE = 20000;
//    public static double NODE_HEIGHT_STDDEV_LEARNING_RATE = 5000; // 1500 too low?
//    public static double POP_SIZE_MEAN_LEARNING_RATE = 10000;
//    public static double POP_SIZE_STDDEV_LEARNING_RATE = 5000;
    public static double NODE_HEIGHT_MEAN_LEARNING_RATE = 20000;
    public static double NODE_HEIGHT_STDDEV_LEARNING_RATE = 500;
    public static double POP_SIZE_MEAN_LEARNING_RATE = 10000;
    public static double POP_SIZE_STDDEV_LEARNING_RATE = 500;
    public static final double RECOMB_RATE_MEAN_LEARNING_RATE = 0.000005;
    public static final double RECOMB_RATE_STDDEV_LEARNING_RATE = 0.000001;
//    public static double RECOMB_RATE_SCALE = 1E-8; // not useful anymore
    public static double SIGMA_LEARNING_RATE = 25000000; // 5000^2
    public static double NODE_HEIGHT_MIN_STDDEV = 10000;
    public static double POP_SIZE_MIN_STDDEV = 3000;

    // --- variational inference settings, learning settings for reparametrized model---
    public static double BRANCH_LENGTH_INIT_STDDEV = 20000;
    public static double BRANCH_LENGTH_MEAN_LEARNING_RATE = 20000;
    public static double BRANCH_LENGTH_STDDEV_LEARNING_RATE = 500;
    public static double BRANCH_LENGTH_MIN_STDDEV = 10000;

    // --- logging ---
    public static long buildingTime = 0;
    public static long likelihoodTime = 0;

    // --- from mcmcseq ---
    public static boolean _PHASING = false;
    public static double EPSILON = 2.220446049250313E-16; // what is it

    // --- testing ---
    public static boolean lastiter = false;

    /**
     * Return the Hadamard product of m1 and m2.
     */
    public static RealMatrix hadamardProduct(RealMatrix m1, RealMatrix m2) {
        RealMatrix result = m1.copy();
        for (int i = 0; i < result.getRowDimension(); i++) {
            for (int j = 0; j < result.getColumnDimension(); j++) {
                result.multiplyEntry(i, j, m2.getEntry(i, j));
            }
        }
        return result;
    }
}
