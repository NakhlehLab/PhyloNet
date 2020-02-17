package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.alignment.Alignment;

import java.util.List;
import java.util.Set;

/**
 * Utils fields and methods for the whole MCMCcoal program
 * Created by Xinhao Liu on 11/4/19
 * Same usage as edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils
 */
public class Utils {
    public static long _SEED = 12345678;

    // --- simulation ---
    public static final int N0 = 10000;  // N0 for ms
    public static int sequenceLength = 500000; // simulated sequence length for ms TODO: just removed final
    public static final int hmmSmoothingParam = 3; // smoothing on transition probability of hmm hidden states

    // --- likelihood calculation ---
    //public static double MUTATION_RATE = 1.0;
//    public static double MUTATION_RATE = 0.001;
    public static double MUTATION_RATE = 0.01; // TODO: this is HCG mutation rate of 2.5e-7 (one order of magnitude larger)
    public static List<Alignment> DATA = null;

    // --- tree ---
    // TODO: how to set default root height?
    public static final int DEFAULT_TREE_ROOT_HEIGHT = 720000; // This is human-orangutan divergence time in unit of generations (18 Mya, 25 y generation time)
    public static final int DEFAULT_TREE_LEAF_HEIGHT = 0;
    public static final double TREE_INTI_SCALE = 0.5;
    // TODO: how to set initial population size?
    public static int POP_SIZE_MEAN = 20000;

    // --- from mcmcseq ---
    public static boolean _PHASING = false;
    public static double EPSILON = 2.220446049250313E-16; //TODO: what is it

}
