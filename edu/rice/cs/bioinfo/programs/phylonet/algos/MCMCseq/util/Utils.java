package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Utils fields and methods for the whole program
 * Created by wendingqiao on 3/3/16.
 */
public class Utils {

    public static final boolean DEBUG_MODE = false;
    public static double EPSILON = 2.220446049250313E-16;

    // --- settings ---
    // MCMC chain
    public static long _CHAIN_LEN = 10000000;
    public static long _BURNIN_LEN = 2000000;
    public static long _SAMPLE_FREQUENCY = 5000;
    public static long _SEED = 12345678;
    public static int _NUM_THREADS = Runtime.getRuntime().availableProcessors();
    public static String _OUT_DIRECTORY = System.getProperty("user.home");
    // MC3
    public static List<Double> _MC3_CHAINS = null;
    // inference
    public static int _NET_MAX_RETI = 4;
    public static Map<String, List<String>> _TAXON_MAP = null;
    // pop size
    public static boolean _ESTIMATE_POP_SIZE = true;
    public static boolean _CONST_POP_SIZE = true;
    // priors
    public static double _POISSON_PARAM = 1.0;
    public static boolean _TIMES_EXP_PRIOR = false;
    public static boolean _DIAMETER_PRIOR = true;
    // Substitution model
    public static String _SUBSTITUTION_MODEL = "JC";
    public static double[] _BASE_FREQS = null;
    public static double[] _TRANS_RATES = null;
    // starting state
    public static double _POP_SIZE_MEAN = 0.036;
    public static String _START_NET = null;
    public static List<String> _START_GT_LIST = null;
    public static boolean _PRE_BURN_IN = true;
    public static int _PRE_BURN_IN_ITER = 10;
    // summary
    public static int _TOPK_NETS = 10;
    // site model
    public static double _MUTATION_RATE = 1.0;
    // diploid phasing
    public static Set<String> _DIPLOID_SPECIES = null;
    public static boolean _PHASING = false;
    // divergence time window size
    public static double _TIME_WINDOW_SIZE = 0.01;
    public static double _POP_SIZE_WINDOW_SIZE = 0.01;

    // --- net ---
    public static final double NET_INTI_SCALE = 0.95;
    public static final double DEFAULT_NET_LEAF_HEIGHT = 0;
    public static final double DEFAULT_NET_ROOT_HEIGHT = 6;
    public static final double NET_MAX_HEIGHT = 1000; // used in debug mode
    // --- tree ---
    public static final double TREE_INTI_SCALE = 1.05;
    public static final double DEFAULT_TREE_LEAF_HEIGHT = 0;
    public static final double ROOT_TIME_UPPER_BOUND = 1000;
    // --- moves ---
    public static final double INVALID_MOVE = Double.NEGATIVE_INFINITY;
    public static enum  MOVE_TYPE {TREE, NETWORK, ALL, PRIOR};
    public static final double TARGET_ACRATE = 0.345;
    public static enum Transform {None, Log, Sqrt};
    // --- MCMC chain ---
    public static final int SWAP_FREQUENCY = 100;
    // --- priors ---
    public static final double EXP_PARAM = 10; // Mr.Bayes
    public static final double GAMMA_SHAPE = 2; // *BEAST
    // --- substitution model ---
    public static final boolean ESTIMATE_SUBSTITUTION = false; // TODO future improvement
    // --- samples ---
    public static enum SampleType {Tree, Network, ArrayParam, DoubleParam};
    // --- move weights ---
    public static final double DIMENSION_CHANGE_WEIGHT = 0.015;
    public static final double[] Tree_Op_Weights = new double[] {
            0.4, 0.2, 0.2, 0.05, 0.05, 0.05, 0.05
    };
    // ChangePopSize ScalePopSize --// ScaleAll
    // ScaleTime ScaleRootTime ChangeTime
    // SlideSubNet SwapNodes MoveTail AddReticulation
    // FlipReticulation MoveHead DeleteReticulation
    // ChangeInheritance
    public static final double[] Net_Op_Weights = new double[] {
            0.03, 0.01,
            0.01, // scaleAll TODO by dw20: sometimes this operator perform poorly
            0.04, 0.05, 0.25,
            0.20, 0.03, 0.05, DIMENSION_CHANGE_WEIGHT,
            0.06 - DIMENSION_CHANGE_WEIGHT, 0.05, DIMENSION_CHANGE_WEIGHT, 0.06 - DIMENSION_CHANGE_WEIGHT
    };
    public static final double[] Net_Tree_Op_Weights = new double[] {
            0.03, 0.01,
            0.01, // scaleAll TODO by dw20: sometimes this operator perform poorly
            0.04, 0.05, 0.30,
            0.27, 0.06, 0.07 - DIMENSION_CHANGE_WEIGHT * 2, DIMENSION_CHANGE_WEIGHT * 2
    };
    public static final double[] PopSize_Op_Weights = new double[]{0.5, 1.0};
    // phasing
    private static Map<Character, String[]> PHASING_NUCLEOTIDES = null;


    public static double[] getOperationWeights(double[] weights, int start, int end) {
        double[] arr = new double[weights.length];
        double sum = 0;
        for(int i = start; i < end; i++) {
            sum += weights[i];
        }
        for(int i = 0; i < weights.length; i++) {
            if (i < start) {
                arr[i] = 0;
            } else if (i >= end - 1) {
                arr[i] = 1;
            } else {
                arr[i] = weights[i] / sum + (i == 0 ? 0 : arr[i-1]);
            }
        }
        return arr;
    }

    public static double[] getOperationWeights(double[] weights) {
        double[] arr = new double[weights.length];
        double sum = 0;
        for(double d : weights) {
            sum += d;
        }
        for(int i = 0; i < weights.length - 1; i++) {
            arr[i] = weights[i] / sum + (i == 0 ? 0 : arr[i-1]);
        }
        arr[weights.length-1] = 1;
        return arr;
    }

    public static double sum(double[] array) {
        double res = 0.0;
        for(double d : array) res += d;
        return res;
    }

    public static double[] copy(double[] array) {
        double[] res = new double[array.length];
        for(int i = 0; i < array.length; i++) {
            res[i] = array[i];
        }
        return res;
    }

    public static double calcDelta(Operator operator, double logAlpha) {
        double target = Utils.TARGET_ACRATE;

        double count = (operator._rejCorrectionCounter + operator._acCorrectionCounter + 1.0);
        switch (operator._transform) {
            case Log:
                count = Math.log(count + 1.0);
                break;
            case Sqrt:
                count = Math.sqrt(count);
                break;
            case None:
            default:
                break;
        }

        double deltaP = (Math.exp(Math.min(logAlpha, 0)) - target) / count;

        return (deltaP > -Double.MAX_VALUE && deltaP < Double.MAX_VALUE) ? deltaP : 0;
    }

    public static String plotNetwork(String netStr) {
        Network net = Networks.readNetwork(netStr);
        for(Object n : Networks.postTraversal(net)) {
            NetNode node = (NetNode) n;
            for(Object p : node.getParents()) {
                NetNode par = (NetNode) p;
                node.setParentSupport(par, node.NO_SUPPORT);
                node.setParentProbability(par, node.NO_PROBABILITY);
            }
        }
        return net.toString();
    }

    public static boolean equals(double a, double b) {
        return Math.abs(a-b) < EPSILON;
    }

    public static boolean varyPopSizeAcrossBranches() {
        return Utils._ESTIMATE_POP_SIZE && !Utils._CONST_POP_SIZE;
    }

    public static void taxonMapPhasing(Set<String> taxa) {
        if(Utils._TAXON_MAP == null) {
            Utils._TAXON_MAP = new HashMap<>();
            for(String key : taxa) {
                Utils._TAXON_MAP.put(key, Arrays.asList(new String[] {key}));
            }
        }
        Map<String, List<String>> taxonMap = new HashMap<>();
        for(String key : Utils._TAXON_MAP.keySet()) {
            List<String> treeTaxa = new ArrayList<>();
            for(String val : Utils._TAXON_MAP.get(key)) {
                if(Utils._DIPLOID_SPECIES.contains(val)) {
                    treeTaxa.add(val + "_1");
                    treeTaxa.add(val + "_2");
                } else {
                    treeTaxa.add(val);
                }
                taxonMap.put(key, treeTaxa);
            }
        }
        Utils._TAXON_MAP = taxonMap;
    }

    public static Map<Character, String[]> getPhasingNucleotides() {
        if(PHASING_NUCLEOTIDES == null) {
            Map<Character, String[]> nucleotides = new HashMap<>();
            nucleotides.put('R', new String[] {"AG", "GA"});
            nucleotides.put('Y', new String[] {"CT", "TC"});
            nucleotides.put('M', new String[] {"AC", "CA"});
            nucleotides.put('W', new String[] {"AT", "TA"});
            nucleotides.put('S', new String[] {"CG", "GC"});
            nucleotides.put('K', new String[] {"GT", "TG"});
            PHASING_NUCLEOTIDES = nucleotides;
        }
        return PHASING_NUCLEOTIDES;
    }
}