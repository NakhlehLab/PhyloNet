package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.io.File;
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
    public static String _OUT_DIRECTORY = System.getProperty("user.dir");
    // MC3
    public static List<Double> _MC3_CHAINS = null;
    // inference
    public static int _NET_MAX_RETI = 4;
    public static Map<String, List<String>> _TAXON_MAP = null;
    // pop size
    public static boolean _ESTIMATE_POP_SIZE = true;
    public static boolean _CONST_POP_SIZE = true;
    public static boolean _ESTIMATE_POP_PARAM = true;
    // priors
    public static double _POISSON_PARAM = 1.0;
    public static boolean _TIMES_EXP_PRIOR = false;
    public static boolean _DIAMETER_PRIOR = true;
    public static boolean _DISABLE_ALL_PRIOR = false;
    public static final double _DIRICHLET_ALPHA = 2.0; //todo: test dirichlet distribution prior for delta exchange
    public static int NUM_LOCI = 0;
    // Substitution model
    public static String _SUBSTITUTION_MODEL = "JC";
    public static double[] _BASE_FREQS = null;
    public static double[] _TRANS_RATES = null;
    // starting state
    public static boolean _FIX_NET = false;
    public static boolean _FIX_NET_TOPOLOGY = false;
    public static boolean _FIX_GENE_TREE_TOPOLOGIES = false;
    public static boolean _FIX_GENE_TREES = false;
    public static boolean _MLE_NETBL = false;
    public static double _POP_SIZE_MEAN = 0.02;
    public static List<String> _START_NET = null;
    public static List<Map<String, String>> _START_GT_LIST = null;
    public static boolean _START_GT_BURN_IN = false;
    public static boolean _START_NET_BURN_IN = false;
    public static boolean _PRE_BURN_IN = false;
    public static int _PRE_BURN_IN_ITER = 0;
    // summary
    public static int _TOPK_NETS = 10;
    // site model
    public static double _MUTATION_RATE = 1.0;
    // diploid phasing
    public static Set<String> _DIPLOID_SPECIES = null;
    public static boolean _PHASING = false;
    // divergence time window size
    public static double _TIME_WINDOW_SIZE = 0.006;
    public static double _POP_SIZE_WINDOW_SIZE = 0.006;

    // --- net ---
    public static final double NET_INTI_SCALE = 0.95;
    public static final double DEFAULT_NET_LEAF_HEIGHT = 0;
    public static final double DEFAULT_NET_ROOT_HEIGHT = 6;
    public static final double NET_MAX_HEIGHT = 1000; // used in debug mode
    public static double RESAMPLE_GENE_TREE_RATE = 0.5;
    public static boolean RESAMPLE_GENE_TREES = false;
    public static boolean SAMPLE_EMBEDDINGS = false; // experimental!!!
    public static boolean ONLY_BACKBONE_OP = false;
    public static boolean PSEUDO_LIKELIHOOD = false;
    public static boolean SAMPLE_MUTATION_RATE = false; //todo by zhen
    public static int NUM_OPERATORS = 14;
    // --- tree ---
    public static final double TREE_INTI_SCALE = 1.05;
    public static final double DEFAULT_TREE_LEAF_HEIGHT = 0;
    public static final double ROOT_TIME_UPPER_BOUND = 10;
    // --- moves ---
    public static final double INVALID_MOVE = Double.NEGATIVE_INFINITY;
    public static enum  MOVE_TYPE {TREE, NETWORK, ALL, PRIOR};
    public static final double TARGET_ACRATE = 0.345;
    public static enum Transform {None, Log, Sqrt};
    // --- delta exchange operator ---
    public static List<Integer> _WEIGHT_VECTOR = new ArrayList<>();
    public static List<Double> _PARAMETER_INPUT = new ArrayList<>();
    public static boolean _IS_INTEGER_OPERATOR = false;
    public static double _DELTA = 1.0;
    public static double _LOWERBOUND = 0.0;
    public static boolean _MUTATION_RATE_PRIOR = false;
    public static boolean MUTATION_RATE_ONLY = false; //for debug only
    // --- MCMC chain ---
    public static final int SWAP_FREQUENCY = 100;
    // --- priors ---
    public static final double EXP_PARAM = 10; // Mr.Bayes
    public static final double GAMMA_SHAPE = 2; // *BEAST
    // --- substitution model ---
    public static final boolean ESTIMATE_SUBSTITUTION = false; // TODO future improvement

    // --- resume running a chain ---
    public static boolean _PreRun;
    public static String _START_STATE = "";
    public static long _START_TIME;
    public static String _CGT_DIRECTORY = "";
    public static int _SubrunIndex;
    public static int _SubrunNum;
    public static int _SAMPLE_NUM;
    public static boolean _TRUE_START;
    public static int resourceNum;
    public static String _IN_DIRECTORY;
    // --- samples ---
    public static enum SampleType {Tree, Network, ArrayParam, DoubleParam};
    // --- move weights ---
    public static final double DIMENSION_CHANGE_WEIGHT = 0.015;
    public static final double[] Tree_Op_Weights = new double[] {
            0.4, 0.2, 0.2, 0.05, 0.05, 0.05, 0.05
    };
    //todo by zhen: just for temp use, need to test
//    public static final double[] Tree_Op_Weights = new double[] {
//            0.4, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05, 0.1
//    };
    // ChangePopSize ScalePopSize --// ScaleAll
    // ScaleTime ScaleRootTime ChangeTime
    // SlideSubNet SwapNodes MoveTail AddReticulation
    // FlipReticulation MoveHead DeleteReticulation
    // ChangeInheritance
    public static final double[] Murate_Net_Op_Weights = new double[] {
            0.03, 0.01,
            0.01, // scaleAll TODO by dw20: sometimes this operator perform poorly
            0.04, 0.05, 0.25,
            0.20, 0.03, 0.05, DIMENSION_CHANGE_WEIGHT,
            0.1, // deltaexchange TODO by zhen, check mixing
            0.06 - DIMENSION_CHANGE_WEIGHT, 0.05, DIMENSION_CHANGE_WEIGHT,
            0.06 - DIMENSION_CHANGE_WEIGHT

    };
    public static final double[] Murate_Net_Tree_Op_Weights = new double[] {
            0.03, 0.01,
            0.01, // scaleAll TODO by dw20: sometimes this operator perform poorly
            0.04, 0.05, 0.30,
            0.27, 0.06, 0.07 - DIMENSION_CHANGE_WEIGHT * 2, DIMENSION_CHANGE_WEIGHT * 2,
            0.1 // deltaexchange TODO by zhen, check mixing
    };

    public static final double[] Net_Op_Weights = new double[] {
            0.03, 0.01,
            0.01, // scaleAll TODO by dw20: sometimes this operator perform poorly
            0.04, 0.05, 0.25,
            0.20, 0.03, 0.05, DIMENSION_CHANGE_WEIGHT,
            0.06 - DIMENSION_CHANGE_WEIGHT, 0.05, DIMENSION_CHANGE_WEIGHT,
            0.06 - DIMENSION_CHANGE_WEIGHT
    };
    public static final double[] Net_Tree_Op_Weights = new double[] {
            0.03, 0.01,
            0.01, // scaleAll TODO by dw20: sometimes this operator perform poorly
            0.04, 0.05, 0.30,
            0.27, 0.06, 0.07 - DIMENSION_CHANGE_WEIGHT * 2, DIMENSION_CHANGE_WEIGHT * 2
    };
    // ChangePopSize ScalePopSize ScaleAll
    // ScaleTime ScaleRootTime ChangeTime
    //  AddReticulation
    // FlipReticulation DeleteReticulation
    // ChangeInheritance
    public static final double BACKBONE_DIMENSION_CHANGE_WEIGHT = 0.005;
    public static final double[] Backbone_Net_Op_Weights = new double[] {
            0.03, 0.01,
            0.01, // scaleAll TODO by dw20: sometimes this operator perform poorly
            0.04, 0.05, 0.25,
            BACKBONE_DIMENSION_CHANGE_WEIGHT,
            0.06 - BACKBONE_DIMENSION_CHANGE_WEIGHT, BACKBONE_DIMENSION_CHANGE_WEIGHT,
            0.06 - BACKBONE_DIMENSION_CHANGE_WEIGHT
    };
    public static final double[] Backbone_Net_Tree_Op_Weights = new double[] {
            0.03, 0.01,
            0.01, // scaleAll TODO by dw20: sometimes this operator perform poorly
            0.04, 0.05, 0.30,
            BACKBONE_DIMENSION_CHANGE_WEIGHT
    };

    public static final double[] PopSize_Op_Weights = new double[]{0.5, 1.0};
    // phasing
    private static Map<Character, String[]> PHASING_NUCLEOTIDES = null;


    public static double[] getOperationWeights(double[] weights, int start, int end, boolean enableTopologyMoves) {
        if(!enableTopologyMoves) {
            for(int i = 6 ; i < Math.min(end, NUM_OPERATORS-1) ; i++) {
                weights[i] = 0.0;
            }
        }

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

    public static double[] getOperationWeights(double[] weights, boolean enableNetTopologyMoves) {
        if(!enableNetTopologyMoves) {
            for(int i = 6 ; i < Math.min(weights.length, NUM_OPERATORS-1) ; i++) {
                weights[i] = 0.0;
            }
        }

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

    public static void makeUltrametric(STITree t) {
        Map<STINode, Double> heights = new HashMap<>();
        for(Object n : t.postTraverse()) {
            STINode node = (STINode) n;
            if(node.isLeaf()) {
                heights.put(node, 0.0);
            } else {
                double h = 0;
                for(Object c : node.getChildren()) {
                    STINode child = (STINode) c;
                    h = Math.max(h, heights.get(child) + child.getParentDistance());
                }
                for(Object c : node.getChildren()) {
                    STINode child = (STINode) c;
                    child.setParentDistance(h - heights.get(child));
                }
                heights.put(node, h);
            }
        }
    }

    public static boolean isExistsDir(String filePath) {
        File f = new File(filePath);
        if (!f.exists()) {
            f.mkdirs();
        }
        return f.exists() && f.isDirectory();
    }
}
