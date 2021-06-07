package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core.StateNode;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.distribution.GeneTreeBrSpeciesNetDistribution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.distribution.GeneTreeBrSpeciesNetDistributionParallel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.distribution.GeneTreeBrSpeciesNetPseudoLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.all.ScaleAll;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.network.dimension.AddReticulation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.network.dimension.AddReticulationToBackbone;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.network.dimension.DeleteReticulation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.network.dimension.DeleteReticulationFromBackbone;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.network.param.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.network.topo.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.TemporalConstraints;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_Rooted;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkWithTheta;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Ultrametric network
 * Created by wendingqiao on 2/25/16.
 */
public class UltrametricNetwork extends StateNode {

    private Network<NetNodeInfo> _network; // not being copied across states
    protected List<UltrametricTree> _geneTrees;
    private List<Tree> _oldGeneTrees = null;
    protected List<TreeEmbedding> _embeddings;
    private List<TreeEmbedding> _oldEmbeddings;
    private List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> _retiEdges = null;

    private int proposes = 0;

    private Map<String, List<String>> _species2alleles;
    private Map<String, String> _alleles2species;

    public double[] _logGeneTreeNetwork = null;
    private double[] _logLtemp = null;

    private double[] _treeSpaceOpWeights;
    private double[] _treeOpWeights;
    private double[] _netOpWeights;
    private Operator[] _operators;

    private int _numThreads = 1;

    public UltrametricNetwork(List<UltrametricTree> gts, Map<String, List<String>> species2alleles) {
        this(null, gts, species2alleles);
    }

    public UltrametricNetwork(String s, List<UltrametricTree> gts) {
        this(s, gts, null);
    }

    public UltrametricNetwork(String s, List<UltrametricTree> gts, Map<String, List<String>> s2a) {
        this._geneTrees = gts;
        this._species2alleles = s2a;
        this._alleles2species = null;
        this._numThreads = Utils._NUM_THREADS;
        if(Utils.ONLY_BACKBONE_OP)
            this._retiEdges = new ArrayList<>();
        if(this._species2alleles != null) {
            this._alleles2species = new HashMap<>();
            for(String key: this._species2alleles.keySet()){
                for(String allele: this._species2alleles.get(key)){
                    this._alleles2species.put(allele, key);
                }
            }
        }
        boolean init = (s == null);
        double popSize = Double.NaN;
        if(!init) {
            if(s.startsWith("[")) {
                popSize = Double.parseDouble(s.substring(1, s.indexOf("]")));
                s = s.substring(s.indexOf("]") + 1);
            }
            this._network = Networks.readNetwork(s); // adopt topology only
            Set<String> taxa = new HashSet<>();
            for(NetNode<NetNodeInfo> leaf : this._network.getLeaves()) {
                taxa.add(leaf.getName());
            }
            if(_species2alleles != null && (!taxa.containsAll(_species2alleles.keySet()) || !_species2alleles.keySet().containsAll(taxa))) {
                System.err.println("The starting network doesn't match the taxaMap");
                System.err.println(_network.toString() + "\n v.s. \n" + Arrays.toString(_species2alleles.keySet().toArray()));
                System.err.println("The starting network will be set to the MDC tree given gene trees");
                init = true;
            }
        }
        if(init) {
            MDCInference_Rooted mdc = new MDCInference_Rooted();
            List<MutableTuple<Tree,Double>> trees = new ArrayList<>();
            for(UltrametricTree t : this._geneTrees) {
                trees.add(new MutableTuple<>(t.getTree(), 1.0));
            }
            Solution sol = (this._alleles2species == null) ?
                    mdc.inferSpeciesTree(trees, false, 1, false, true, -1).get(0) :
                    mdc.inferSpeciesTree(trees, this._alleles2species, false, 1, false, true, -1).get(0);
            Tree startingTree= Trees.generateRandomBinaryResolution(sol._st);
            s = startingTree.toNewick();
            this._network = Networks.readNetwork(s); // adopt topology only
        }
        Map<String, Double> constraints = TemporalConstraints.getTemporalConstraints(gts, _species2alleles, _alleles2species);
        initNetHeights(popSize, constraints);
        setOperators();
        initializeEmbeddings();
        NetBurnin();
    }

    private void setOperators() {
        if(Utils.ONLY_BACKBONE_OP) {
            this._operators = new Operator[]{
                    new ChangePopSize(this),
                    new ScalePopSize(this),
                    new ScaleAll(_geneTrees, this), // TODO by dw20: sometimes this operator perform poorly
                    new ScaleTime(this), new ScaleRootTime(this), new ChangeTime(this),
                    new AddReticulationToBackbone(this),
                    new FlipReticulationOnBackbone(this),
                    new DeleteReticulationFromBackbone(this),
                    new ChangeInheritance(this)
            };
            if (Utils._ESTIMATE_POP_SIZE) {
                this._treeSpaceOpWeights = Utils.getOperationWeights(
                        Utils.Backbone_Net_Tree_Op_Weights, 0, Utils.Backbone_Net_Tree_Op_Weights.length - 1, !Utils._FIX_NET_TOPOLOGY);
                this._treeOpWeights = Utils.getOperationWeights(Utils.Backbone_Net_Tree_Op_Weights, !Utils._FIX_NET_TOPOLOGY);
                this._netOpWeights = Utils.getOperationWeights(Utils.Backbone_Net_Op_Weights, !Utils._FIX_NET_TOPOLOGY);
            } else {
                this._treeSpaceOpWeights = Utils.getOperationWeights(
                        Utils.Backbone_Net_Tree_Op_Weights, 3, Utils.Backbone_Net_Tree_Op_Weights.length - 1, !Utils._FIX_NET_TOPOLOGY);
                this._treeOpWeights = Utils.getOperationWeights(
                        Utils.Backbone_Net_Tree_Op_Weights, 3, Utils.Backbone_Net_Tree_Op_Weights.length, !Utils._FIX_NET_TOPOLOGY);
                this._netOpWeights = Utils.getOperationWeights(
                        Utils.Backbone_Net_Op_Weights, 3, Utils.Backbone_Net_Op_Weights.length, !Utils._FIX_NET_TOPOLOGY);
            }
        } else {

            this._operators = new Operator[]{
                    new ChangePopSize(this),
                    new ScalePopSize(this),
                    new ScaleAll(_geneTrees, this), // TODO by dw20: sometimes this operator perform poorly
                    new ScaleTime(this), new ScaleRootTime(this), new ChangeTime(this),
                    new SlideSubNet(this), new SwapNodes(this), new MoveTail(this),
                    new AddReticulation(this),
                    new FlipReticulation(this), new MoveHead(this),
                    new DeleteReticulation(this),
                    new ChangeInheritance(this)
            };
            if (Utils._ESTIMATE_POP_SIZE) {
                this._treeSpaceOpWeights = Utils.getOperationWeights(
                        Utils.Net_Tree_Op_Weights, 0, Utils.Net_Tree_Op_Weights.length - 1, !Utils._FIX_NET_TOPOLOGY);
                this._treeOpWeights = Utils.getOperationWeights(Utils.Net_Tree_Op_Weights, !Utils._FIX_NET_TOPOLOGY);
                this._netOpWeights = Utils.getOperationWeights(Utils.Net_Op_Weights, !Utils._FIX_NET_TOPOLOGY);
            } else {
                this._treeSpaceOpWeights = Utils.getOperationWeights(
                        Utils.Net_Tree_Op_Weights, 3, Utils.Net_Tree_Op_Weights.length - 1, !Utils._FIX_NET_TOPOLOGY);
                this._treeOpWeights = Utils.getOperationWeights(
                        Utils.Net_Tree_Op_Weights, 3, Utils.Net_Tree_Op_Weights.length, !Utils._FIX_NET_TOPOLOGY);
                this._netOpWeights = Utils.getOperationWeights(
                        Utils.Net_Op_Weights, 3, Utils.Net_Op_Weights.length, !Utils._FIX_NET_TOPOLOGY);
            }
        }
    }

    // used only for debug only
    public UltrametricNetwork(String s) {
        this._network = Networks.readNetwork(s);
        this._geneTrees = null;
        initNetHeights(null);
    }

    // used only for debug only
    public UltrametricNetwork(Network network, List<UltrametricTree> gts) {
        this._network = network;
        this._geneTrees = gts;
        initNetHeights();
        initializeEmbeddings();
    }

    public void setNetwork(String s) {
        this._network = Networks.readNetwork(s);
        initNetHeights();
    }

    private void NetBurnin() {
        if(Utils._START_NET_BURN_IN) {
            int iter = 0;
            double prevLogL = this.logDensity();
            while (iter < 100000) {
                double logHR = this.propose();
                if(this._operator.getName().contains("Add-Reticulation")) {
                    logHR = Utils.INVALID_MOVE;
                }
//              else if(!this._operator.getName().contains("Time") /*&& !this._operator.getName().contains("Pop")*/) {
//                    logHR = Utils.INVALID_MOVE;
//                }
                else if(!isValid()) {
                    logHR = Utils.INVALID_MOVE;
                }
                this.setDirty(true);
                if(logHR != Utils.INVALID_MOVE) {
                    double newLogL = this.logDensity();
                    System.out.println("Proposed: " + this._operator.getName());
                    System.out.println(Networks.getDendroscopeCompatibleString(this._network));
                    if(newLogL > prevLogL) {
                        prevLogL = newLogL;
                        this.accept();
                        System.out.println("Accepted");
                    } else {
                        this.undo();
                        this.reject();
                    }
                } else {
                    this.undo();
                    this.reject();
                }
                iter++;
            }

            System.out.println("Optimized starting network ");
            System.out.println(Networks.getDendroscopeCompatibleString(this._network));
        }
    }

    public Network<NetNodeInfo> getNetwork() {
        return this._network;
    }

    public int getReticulationCount() {
        return this._network.getReticulationCount();
    }

    public int getInternalNodeCount() {
        return 2 * this._network.getReticulationCount() + this._network.getLeafCount() - 1;
    }

    public boolean isUltrametric() {
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(_network)) {
            double height = node.getData().getHeight();
            if(node.isLeaf()) {
                if(height != Utils.DEFAULT_NET_LEAF_HEIGHT) {
                    System.err.println(height + " vs " + Utils.DEFAULT_NET_LEAF_HEIGHT);
                    return false;
                }
            } else {
                for(NetNode<NetNodeInfo> child : node.getChildren()) {
                    double temp = child.getData().getHeight() + child.getParentDistance(node);
                    if(Math.abs(temp - height) > 0.000001) {
                        System.err.println(node.getName() + " - " + height + " vs " + temp);
                        return false;
                    }
                }
            }
        }
        return true;
    }

    /************ State node methods ************/

    @Override
    public double propose() {
        if(Utils._PRE_BURN_IN) {
            this._operator = getOp(_operators, _treeSpaceOpWeights);
        } else if(_network.getReticulationCount() == 0) {
            this._operator = getOp(_operators, _treeOpWeights);
        } else {
            this._operator = getOp(_operators, _netOpWeights);
        }
        proposes++;

        double logHR = this._operator.propose();
        if(logHR != Utils.INVALID_MOVE && !this._operator.getName().equals("Scale-All") && Utils.RESAMPLE_GENE_TREES ) {
            double rand = Randomizer.getRandomDouble();
            if(rand < Utils.RESAMPLE_GENE_TREE_RATE) {
                // experimental!
                logHR += rebuildGeneTrees();
            }
        }

        return logHR;
    }

    @Override
    public void undo() {
        Set<NetNode> s1 = new HashSet<>();
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(_network)) {
            s1.add(node);
        }
        if(this._operator != null) {
            if (this._operator == null) throw new IllegalArgumentException("null operator");
            this._operator.undo();
        }

        Set<NetNode> s2 = new HashSet<>();
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(_network)) {
            s2.add(node);
        }

        //if(!s1.containsAll(s2)) {
        //    throw new RuntimeException("!!!!!");
        //}
    }

    public List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> getRetiEdges() {
        if(Utils.ONLY_BACKBONE_OP) return _retiEdges;
        else return null;
    }

    @Override
    public double logDensity() {
        if(_logGeneTreeNetwork == null) {
            _logLtemp = computeLikelihood();
            return Utils.sum(_logLtemp);
        }
        // network changed
        if(_dirty || true) {
            _logLtemp = computeLikelihood();
            return Utils.sum(_logLtemp);
        }
        // see if gene tree changed
        else {
            _logLtemp = Utils.copy(_logGeneTreeNetwork);
            GeneTreeBrSpeciesNetDistribution likelihoodCal = new GeneTreeBrSpeciesNetDistribution(_network, _species2alleles);
            for(int i = 0; i < _geneTrees.size(); i++) {
                if(_geneTrees.get(i).isDirty()) {
                    if(Utils.SAMPLE_EMBEDDINGS) {
                        _logLtemp[i] = likelihoodCal.calculateGTDistribution(_geneTrees.get(i), _embeddings.get(i));
                    } else {
                        _logLtemp[i] = likelihoodCal.calculateGTDistribution(_geneTrees.get(i), null);
                    }
                }
            }
        }
        return Utils.sum(_logLtemp);
    }

    public void checkLogDensity() {
        if(_logGeneTreeNetwork != null) {
            double[] l1s = logDensity1();
            double l1 = Utils.sum(l1s);
            double[] l2s = _logGeneTreeNetwork;
            double l2 = Utils.sum(_logGeneTreeNetwork);
            if (Math.abs(l1 - l2) > 0.001) {
                for(int i = 0 ; i < l1s.length ; i++) {
                    if(Math.abs(l1s[i] - l2s[i]) > 0.001) {
                        System.out.print(i + " " + Math.abs(l1s[i] - l2s[i]) + " ");
                    }
                }
                System.out.println();
                throw new RuntimeException("!!!!! " + proposes);
            }
        }
    }

    public double[] logDensity1() {
        return computeLikelihood();
    }

    @Override
    public boolean mayViolate() {
        return _operator.mayViolate();
    }

    @Override
    public void accept() {
        if(Utils.SAMPLE_EMBEDDINGS) {
            for (int i = 0; i < _geneTrees.size(); i++) {
                _oldEmbeddings.get(i).clear();
            }
        }

        if(Utils.RESAMPLE_GENE_TREES) {
            _oldGeneTrees = null;
        }

        _dirty = false;
        if(_logLtemp != null) _logGeneTreeNetwork = _logLtemp;
        _logLtemp = null;

        if(Utils.DEBUG_MODE) {
            checkLogDensity();
        }
    }

    @Override
    public void reject() {
        if(Utils.SAMPLE_EMBEDDINGS) {
            for (int i = 0; i < _geneTrees.size(); i++) {
                if (_oldEmbeddings.get(i).initialized()) {
                    _embeddings.get(i).setTo(_oldEmbeddings.get(i));
                    _oldEmbeddings.get(i).clear();
                }
            }
        }

        if(Utils.RESAMPLE_GENE_TREES) {
            if(_oldGeneTrees != null) {
                System.out.println("Undo rebuilding gene trees");
                for (int i = 0; i < _geneTrees.size(); i++) {
                    _geneTrees.get(i).resetTree(_oldGeneTrees.get(i));
                }
                _oldGeneTrees = null;
            }
        }

        if(Utils.DEBUG_MODE) {
            checkLogDensity();
        }

        _dirty = false;
        _logLtemp = null;
    }

    @Override
    public boolean isValid() {
        if(_network.getRoot().getData().getHeight() > Utils.NET_MAX_HEIGHT) return false;

        Map<String, Double> constraints = TemporalConstraints.getTemporalConstraints(_geneTrees, _species2alleles, _alleles2species);
        Map<String, Double> lowerBound = TemporalConstraints.getLowerBounds(_network);
        for(String key : lowerBound.keySet()) {
            if(constraints.get(key) < lowerBound.get(key)) {
                return false;
            }
        }
        return true;
    }

    /************ moves *************/

    public double[] getLowerAndUpperBound(NetNode<NetNodeInfo> node) {
        double[] bounds = new double[] {Double.MIN_VALUE, Double.MAX_VALUE};
        for(NetNode<NetNodeInfo> child: node.getChildren()) {
            bounds[0] = Math.max(bounds[0], child.getData().getHeight());
        }
        for(NetNode<NetNodeInfo> par: node.getParents()) {
            bounds[1] = Math.min(bounds[1], par.getData().getHeight());
        }
        return bounds;
    }

    // for debug only
    public List<Double> getOldHeights() {
        List<Double> heights = new ArrayList<>();
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(_network)) {
            heights.add(node.getData().getHeight());
        }
        return heights;
    }

    /************ Likelihood computation **************/

    private double[] computeLikelihood() {
        //System.out.println(_network.getReticulationCount());

        double[] likelihoodArray;

        if(!Utils.PSEUDO_LIKELIHOOD) {
            // Compute full likelihood
            GeneTreeBrSpeciesNetDistribution likelihoodCal = new GeneTreeBrSpeciesNetDistributionParallel(
                    _network, _geneTrees, _embeddings, _species2alleles);
            likelihoodArray = new double[_geneTrees.size()];

            Thread[] myThreads = new Thread[_numThreads];
            for (int i = 0; i < _numThreads; i++) {
                myThreads[i] = new MyThreadFromScratch(likelihoodCal, likelihoodArray);
                myThreads[i].start();
            }
            for (int i = 0; i < _numThreads; i++) {
                try {
                    myThreads[i].join();
                } catch (InterruptedException ex) {
                    ex.printStackTrace();
                }
            }
        } else {
            // Compute pseudo-likelihood (experimental)

            likelihoodArray = new double[1];

            GeneTreeBrSpeciesNetPseudoLikelihood likelihood = new GeneTreeBrSpeciesNetPseudoLikelihood();
            likelihoodArray[0] = likelihood.calculateGTDistribution(this, _geneTrees);
        }
        return likelihoodArray;
    }

    private class MyThreadFromScratch extends Thread{

        private GeneTreeBrSpeciesNetDistributionParallel _gtp;
        private double[] _ls;

        public MyThreadFromScratch(GeneTreeBrSpeciesNetDistribution gtp, double[] ls){
            _gtp = (GeneTreeBrSpeciesNetDistributionParallel) gtp;
            _ls = ls;
        }

        public void run() {
            _gtp.calculateGTDistribution(_ls);
        }
    }

    /************** init nodes **************/
    private void initNetHeights(double popSize, Map<String, Double> constraints) {
        if(Double.isNaN(popSize)) {
            initNetHeights(constraints);
            return;
        }
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(_network)) {
            if(node.isLeaf()) {
                node.setData(new NetNodeInfo(Utils.DEFAULT_NET_LEAF_HEIGHT));
                continue;
            }
            double height = Double.MAX_VALUE;
            for(NetNode<NetNodeInfo> child : node.getChildren()) {
                height = Math.min(height, child.getParentDistance(node) + child.getData().getHeight());
            }
            node.setData(new NetNodeInfo(height));
        }

        boolean setPopSize = Utils.varyPopSizeAcrossBranches();
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(_network)) {
            for(NetNode<NetNodeInfo> par : node.getParents()) {
                node.setParentDistance(par, par.getData().getHeight() - node.getData().getHeight());
                if(node.getParentSupport(par) == node.NO_POP_SIZE && setPopSize) {
                    node.setParentSupport(par, popSize);
                }
            }
        }
        _network.getRoot().setRootPopSize(popSize);
    }

    private void initNetHeights(Map<String, Double> constraints) {
        Map<NetNode<NetNodeInfo>, Set<String>> restrictedNodes = TemporalConstraints.getNodeRestriction(this._network);
        Stack<Tuple<NetNode<NetNodeInfo>, Double>> stack = new Stack<>();
        stack.add(new Tuple<>(_network.getRoot(), Utils.DEFAULT_NET_ROOT_HEIGHT));
        while(!stack.isEmpty()) {
            Tuple<NetNode<NetNodeInfo>, Double> tup = stack.pop();
            double height = tup.Item2;
            if(restrictedNodes != null && restrictedNodes.containsKey(tup.Item1) && constraints != null) {
                if(tup.Item1.isNetworkNode()) throw new RuntimeException("Invalid netwok node");
                for(String key : restrictedNodes.get(tup.Item1)) {
                    height = Math.min(height, constraints.get(key));
                }
            }
            height *= Utils.NET_INTI_SCALE;
            if(tup.Item1.getData()!=null) {
                tup.Item1.setData(new NetNodeInfo(Math.min(tup.Item1.getData().getHeight(), height)));
            } else {
                tup.Item1.setData(new NetNodeInfo(height));
            }
            for(NetNode<NetNodeInfo> child : tup.Item1.getChildren()) {
                if(child.isLeaf()) {
                    child.setData(new NetNodeInfo(Utils.DEFAULT_NET_LEAF_HEIGHT));
                } else {
                    stack.add(new Tuple<>(child, height));
                }
            }
        }
        boolean setPopSize = Utils.varyPopSizeAcrossBranches();
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(_network)) {
            for(NetNode<NetNodeInfo> par : node.getParents()) {
                node.setParentDistance(par, par.getData().getHeight() - node.getData().getHeight());
                if(node.isNetworkNode()) {
                    node.setParentProbability(par, 0.5);
                }
                if(node.getParentSupport(par) == node.NO_POP_SIZE && setPopSize) {
                    node.setParentSupport(par, Utils._POP_SIZE_MEAN);
                }
            }
        }
        _network.getRoot().setRootPopSize(Utils._POP_SIZE_MEAN);
    }

    // used only for debug only
    private void initNetHeights() {
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(_network)) {
            if(node.getData() == null) {
                node.setData(new NetNodeInfo(0.0));
            }
            for(NetNode<NetNodeInfo> par : node.getParents()) {
                double dist = node.getParentDistance(par);
                if(par.getData() == null) {
                    par.setData(new NetNodeInfo(node.getData().getHeight() + dist));
                }
            }
        }
    }

    private static void updateWeights(double[] arr, int start) {
        double cutoff = arr[start-1];
        for(int i = 0; i < arr.length; i++) {
            arr[i] = (i < start)? 0 : (arr[i] - cutoff) / (1.0 - cutoff);
        }
        if(Math.abs(arr[arr.length-1] - 1.0) > 0.0000001) {
            throw new IllegalArgumentException(Arrays.toString(arr));
        }
    }

    public void initializeEmbeddings() {
        if(Utils.SAMPLE_EMBEDDINGS) {
            _oldEmbeddings = new ArrayList<>();
            _embeddings = new ArrayList<>();
            for (int i = 0; i < _geneTrees.size(); i++) {
                _oldEmbeddings.add(new TreeEmbedding());
                _embeddings.add(new TreeEmbedding());
            }

            for (int i = 0; i < _geneTrees.size(); i++) {
                EmbeddingRebuilder rebuilder = new EmbeddingRebuilder();
                double temp = rebuilder.rebuild(_geneTrees.get(i), _embeddings.get(i));
                if (temp == Utils.INVALID_MOVE) {
                    i--;
                }
            }
        }
    }

    public double rebuildGeneTrees() {
        double geneTreeLogHR = 0.0;

        double oldLogLikelihood = Utils.sum(_logGeneTreeNetwork);

        SimGTInNetworkWithTheta simulator = new SimGTInNetworkWithTheta();
        simulator.setSeed((long) Randomizer.getRandomInt(99999999));
        int numGT = _geneTrees.size();
        Network<NetNodeInfo> net = _network.clone();

        for(NetNode<NetNodeInfo> node : net.dfs()){
            for(NetNode<NetNodeInfo> parent : node.getParents()) {
                if(node.getParentSupport(parent) == NetNode.NO_SUPPORT) {
                    node.setParentSupport(parent, net.getRoot().getRootPopSize());
                }
            }
        }

        List<Tree> simulatedGTs = new ArrayList<>();
        for(int i = 0 ; i < numGT ; i++) {
            List<Tree> onegt = simulator.generateGTs(net, _species2alleles, 1.0, 1);
            Trees.autoLabelNodes((STITree)onegt.get(0));
            simulatedGTs.add(onegt.get(0));
        }

        _oldGeneTrees = new ArrayList<>();
        double l1s = 0.0;
        double l2s = 0.0;
        for(int i = 0 ; i < numGT ; i++) {
            double l1 = _geneTrees.get(i).logDensity();
            l1s += l1;
            _oldGeneTrees.add(_geneTrees.get(i).getTree());
            _geneTrees.get(i).resetTree(simulatedGTs.get(i));
            //_geneTrees.get(i).resetTree(_oldGeneTrees.get(i));
            double l2 = _geneTrees.get(i).logDensity();
            l2s += l2;
            //if(Math.abs(l1 - l2) > 0.001)
            //    System.out.print("");
            System.out.print("");
        }

        double newLogLikelihood = Utils.sum(computeLikelihood());

        System.out.println("Rebuild gene trees: old " + oldLogLikelihood + " new " + newLogLikelihood);
        System.out.println("old Log(L(gts)) " + l1s + " new Log(L(gts)) " + l2s);

        return oldLogLikelihood - newLogLikelihood;
    }

    public double rebuildEmbeddings() {
        double embeddingLogHR = 0.0;

        if(Utils.SAMPLE_EMBEDDINGS) {
            for (int i = 0; i < _geneTrees.size(); i++) {
                _oldEmbeddings.get(i).setTo(_embeddings.get(i));
                if (_dirty || _geneTrees.get(i).isDirty() /*|| true*/) {
                    EmbeddingRebuilder rebuilder = new EmbeddingRebuilder();
                    double temp = rebuilder.rebuild(_geneTrees.get(i), _embeddings.get(i));
                    if (temp == Utils.INVALID_MOVE) {
                        return Utils.INVALID_MOVE;
                    }
                    embeddingLogHR += temp;
                }
            }
        }
        return embeddingLogHR;
    }

    class EmbeddingRebuilder {
        private Multimap<TNode, TNode> geneTreeNodeHeirs;
        private Multimap<NetNode, TNode> networkNodeHeirs;

        private String gtTaxa[];
        private Map<TNode, STITreeCluster> clusterMap;

        private String getSpecies(String allele) {
            if(_alleles2species == null) return allele;
            return _alleles2species.get(allele);
        }

        private double rebuild(UltrametricTree tree, TreeEmbedding embedding) {
            gtTaxa = tree.getTree().getLeaves();
            clusterMap = new HashMap<>();
            for (TNode node : tree.getTree().postTraverse()) {
                STITreeCluster cl = new STITreeCluster(gtTaxa);
                if (node.isLeaf()) {
                    cl.addLeaf(node.getName());
                } else {
                    for(TNode child : node.getChildren()) {
                        cl = cl.merge(clusterMap.get(child));
                    }
                }
                clusterMap.put(node, cl);
            }

            geneTreeNodeHeirs = HashMultimap.create();
            networkNodeHeirs = HashMultimap.create();
            for (TNode node : tree.getTree().postTraverse()) {
                if (node.isLeaf()) {
                    geneTreeNodeHeirs.put(node, node);
                    networkNodeHeirs.put(_network.findNode(getSpecies(node.getName())), node);
                } else {
                    for(TNode child : node.getChildren()) {
                        geneTreeNodeHeirs.putAll(node, geneTreeNodeHeirs.get(child));
                    }
                }
            }

            for(NetNode<NetNodeInfo> node : Networks.postTraversal(_network)) {
                for(NetNode<NetNodeInfo> child : node.getChildren()) {
                    networkNodeHeirs.putAll(node, networkNodeHeirs.get(child));
                }
            }

            double logHR = 0.0;
            logHR += Math.log(embedding.prob) - Math.log(embedding.probsum);

            TreeEmbedding result = recursivelyRebuild(_network.getRoot(), tree.getTree().getRoot());
            if(result == null) {
                return Utils.INVALID_MOVE;
            }
            embedding.setTo(result);

            logHR -= Math.log(embedding.prob) - Math.log(embedding.probsum);

            return logHR;
        }

        private TreeEmbedding recursivelyRebuild(NetNode<NetNodeInfo> netnode, TNode treenode) {
            if(treenode.getNodeHeight() < netnode.getData().getHeight()) {
                TreeEmbedding[] altEmbeddings = new TreeEmbedding[2];
                Collection<TNode> requiredHeirs = geneTreeNodeHeirs.get(treenode);
                int i = 0;
                double probsum = 0.0;
                for(NetNode<NetNodeInfo> child : netnode.getChildren()) {
                    if(networkNodeHeirs.get(child).containsAll(requiredHeirs)) {
                        altEmbeddings[i] = recursivelyRebuild(child, treenode);
                        if(altEmbeddings[i] == null) return null;
                        altEmbeddings[i].setDirection(treenode, netnode, child);

                        if(child.isNetworkNode()) {
                            double gamma = child.getParentProbability(netnode);
                            altEmbeddings[i].prob *= gamma;
                            altEmbeddings[i].probsum *= gamma;
                        }

                        probsum += altEmbeddings[i].probsum;
                        i++;
                    }
                }
                if(i == 0 || probsum == 0.0)
                    return null;

                double u = Randomizer.getRandomDouble() * probsum;
                if(u < altEmbeddings[0].probsum) {
                    altEmbeddings[0].probsum = probsum;
                    return altEmbeddings[0];
                } else {
                    altEmbeddings[1].probsum = probsum;
                    return altEmbeddings[1];
                }
            } else if(treenode.isLeaf()) {
                return new TreeEmbedding(gtTaxa, clusterMap);
            } else {
                TreeEmbedding embedding = null;
                for(TNode childTreeNode : treenode.getChildren()) {
                    TreeEmbedding childEmbedding = recursivelyRebuild(netnode, childTreeNode);
                    if(childEmbedding == null) {
                        return null;
                    } else if(embedding == null) {
                        embedding = childEmbedding;
                    } else {
                        embedding.merge(childEmbedding);
                    }
                }
                return embedding;
            }
        }
    }
}
