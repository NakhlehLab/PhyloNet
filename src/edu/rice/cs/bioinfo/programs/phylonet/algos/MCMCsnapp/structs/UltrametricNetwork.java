package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.core.StateNode;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SNAPPLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SNAPPLikelihoodSampling;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.MarkerSeq;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.all.ScaleAll;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.dimension.AddReticulation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.dimension.DeleteReticulation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.param.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.topo.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.TemporalConstraints;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.RPattern;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_Rooted;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Ultrametric network
 */
public class UltrametricNetwork extends StateNode {

    private Network<NetNodeInfo> _network;
    private List<UltrametricTree> _geneTrees;
    private List<MarkerSeq> _markers; // This list always has size 1.
    private Map<String, List<String>> _species2alleles;
    private Map<String, String> _alleles2species;

    private double[] _logGeneTreeNetwork = null;
    private double[] _logLtemp = null;

    private double[] _treeOpWeights;
    private double[] _netOpWeights;
    private Operator[] _operators;
    private BiAllelicGTR _BAGTRModel;

    private List<Splitting> _splittings;
    private int _numSites = 0;
    private Map<String, Double> _abcData;
    private Map<String, Double> _abcDataPrev;

    private int _numThreads;

    public UltrametricNetwork(String s, Map<String, List<String>> species2alleles) {
        this(s, null, null, species2alleles, null);
    }

    public UltrametricNetwork(List<UltrametricTree> gts, Map<String, List<String>> species2alleles) {
        this(null, gts, null, species2alleles, null);
    }

    public UltrametricNetwork(String s, List<UltrametricTree> gts) {
        this(s, gts, null, null, null);
    }

    public UltrametricNetwork(String s, List<UltrametricTree> gts, List<MarkerSeq> markerSeqs, Map<String, List<String>> s2a, BiAllelicGTR BAGTRModel) {
        this._markers = markerSeqs;
        this._numSites = markerSeqs.get(0).getSiteCount();
        this._geneTrees = gts;
        this._species2alleles = s2a;
        this._alleles2species = null;
        this._BAGTRModel = BAGTRModel;
        this._numThreads = Utils._NUM_THREADS;
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
                //Utils._POP_SIZE_MEAN = popSize;
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

        if(Utils.SAMPLE_SPLITTING) {
            initSplitting();
        } else {
            _splittings = null;
        }

        if(SNAPPLikelihood.useApproximateBayesian) {
            _abcData = new HashMap<>();
        }

        setOperators();
    }

    private void setOperators() {
        if(Utils._MCMC) {
            this._operators = new Operator[]{
                    new ChangePopSize(this),
                    new ScalePopSize(this),
                    new ScaleAll(_geneTrees, this),
                    new ScaleTime(this), new ScaleRootTime(this), new ChangeTime(this),
                    new SlideSubNet(this), new SwapNodes(this), new MoveTail(this),
                    new AddReticulation(this),
                    new FlipReticulation(this), new MoveHead(this),
                    new DeleteReticulation(this),
                    new ChangeInheritance(this)
            };
            if (Utils._ESTIMATE_POP_SIZE) {
                this._treeOpWeights = Utils.getOperationWeights(Utils.Net_Tree_Op_Weights);
                this._netOpWeights = Utils.getOperationWeights(Utils.Net_Op_Weights);
            } else {
                this._treeOpWeights = Utils.getOperationWeights(
                        Utils.Net_Tree_Op_Weights, 3, Utils.Net_Tree_Op_Weights.length);
                this._netOpWeights = Utils.getOperationWeights(
                        Utils.Net_Op_Weights, 3, Utils.Net_Op_Weights.length);
            }
        } else {
            this._operators = new Operator[]{
                    new ChangePopSize(this),
                    new ScalePopSize(this),
                    new ScaleAll(_geneTrees, this),
                    new ScaleTime(this), new ScaleRootTime(this), new ChangeTime(this),
                    new SlideSubNet(this), new SwapNodes(this), new MoveTail(this),
                    new AddReticulation(this),
                    new FlipReticulation(this), new MoveHead(this),
                    new DeleteReticulation(this),
                    new ChangeInheritance(this),
                    new ReplaceReticulation(this)
            };
            if (Utils._ESTIMATE_POP_SIZE) {
                this._treeOpWeights = Utils.getOperationWeights(Utils.Search_Net_Tree_Op_Weights);
                this._netOpWeights = Utils.getOperationWeights(Utils.Search_Net_Op_Weights);
            } else {
                this._treeOpWeights = Utils.getOperationWeights(
                        Utils.Search_Net_Tree_Op_Weights, 3, Utils.Search_Net_Tree_Op_Weights.length);
                this._netOpWeights = Utils.getOperationWeights(
                        Utils.Search_Net_Op_Weights, 3, Utils.Search_Net_Op_Weights.length);
            }
        }

    }

    // used only for debug only
    public UltrametricNetwork(String s, List<MarkerSeq> markerSeqs, Map<String, List<String>> species2alleles, BiAllelicGTR BAGTRModel) {
        this._network = Networks.readNetworkWithRootPop(s);
        this._species2alleles = species2alleles;
        if(this._species2alleles != null) {
            this._alleles2species = new HashMap<>();
            for(String key: this._species2alleles.keySet()){
                for(String allele: this._species2alleles.get(key)){
                    this._alleles2species.put(allele, key);
                }
            }
        }
        this._geneTrees = null;
        this._markers = markerSeqs;
        this._numSites = markerSeqs.get(0).getSiteCount();
        initNetHeights();
        this._BAGTRModel = BAGTRModel;

        if(Utils.SAMPLE_SPLITTING) {
            initSplitting();
        } else {
            _splittings = null;
        }
    }

    public UltrametricNetwork(String s) {
        this._network = Networks.readNetworkWithRootPop(s);
        initNetHeights();

    }

    public Network<NetNodeInfo> getNetwork() {
        return this._network;
    }

    public void setNetwork(Network<NetNodeInfo> network) {
        this._network = network;
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
        double logHR = 0.0;

        if(_network.getReticulationCount() == 0) {
            this._operator = getOp(_operators, _treeOpWeights);
        } else {
            this._operator = getOp(_operators, _netOpWeights);
        }
        logHR = this._operator.propose();

        if(Utils.SAMPLE_SPLITTING) {
            // experimental!
            Networks.autoLabelNodes(_network);
            for(int i = 0 ; i < _splittings.size() ; i++) {
                _splittings.get(i).propose();
            }
        }

        if(SNAPPLikelihood.useApproximateBayesian) {
            _abcDataPrev = new HashMap<>(_abcData) ;
            _abcData = new HashMap<>();
        }

        return logHR;
    }

    @Override
    public void undo() {
        if(this._operator == null) throw new IllegalArgumentException("null operator");
        this._operator.undo();
        if(Utils.SAMPLE_SPLITTING) {
            // experimental!
            for(int i = 0 ; i < _splittings.size() ; i++) {
                _splittings.get(i).undo();
            }
        }

        if(SNAPPLikelihood.useApproximateBayesian) {
            _abcData = new HashMap<>(_abcDataPrev);
            _abcDataPrev.clear();
            _logGeneTreeNetwork = computeLikelihood(); // approximate bayesian
        }
    }

    @Override
    public double logDensity() {
        if(_logGeneTreeNetwork == null) {
            _logLtemp = computeLikelihood();
            return Utils.sum(_logLtemp);
        }
        // network changed
        if(_dirty) {
            if(_network.getReticulationCount() > Utils._NET_MAX_RETI)
                return Utils.INVALID_MOVE;
            else {
                _logLtemp = computeLikelihood();
                return Utils.sum(_logLtemp);
            }
        }
        else {
            _logLtemp = Utils.copy(_logGeneTreeNetwork);
        }

        return Utils.sum(_logLtemp);
    }

    public double recomputeLogDensity() {

        return Utils.sum(computeLikelihood());
    }

    @Override
    public boolean mayViolate() {
        return _operator.mayViolate();
    }

    @Override
    public void accept() {
        _dirty = false;
        if(_logLtemp != null) {
            _logGeneTreeNetwork = _logLtemp;
        }
        _logLtemp = null;

        if(Utils.SAMPLE_SPLITTING) {
            for(int i = 0 ; i < _splittings.size() ; i++) {
                _splittings.get(i).accept();
            }
        }

        if(SNAPPLikelihood.useApproximateBayesian) {
            _abcDataPrev.clear();
        }
    }

    @Override
    public void reject() {
        _dirty = false;
        _logLtemp = null;

        if(Utils.SAMPLE_SPLITTING) {
            for(int i = 0 ; i < _splittings.size() ; i++) {
                _splittings.get(i).reject();
            }
        }
    }

    @Override
    public boolean isValid() {
        if(_network.getRoot().getData().getHeight() > Utils.NET_MAX_HEIGHT) return false;

        if(Utils._Forbid_Net_Net) {
            for(Object nodeObj : _network.bfs()) {
                NetNode node = (NetNode) nodeObj;
                if(node.isNetworkNode()) {
                    for(Object parentObj : node.getParents()) {
                        NetNode parent = (NetNode) parentObj;
                        if(parent.isNetworkNode()) {
                            return false;
                        }
                    }
                }
            }
        }

        if(Utils._Forbid_Net_Triangle) {
            for(Object nodeObj : _network.bfs()) {
                NetNode node = (NetNode) nodeObj;
                if(node.isNetworkNode()) {
                    Iterator it = node.getParents().iterator();
                    NetNode parent1 = (NetNode) it.next();
                    NetNode parent2 = (NetNode) it.next();
                    if(parent1.hasParent(parent2) || parent2.hasParent(parent1)) {
                        return false;
                    }
                }
            }
        }
        /*Map<String, Double> constraints = TemporalConstraints.getTemporalConstraints(_geneTrees, _species2alleles, _alleles2species);
        Map<String, Double> lowerBound = TemporalConstraints.getLowerBounds(_network);
        for(String key : lowerBound.keySet()) {
            if(constraints.get(key) < lowerBound.get(key)) {
                return false;
            }
        }*/
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
        double[] likelihoodArray = new double[1];
        //likelihoodArray[0] = SNAPPLikelihood.computeSNAPPLikelihood(_network, _alleles2species, _markers, _BAGTRModel);

        if(!Utils.SAMPLE_SPLITTING) {

            if (SNAPPLikelihood.usePseudoLikelihood) {
                likelihoodArray[0] = SNAPPLikelihood.computeSNAPPPseudoLikelihood(_network, _alleles2species, _markers, _BAGTRModel); // pseudo likelihood
            } else if (SNAPPLikelihood.useApproximateBayesian) {
                likelihoodArray[0] = SNAPPLikelihood.computeApproximateBayesian(_network, _alleles2species, _markers, _BAGTRModel, _abcData); // approximate bayesian
            } else {
                likelihoodArray[0] = SNAPPLikelihood.computeSNAPPLikelihood(_network, _markers.get(0)._RPatterns, _BAGTRModel); // normal likelihood
            }
        } else {
            likelihoodArray[0] = SNAPPLikelihoodSampling.computeSNAPPLikelihoodST(_network, _splittings, _markers.get(0)._RPatterns, _BAGTRModel);
        }

        return likelihoodArray;
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
            tup.Item1.setData(new NetNodeInfo(height));
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

    private void initSplitting() {
        _splittings = new ArrayList<>();
        for(RPattern pattern : this._markers.get(0)._RPatterns.keySet()) {
            double count = this._markers.get(0)._RPatterns.get(pattern)[0];
            for (int k = 0; k < count; k++) {
                _splittings.add(new Splitting(this, _species2alleles, pattern, _BAGTRModel));
            }
        }

        if(_splittings.size() != _numSites) {
            throw new RuntimeException("Number of splittings != number of sites.");
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
}
