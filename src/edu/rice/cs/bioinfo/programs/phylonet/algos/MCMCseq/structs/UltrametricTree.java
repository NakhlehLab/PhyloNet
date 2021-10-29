package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs;

import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core.StateNode;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.likelihood.BeagleTreeLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.sitemodel.SiteModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.substitution.Frequencies;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.substitution.GTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.substitution.JukesCantor;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.substitution.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.tree.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.start.UPGMATree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.start.distance.JCDistance;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.StringReader;
import java.util.*;

/**
 * Ultrametric tree
 * Created by wendingqiao on 1/27/16.
 */
public class UltrametricTree extends StateNode implements Comparable<UltrametricTree> {

    // ----- sequences -----
    private SubstitutionModel _substModel;
    private Alignment _alignment;

    // ----- tree states -----
    private Tree _tree;
    private Tree _oldTree;
    private List<TNode> _nodes;
    private List<TNode> _internalNodes; // contains tree root
    private double _clockRate = 1.0; // TODO can be inferred from gamma distribution
    private double _mutationRate = 1.0; //Todo by zhen: for mutation rate

    // ----- likelihood -----
    private BeagleTreeLikelihood _beagle;
    private double _logL = Utils.INVALID_MOVE;
    private double _logLtemp = Utils.INVALID_MOVE;

    // ---- moves -----
    private double[] _opWeights;
    private TreeOperator[] _operators;


    public UltrametricTree(Alignment aln) {
        this._alignment = aln;
        if(Utils._SUBSTITUTION_MODEL.equals("GTR")) {
            Frequencies freq = Utils._BASE_FREQS == null ?
                    new Frequencies(aln, Utils.ESTIMATE_SUBSTITUTION) :
                    new Frequencies(aln, Utils._BASE_FREQS);
            this._substModel = Utils._TRANS_RATES == null ?
                    new GTR(freq) : new GTR(freq, Utils._TRANS_RATES);
        } else {
            Frequencies freq = new Frequencies(aln, false);
            this._substModel = new JukesCantor(freq);
        }

        UPGMATree temp = new UPGMATree(new JCDistance( this._alignment.getAlignment() ));
        this._tree = temp.getTree();
        this._nodes = IterableHelp.toList(_tree.getNodes());
        for(TNode node : _tree.getNodes()) {
            if(node.getParent() == null) continue;
            if(node.isLeaf()) {
                node.setParentDistance(getNodeHeight(node.getParent()));
            } else {
                node.setParentDistance(getNodeHeight(node.getParent()) - getNodeHeight(node));
            }
        }
        this._internalNodes = Trees.getInternalNodes(_tree);

        this._beagle = new BeagleTreeLikelihood(aln, this, new SiteModel(this._substModel), null);

        setOperators();

        GTBurnin();
    }


    public UltrametricTree(String newick, Alignment aln) {
        try {
            this._alignment = aln;
            if(Utils._SUBSTITUTION_MODEL.equals("GTR")) {
                Frequencies freq = Utils._BASE_FREQS == null ?
                        new Frequencies(aln, Utils.ESTIMATE_SUBSTITUTION) :
                        new Frequencies(aln, Utils._BASE_FREQS);
                this._substModel = Utils._TRANS_RATES == null ?
                        new GTR(freq) : new GTR(freq, Utils._TRANS_RATES);
            } else {
                Frequencies freq = new Frequencies(aln, false);
                this._substModel = new JukesCantor(freq);
            }

            NewickReader nr = new NewickReader(new StringReader(newick));
            STITree tree = new STITree<>();
            nr.readTree(tree);
            this._tree = tree;
            if(_tree.getRoot().getChildren().iterator().next().getParentDistance() != TNode.NO_DISTANCE) {
                buildNodeHeightMap();
            } else {
                UPGMATree temp = new UPGMATree(new JCDistance( this._alignment.getAlignment() ));

                for(TNode node : temp.getTree().postTraverse()) {
                    if(node.isLeaf()) {
                        node.setNodeHeight(Utils.DEFAULT_TREE_LEAF_HEIGHT);
                    } else {
                        double height = 0;
                        for(TNode child : node.getChildren()) {
                            height = Math.max(height, getNodeHeight(child) + child.getParentDistance());
                        }
                        node.setNodeHeight(height);
                    }
                }

                double height = temp.getTree().getRoot().getNodeHeight();
                double step = (height - 1e-4) / (this._tree.getLeafCount() -1 );

                double hi = 1e-4;
                for(Object nodeObj : this._tree.postTraverse()) {
                    STINode node = (STINode) nodeObj;
                    if(node.isLeaf()) node.setNodeHeight(Utils.DEFAULT_TREE_LEAF_HEIGHT);
                    else {
                        node.setNodeHeight(hi);
                        hi+=step;
                    }
                }

                for(Object nodeObj : this._tree.postTraverse()) {
                    STINode node = (STINode) nodeObj;
                    if(!node.isRoot()) {
                        node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
                    }
                }
            }
            this._nodes = IterableHelp.toList(this._tree.getNodes());
            this._internalNodes = Trees.getInternalNodes(_tree);

            this._beagle = new BeagleTreeLikelihood(aln, this, new SiteModel(this._substModel), null);

            setOperators();

        } catch (Exception ex) {
            ex.printStackTrace();
        }

        GTBurnin();


    }

    // This constructor is only used for debug only
    public UltrametricTree(Tree tree) {
        this._tree = tree;
        buildNodeHeightMap();
        _nodes = IterableHelp.toList(_tree.getNodes());
        this._internalNodes = Trees.getInternalNodes(_tree);
    }
    public UltrametricTree(UltrametricTree temp) {
        this._tree = new STITree<>(temp._tree.getRoot().getName(), true);
        copyNode(this._tree.getRoot(), temp._tree.getRoot());
        _nodes = IterableHelp.toList(_tree.getNodes());
        this._internalNodes = Trees.getInternalNodes(_tree);
    }

    private void GTBurnin() {
        if(Utils._START_GT_BURN_IN) {
            int iter = 0;
            while (iter < 10000) {
                double logHR = this.propose();
                this.setDirty(true);
                if(logHR != Utils.INVALID_MOVE) {
                    double newLogL = this.logDensity();
                    if(newLogL > _logL) {
                        this.accept();
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

            System.out.println("Optimized starting gene tree for " + _alignment.getName());
        }
    }

    public void resetTree(Tree tree) {
        this._tree = new STITree<>(tree.getRoot().getName(), true);
        copyNode(this._tree.getRoot(), tree.getRoot());
        getNodeArray();
        buildNodeHeightMap();
        _dirty = true;
        _nodes = IterableHelp.toList(_tree.getNodes());
        this._internalNodes = Trees.getInternalNodes(_tree);
        this._beagle.reset();
        //this._beagle = new BeagleTreeLikelihood(_markerSeq, this, new SiteModel(this._substModel), null);

    }

    private void setOperators() {
        this._operators = new TreeOperator[] {
                new TNodeReheight(this),
                new SubtreeSlide(this), new NarrowNNI(this),
                new TreeScaler(this), new TreeRootScaler(this),
                new WilsonBalding(this), new WildNNI(this)
        };

        double[] Tree_Op_Weights = Utils.Tree_Op_Weights.clone();

        if(Utils._FIX_GENE_TREE_TOPOLOGIES) {
            Tree_Op_Weights[1] = 0.0;
            Tree_Op_Weights[2] = 0.0;
            Tree_Op_Weights[5] = 0.0;
            Tree_Op_Weights[6] = 0.0;
        }

        this._opWeights = Utils.getOperationWeights(Tree_Op_Weights, true);
    }

    @Override
    public int compareTo(UltrametricTree o) {
        Map<TNode, TNode> nodeMap = Trees.mapTwoTopologies(this._tree, o._tree);
        if(nodeMap == null) return -1;
        for(TNode key : nodeMap.keySet()) {
            if(key.isLeaf()) continue;
            if(getNodeHeight(key) == TNode.NO_NODE_HEIGHT ||
                    getNodeHeight(nodeMap.get(key)) == TNode.NO_NODE_HEIGHT) {
                return -1;
            }
            if(Math.abs(getNodeHeight(key) - getNodeHeight(nodeMap.get(key))) > 0.000001) {
                return -1;
            }
        }
        return 0;
    }

    public List<TNode> getNodes () {
        return _nodes;
    }

    public List<TNode> getInternalNodes () {
        return _internalNodes;
    }

    public Tree getTree() {
        return _tree;
    }

    public double getNodeHeight(TNode node) {
        return node.getNodeHeight();
    }

    public void setNodeHeight(TNode node, double height) {
        if(node.isLeaf()) {
            System.err.println("The height of a leaf node cannot be modified!");
            return;
        }
        double[] bound = getNodeHeightBounds(node); // lower bound, upper bound = parent's height
        if(height < bound[0] || height > bound[1]) {
            System.err.printf("The new height %.4f is out of bound %s", height, Arrays.toString(bound));
            return;
        }
        if(!node.isRoot()) {
            node.setParentDistance(bound[1] - height);
        }
        for(TNode child : node.getChildren()) {
            child.setParentDistance(height - getNodeHeight(child));
        }
        node.setNodeHeight(height);
    }

    /**
     * get the lower and upper bound of node height
     * @param node
     * @return  [lower, upper] bound
     */
    public double[] getNodeHeightBounds(TNode node) {
        if(node.isLeaf()) {
            System.err.println("The height of a leaf node cannot be modified!");
            return null;
        }
        double[] bound = new double[2];
        bound[1] = node.isRoot() ? Double.NaN : getNodeHeight(node.getParent());
        for(TNode child : node.getChildren()) {
            bound[0] = Math.max(getNodeHeight(child), bound[0]);
        }
        return bound;
    }

    private double[] _prevHeights = null;

    public void scale(double scaleFactor, boolean undo) {
        if(undo) {
            for(int i = 0; i < _internalNodes.size(); i++) {
                _internalNodes.get(i).setNodeHeight(_prevHeights[i]);
            }
            _prevHeights = null;
        } else {
            _prevHeights = new double[_internalNodes.size()];
            for(int i = 0; i < _internalNodes.size(); i++) {
                TNode node = _internalNodes.get(i);
                double nh = getNodeHeight(node);
                _prevHeights[i] = nh;
                node.setNodeHeight(nh * scaleFactor);
            }
        }
        for(TNode node : _tree.getNodes()) {
            if(node.isRoot()) continue;
            node.setParentDistance(getNodeHeight(node.getParent()) - getNodeHeight(node));
        }
    }

    public boolean checkUltrametric() {
        double diff;
        for(TNode node : _tree.getNodes()) {
            if(node.isRoot()) continue;
            diff = Math.abs(getNodeHeight(node) + node.getParentDistance() - getNodeHeight(node.getParent()));
            if(diff > 0.000001) {
                System.err.println(_tree.toString());
                System.err.println(diff + " " + node.getName() + "<-" + node.getParent().getName());
                throw new IllegalArgumentException();
            }
        }
        return true;
    }

    private void buildNodeHeightMap() {
        for(TNode node : _tree.postTraverse()) {
            if(node.isLeaf()) {
                node.setNodeHeight(Utils.DEFAULT_TREE_LEAF_HEIGHT);
            } else {
                double height = 0;
                for(TNode child : node.getChildren()) {
                    height = Math.max(height, getNodeHeight(child) + child.getParentDistance());
                }
                node.setNodeHeight(height);
            }
        }
    }

    public void resetNodeHeights() {
        for(Object o : _tree.postTraverse()) {
            TNode node = (TNode) o;
            if(node.isLeaf()) {
                node.setNodeHeight(Utils.DEFAULT_TREE_LEAF_HEIGHT);
                continue;
            }
            double height = 0;
            double update;
            for(TNode child : node.getChildren()) {
                update = child.getParentDistance() + getNodeHeight(child);
                height = Math.max(height, update);
            }
            node.setNodeHeight(height);
            for(TNode child : node.getChildren()) {
                child.setParentDistance(height - getNodeHeight(child));
            }
        }
    }

    private void copyNode(TNode node, TNode temp) {
        node.setNodeHeight(getNodeHeight(temp));
        for(TNode tempChild : temp.getChildren()) {
            STINode child =((STINode) node).createChild(tempChild.getName());
            child.setParentDistance(tempChild.getParentDistance());
            copyNode(child, tempChild);
        }
    }

    @Override
    public double logDensity() {
        if(_logL == Utils.INVALID_MOVE) {
            _logL = _beagle.calculateLogP();
            return _logL;
        }
        if(!_dirty) return _logL;
        _logLtemp = _beagle.calculateLogP();
        return _logLtemp;
    }

    @Override
    public boolean mayViolate() {
        return _operator.mayViolate();
    }

    @Override
    public void accept() {
        _dirty = false;
        if(_logLtemp != Utils.INVALID_MOVE) _logL = _logLtemp;
        _logLtemp = Utils.INVALID_MOVE;
    }

    @Override
    public void reject() {
        _dirty = false;
        _logLtemp = Utils.INVALID_MOVE;
    }

    @Override
    public boolean isValid() {
        return checkUltrametric();
    }

    @Override
    public double propose() {
        _oldTree = Trees.readTree(_tree.toNewick());
        this._operator = getOp(_operators, _opWeights);
        return this._operator.propose();
    }

    @Override
    public void undo() {
        if(this._operator == null) throw new IllegalArgumentException("null operator");
        this._operator.undo();
        //if(!_oldTree.toNewick().equals(_tree.toNewick())) {
        //    throw new RuntimeException("!!!!!! " + _oldTree.toNewick() + "\n" +_tree.toNewick() );
        //}
    }

    public String toString() {
        return this.getTree().toNewick();
    }

    private Map<TNode, Integer> nodeLabel;

    public TNode[] getNodeArray() {
        TNode[] array = new TNode[_tree.getNodeCount()];
        nodeLabel = new HashMap<>();

        String[] taxon = _tree.getLeaves();
        Arrays.sort(taxon);
        for(int i = 0; i< taxon.length; i++) {
            array[i] = _tree.getNode(taxon[i]);
            nodeLabel.put(array[i], i);
        }
        int internalIdx = taxon.length;
        for(TNode node : _tree.postTraverse()) {
            if(node.isLeaf()) continue;
            nodeLabel.put(node, internalIdx);
            array[internalIdx++] = node;
        }
        return array;
    }

    public int getNodeLabel(TNode node) {
        return nodeLabel.get(node);
    }

    public void updateMutationRate(double mu){
        this._mutationRate = mu;
        this._beagle.updateMutationRate(_mutationRate);
        this._beagle.reset();
    }

    public double get_mutationRate(){
        return _mutationRate;
    }
}
