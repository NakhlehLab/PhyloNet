package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.tree;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.start.UPGMATree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.start.distance.JCDistance;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Moves the height of an internal node along the branch.
 * If it moves up, it can exceed the root and become a new root.
 * If it moves down, it may need to make a choice which branch to slide down into.
 * Created by wendingqiao on 2/15/16.
 * Adapted from
 * https://github.com/CompEvol/beast2/blob/master/src/beast/evolution/operators/SubtreeSlide.java
 */
public class SubtreeSlide extends TreeOperator {

    private TNode _target; // the target node
    private double _oldHeight; // the original height of target's parent
    private double _logHastingsRatio;
    private TNode _oldParent;
    private TNode _oldChild;
    private double _windowSize = 0.1;

    public SubtreeSlide(UltrametricTree tree) {
        super(tree);
    }

    @Override
    public double propose() {
        _violate = false;

        List<TNode> nodes = _tree.getInternalNodes();
        _target = nodes.get(Randomizer.getRandomInt(nodes.size()));
        // _target should be internal node not root node
        while(_target.isRoot()) {
            _target = nodes.get(Randomizer.getRandomInt(nodes.size()));
        }
        final double delta = getDelta();
        _oldHeight = _tree.getNodeHeight(_target.getParent());
        final double newHeight = _oldHeight + delta;

        TNode parent = _target.getParent();
        _oldParent = parent.getParent();
        _oldChild = super.getOtherChild(parent, _target);

        // invalid move : _parent's height cannot exceed _target's height
        if (_tree.getNodeHeight(_target) > newHeight) {
            _logHastingsRatio = Utils.INVALID_MOVE;
        }
        // preserve current topology, just change node height
        else if( ( newHeight > _oldHeight && (_oldParent == null || _tree.getNodeHeight(_oldParent) > newHeight) )
                || ( (newHeight < _oldHeight) && _tree.getNodeHeight(_oldChild) < newHeight ) )  {
            _tree.setNodeHeight(parent, newHeight);
            if(newHeight < _oldHeight) _violate = true;
            _logHastingsRatio = 0.0;
        }
        // topology will change
        // find new parent, new child for _parent to replace grandpar and sibling
        else {
            TNode newParent, newChild;
            if(newHeight > _oldHeight) {
                newParent = _oldParent;
                newChild = parent;
                while (_tree.getNodeHeight(newParent) < newHeight) {
                    newChild = newParent;
                    newParent = newParent.getParent();
                    if (newParent == null) break;
                }
                move(_target, newParent, newChild, newHeight);
                // count the hypothetical sources of this destination.
                _logHastingsRatio = -Math.log( intersectingEdges(newChild, _oldHeight, null) ); // possible sources

            } else {
                _violate = true;
                final List<TNode> newChildren = new ArrayList<>();
                _logHastingsRatio = Math.log( intersectingEdges(_oldChild, newHeight, newChildren) ); // possible dests
                if (newChildren.size() == 0) {
                    throw new IllegalArgumentException("the height of leaf should be smaller than " + newHeight);
                }
                final int idx = Randomizer.getRandomInt(newChildren.size());
                newChild = newChildren.get(idx);
                newParent = newChild.getParent();
                move(_target, newParent, newChild, newHeight);
            }
        }
        return _logHastingsRatio;
    }

    @Override
    public void undo() {
        if (_logHastingsRatio == Utils.INVALID_MOVE) return; // do nothing
        if (_logHastingsRatio == 0.0) {
            _tree.setNodeHeight(_target.getParent(), _oldHeight);
        } else {
            move(_target, _oldParent, _oldChild, _oldHeight);
        }
    }

    @Override
    public String getName() {
        return "SubtreeSlide";
    }

    private void move(TNode tar, TNode newParent, TNode newChild, final double newH) {
        TNode parent = tar.getParent();
        TNode grandpar = parent.getParent();
        TNode sibling = super.getOtherChild(parent, tar);

        // _parent's children {_target, sibling} => {_target, newChild}
        // _grandparent's children {_parent, ?} => {sibling, ?}
        // newParent's children {?, newChild} => {?, _parent}
        ((STINode) newChild).setParent2((STINode) parent);
        ((STINode) sibling).setParent2((STINode) grandpar);
        ((STINode) parent).setParent2((STINode) newParent);

        _tree.setNodeHeight(parent, newH); // reset br of {_target->_parent}, {newChild->_parent}
        if(grandpar != null) {
            _tree.setNodeHeight(grandpar, _tree.getNodeHeight(grandpar)); // reset br of {sibling->grandparent}
        }
        if(newParent != null) {
            _tree.setNodeHeight(newParent, _tree.getNodeHeight(newParent)); // reset br of {_parent->newParent}
        }
    }

    private int intersectingEdges(TNode node, double height, List<TNode> directChildren) {
        final TNode parent = node.getParent();

        if (_tree.getNodeHeight(parent) < height) return 0;

        if (_tree.getNodeHeight(node) < height) {
            if (directChildren != null) directChildren.add(node);
            return 1;
        }
        if (node.isLeaf()) {
            throw new IllegalArgumentException("the height of leaf should be smaller than " + height);
        }

        int count = 0;
        for(TNode child : node.getChildren()) {
            count += intersectingEdges(child, height, directChildren);
        }
        return count;
    }

    private double getDelta() {
        return (Randomizer.getRandomDouble() - 0.5) * _windowSize;
    }

    @Override
    public void optimize(final double logAlpha) {
        _windowSize *= Math.exp(Utils.calcDelta(this, logAlpha));
    }

    public static void main(String[] args) {

        Map<String, String> locus = new HashMap<>();
        {
            locus.put("A", "CTTCGTGACGGGCTCGGCTCGTACGGTCAAGGGCACCTGAGCTAGGCAACTCAAGACGGGCGAGAGTCCCGCTACTACGGAAAAGGGTTTGCAGCTTCGAACATAGAACCACGGGACCTCCGGGACTCGCCCAATTGCGGTCGATGGCACTACGATTGGACAGGCGCTTACGCCAATTTACAGGTAGCGAGAGGTCTGCGCAGCAAAAGACTCCGCTTCACCCGGTGGTACCTGACTCGCGGCGCTGACCGTTCGCCGTAAGCGTTGACAATTTTCGGAAGCATCGCCAAAGTGCTGGAGTACGGGCTTCACACCCGCAACCACCGACAGCCCAACGAATCACTTCGCCGGTTGTCCTTCACCCCTGTTGGCAGAGGATGCTTGCGTGACTTTCATCCTTCTGTCCTTTCGTGCGGACTCGCACAGATCCTCCAAACAAGCGAGATCCGACCGATACTCTGCCCCCAGCAAGGCGGGGTTTCAGCGTCCCGTAGCTAGTAGGTATCCAGTCCAGAGGCACGTAGAGCACTCACCGTCCCCAGCCCCCCTCCACTCATCGCTGTGTTGAGGTCAAGTTCCCCGAGATCGTACTAACCGGTATGCAACCCACATCTGCCCATAGTCTACCCATCATACTACAGTCAGTAAACCCGGAGTGTATGGCTCGACTAATGTCCGTACAGCCAGCGCTCGATAAGTGCGTGCTCGAACGCTTACACCCCTGCTGGCACCAGGTAGCGCATGATCCACCAATGCAGCCCGAGCGCTTGCTGTGACTTCCGTCCGAAAGTATCGAGTTAAGCCTCAGCAATATGGGCCACCTTGACTGAGCTACACTCCCTCCCGTACGGGTGAAGTCCCCGCCGGCACAGGGGGGACACTTCTATTGAACATTTCGCCTCACCGAAGTCGCCCCCGAGTACGCCAGTGACGTACAGTCGCACGCGGGTTTGGTAAGGAGGAGCCAGATCGCGTCTACATGTTGAGTAGGTCCCGGCGC");
            locus.put("C", "CTACGTGACGGGCCCGGCTCGCACGGTCAATGGCACCTGCGCTAGGCAACTCGAGACGGGCGAGGGTCCCGCTACTACGGACAAGGGTGTGCAGCGTCGAACATAGAACCACGGGACCTCCGGTACTCGCTCAACCGCGGTCGGTGGTACTACGATTGGACGGGCACTTTCGGCAATTTACAGTTAGCGAGAGGTCTGCGCGGCAAAAGACTCCAATTCACCCGGCGTTACCTGACTCGCGGCACTGACCGTTCGCCGTAAGCGTTGACAATTTTCGGTAGCATCGCCAAAGTGCTGGAGTACGGGTTTCACACTTGTAGCCACCGACAGCCCAACGGATCAATTCGCCAGTTGCCCTTCACTCCTGTTAGCAGAGGATGCTAGCGTGAGTTTCATCCATCTGTCCTGTCGAGCGGACTCGCACAGATACTCCAAACAAGTGAGATCCGTCCGGTACTCTACCCCTCACAAAGCGGGATTTCAGCGTCTCGTAGTGAGTAGGCATCCGGTCCAGGGGCACGTAGAGCACTCACCGTCCTCAGCCCCCTTCTATTCATTGCTGTGTTGAAGCGAAGTTCCCCGAAATCGTACTAACCGGTATGCAACCCAGATCTGCCCATAGTTTCCACATCATACTACAGCCAGAAAACTCGGAGTGTATGGCTTGACTAATGTCCGTACTGCCAGCGCTAGGTAAGTACGTGCTCGAACGCTTACACCCCTGCTGGCACCTGGTTGCGCATGATCCACCAATGCAGCCCGAGCGCTTGCTGTGACTTCCGTCCAAAAGTATTGAATTAAGCCTCAGCAATATGGGCCACCCAGACTGAGCCACGCTCCCCTCCGTACGGATGAAGTCCCCACCGGGGCAGGGGGGACGCTTCTATTGAACACTTCGCCTCACCGAAGTCGCCCCCGAGTACGCCAGTGACGTACGGTCGCACGCGGGTTCGGTAAGGAGGGGTCAGATCGCGTCTGCATGTTGAGTAGGTCCCTGCGC");
            locus.put("G", "CTACGTGACGGGCCTGGCTCGCACGGTCAATGGCACTTGAGCTAGGCAACCCAAGACGGGCGAGGGTCCCGCTACTACGGAAAAGGGTGTGCAGCTTCGAACATAGAACCACGGGACCTCCGGGACTCGCTCAATTGCGGTCGGTGGTACTACGATTGGACGGGCACTTTCGGCAATTTACAATTAGCGAGAGGTCTGCGCAGCAAAAGACTCCAATTCACCCGGCGCTACCTGATTCGCGGCACTGACCGTTCGCCGTAAGCGTTGACAAATTTCGGTAGCATCGCCAAAGTGCTGGAGTACGGGCTTCACACTTGTAGCCACCGACAGCCCAACGAATCAATTCGCCAGTTGCCCTTCACTCCTGTTAGCAGAGGATGCCAGCATGAGATTCATCCATCTGTCCTTTCGAGCGGACTCGCACAGATACTCCAAACAAGCGAGATCCGACCGGTACTCTACCCCTCACAAAGCGGGATTTCAGCGTCTCGTAGTGAGTAGGTATCCGGTCCAGAGGCACGTAGAGCCCTCACTGTCCTCAGCCCCCCTCTACTCATTGCTGTGTTGAGGTGAAGTTCCCCGAAATCGTACTAACCGGTATGCAACCCAGATCTGCCCATAGGCTACACATCATACTACAGCCAGTAAACCCGGAGTGTATGGCTTGACTAAGGTCCGTACAGCCAGCGCTCGATAAGTACGTGCTCGAACGCTTACACCCCTGCTGGCACCAGGTAGCGTATGATCCACCGATGCAGCCCGAGCGCTTGCTGTGACTTCCGTCCGAAAGTATTGAGTTAAGCCTCAGCAATATGGGCCACCCTGACTGAGCTACGCTCCCTTCCGTACGGATGAAGTCCCCACCGGGACAGGGGGGACGCTTCTATTGAACACTTCACCTCACCGAAGTCGCCCCCGAGTACACCAGTGACGTACAGTCGCACGCGGGTTCGATAAGGAGGGGCCAGATCGCGTCTACATGTTGAGTAGGTCCCTGCGC");
            locus.put("R", "TTGTGCGACGGGCGCGGCCCGGGCGATCAAGGGTACCAGAGTCAGGCAACTTAAGATGGGCGAGAGCCCCGCTAATACGGAAAGGAGTATGCAGCTCCGGACATGGGATCACTGGATCTCCAGGACTCGGTCGACTGCGGTCGGCGTTACTACAATCGGAGAGGCACATTCGCTAACTCATAGTTTGCGAGAGATCTGCGGAGCAAAGGATTCCCCTCTACTCGGCGCTACCTGACTCGCTACGCAGACCGTTCGCCGTAAGTGTTGCCAATCCCCGGAGGCATCGCCAAAGTACTGGAGTACGGGCTTAACACCTACAGCCACCGACAGCGCAACGAATCATTTCACCAGTTGCCTTTTACCCTTGTTAGCTGAGGATGCTAGCTTAACTTTCATCCTATGTTCCTCTCGGGCGGACTCTAACTGATCCCCCGAATGAGCGAGGTCTGACCGATACTCTGCCCCAGGCAAGGGGGGAGCCCATCCCCTTATAGTAAGCAAATCCCCAGTTCAGAGGCACATAGAGCACCCACCGTCCACAGCCCCTTTCCACTCACTGGTGCGCTGAGGTGAAAGTGCCCGAAATCCTACGAATTGGTATGCAACCCAGATCTGTAGGCAGGCTACATGTCATACTACAGCTCGTAAACTTGGAGTGTATGGCTAGACTGATATCCGAACAAACAACGCTCGACAAGCGCGACCTCGACCGCTCACACCCTTGCTGACACCAAGCAACACATGATCCATCAGTGCAGCCCCAACGTTTTTTGTGACCTCCGTCCGAAAGTATGGATTTGAGCCTCAGCAATGTGGCCCACCATGGCCGAGCTACGCTCCCCTACGTACGGATGATTTCCCCGCCGGGACAGGCGGGACGGTTCTATTAAACATCTCACCTTACTGATGTCGCCCCCGGGTACGGCAGCGACGTACAGCCGCACGCGAGCTTGGTAAGGAGGAGCCAGATCGTGCCTACATGTTGAGTAGGTCCCTACTC");
            locus.put("Q", "CCATACGATGGGCTCGGCTCGTATGATTGAGGGCACCGAAGCTAGGCGACTCAAGATGGGCGAGGGCTCCGCGAATACGGAAAAGGGTATGCAGCTTCGGCCATAGGACCACGTGATCTCCGGGACTCGCCCAATTGAGATCGGCGTTACTACAATTAGACGAGCACATTCGCCAATTTATAGTTAACAAGAGATCTGCGTAGCAAAAGATTCAACGTTACCCGGCGCTACCAGACCCGCGGCACAGGCCGTTTGCCCTAAGCGTTGACAATCTTCGGATTCATCGCCAAAGTGCTGGAGTACAGGCTTCACACCTACAGCCACTGACAGCCTAACGAATCACTTCACCAATTGCCTTTCACCCCTGTTAGCGGAGGGCGCTAGCATAACTTTCGTCCTACGTTCCTCTCGTACGGATTCGGACAGATCCTCCGAGCAAGCGAAGTCCGGCCGATACTCTGCCCCTAGCAAGGCGGGATCCCGTCGCCTTGTAGTGAGCAGATATCCAGTTTGGGAGCACATAGAGCACCGACCGTCCACAGTCCCCTTCTCTTCATTGGTGCGTTGAGGTGAAAGTTCCCGAAATCCTACAAACTGGTATGCAAACCGGATCTGCAGGTTGGCTACGTATCATACTACAGCCCGTAAACTCGGAGTGCATGGCTTGACTAACATCCGTACAAACAGCGCTCGATGAGTGCGACCTCGGCAGCTTACACCCTTGCTGACACCAAATAGCGCATGATCCACCAGTACAGCCCAAGCGCCTTTCGCGTCCTCCGCCCGAATGTATGGATGTAAGCCTGAGTCACGTGGACCACCGTGCCCGAGCTACGCTCCCTTACGTGCGGATGATGTCCCCGCCGGGACAGGCGGAACGCTTCTATTGAACATTTCACCTCGCTGAAGTCGCCCCTGAGTACGCCAGCAACGTACAGTCGCACGCGAGCTCAGTAAGGAGGAGCCATATCGCGTCTACATGATGCGTAGGTCCCTGCGC");
            locus.put("L", "CTACGCGACGGGCTCGGCTCGTACAATAAAGGGCACCAGAGCTAGGTAACCCAAGATAGGCGAGGGCTTCGCTAATACGGAGAAGGGTATTCAGCTTCGGTCACAGGGCCACGTGATCTCCGGGACTCGCCCAGTTGCGATCGGCGTTACTACAATTGGAGGGGCACATTCGCCAACTTACAGTTAGTGAGGGATCTGCGTATCAAAAAACTCCACACTACTCGGCGCGACCTGGCTCGTGGCGCAGGCCGTTTGCCGTAAGTGTTGAGAATTTTCGGAAGAATGGCCGGAGTGTTAGAGTCCAGGCTTCACACCTACAGCCATTGACAGCTCAACGAGTCATTTCACCAGTTGCCTTTCACCCCTGTTAGCAGAGGGAGCTAGCGCAGCTTTCATCCTATGTCCATCTCGGGCGGACTCGGACAGATCCTCCGAACAAGCGAGGTCCGGTCGATACTCTGCCCCTAACAAGGCAGGATCCCATCGCCTTGTAGTGAGCAGATATCCAGTTTAGGGGCGCATAGAGCACCCACCGTCCACAGCCCCCTTGTACTTATTGGTGTGTTGAGGTGAAAGTCCCCGAAATCCCAGGAACTTGTATGCAACCCAGATCTGCAGGTAGGCTACGTATCATACTGTACCCCGCAAACTCGGAGTGTGTGACTGGACTAACGTCCGTACAAACAGCGCTCGATAGGTGCGACCTCGACAGCTTACACCCTTGCTGGCACCAAATAGCGCGTGATCCACCAGTGCAGACCAAACGTTTTCTGCGCCCTCCGTCCGAAAGTACGGGTTTAAGTCTCAGTAACATGGTCCACTATGACCGAACTATACCCCCTTACGTACGGGTGATGTCCCCGCCGGGACGGGCGGGACGCTTCTACTGAACGCTTCACTTCACTGAAGTCGCCCCTGAGTATGCCAGCAACGTACAGTCGCACGCGAGCTCGGTAAAGAGGAACCAGATCCCGTCGACATGTTGAGTTGGTCCCTATGC");
        }
        Alignment seq = new Alignment(locus);
        UPGMATree upgma = new UPGMATree(new JCDistance( seq.getAlignment() ));
        UltrametricTree template = upgma.getUltrametricTree();
        System.out.println(template.getTree().toNewick());
        {
            UltrametricTree tree = new UltrametricTree(template);
            System.out.println(tree.getTree().toNewick());

            int runs = 10000;
            int counter = 0;
            int test = 0;
            for(int i = 0; i < runs; i++) {
                Operator op = new SubtreeSlide(tree);
                double logHR = op.propose();
                if(logHR != Utils.INVALID_MOVE && logHR != 0.0) counter++;
                if(tree.checkUltrametric()) test++;
            }
            System.out.println(test == runs);
            System.out.println(tree.getTree().toNewick());
            System.out.printf("%d out of %d\n", counter, runs);
        }
        {
            UltrametricTree tree = new UltrametricTree(template);
            System.out.println(tree.getTree().toNewick());

            int runs = 10000;
            int counter = 0;
            int test = 0;
            for(int i = 0; i < runs; i++) {
                Operator op = new SubtreeSlide(tree);
                double logHR = op.propose();
                if(logHR != Utils.INVALID_MOVE && logHR != 0.0) counter++;
                op.undo();
                if(tree.compareTo(template) == 0) {
                    test++;
                } else {
                    System.out.println(tree.getTree().toNewick());
                }
            }
            System.out.println(test == runs);
            System.out.printf("%d out of %d\n", counter, runs);
            System.out.println(tree.getTree().toNewick());
        }
    }
}
