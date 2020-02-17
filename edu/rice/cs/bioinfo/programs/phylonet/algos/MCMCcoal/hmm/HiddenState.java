package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Hidden state of coal-HMM.
 * Created by Xinhao Liu on 12/30/19.
 */
public class HiddenState {
    // Each coalescent history has its own state. Branch lengths of coalescent trees are averaged across each
    // coalescent history. Branch lengths are in unit of 4N0 generations.
    // ATTENTION: for fsc, branch lengths are in unit of generations (fractions of generations)
    private Tree coalescentTree;
    private int index;

    public HiddenState(Tree tree, int index) {
        this.coalescentTree = tree;
        this.index = index;
    }

    public Tree getTree() {
        return coalescentTree;
    }

    public int getIndex() {
        return index;
    }

    //TODO: beagle likelihood?

    //TODO: List of Alignment in util?
    // what is clock rate?

    private Map<TNode, Integer> nodeLabel;

    public TNode[] getNodeArray() {
        TNode[] array = new TNode[coalescentTree.getNodeCount()];
        nodeLabel = new HashMap<>();

        String[] taxon = coalescentTree.getLeaves();
        Arrays.sort(taxon);
        for(int i = 0; i< taxon.length; i++) {
            array[i] = coalescentTree.getNode(taxon[i]);
            nodeLabel.put(array[i], i);
        }
        int internalIdx = taxon.length;
        for(TNode node : coalescentTree.postTraverse()) {
            if(node.isLeaf()) continue;
            nodeLabel.put(node, internalIdx);
            array[internalIdx++] = node;
        }
        return array;
    }

    public int getNodeLabel(TNode node) {
        return nodeLabel.get(node);
    }
}
