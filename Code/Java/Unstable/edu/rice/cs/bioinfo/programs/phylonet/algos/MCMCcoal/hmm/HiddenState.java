package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.IOException;
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

    /**
     * ***************ONLY FOR TESTING HCG!***************
     * THIS IS HARD CODED!
     * Convert gene tree to its name in one of {HC1, HC2, HG, CG}
     */
    public String getStateName() {
        try {
            if (Trees.haveSameRootedTopology(coalescentTree, new STITree<>("((H,G), C);"))) {
                return "HG";
            } else if (Trees.haveSameRootedTopology(coalescentTree, new STITree<>("((C,G), H);"))) {
                return "CG";
            } else {
                if (Trees.getInternalNodes(coalescentTree).get(0).getNodeHeight() >= 5.5) {
                    return "HC2";
                } else if (Trees.getInternalNodes(coalescentTree).get(0).getNodeHeight() < 5.5) {
                    return "HC1";
                }
                System.out.println("ERROR!!!!!!!");
                System.exit(1);
            }
        } catch (IOException | ParseException e) {
            e.printStackTrace();
        }
        return null;
    }

    /**
     * ***************ONLY FOR TESTING cydno-timareta-numata!***************
     * THIS IS HARD CODED!
     * Convert gene tree to its name in one of {CT1, CT2, CN, TN}
     *
     * rootHeight is in coalescent units!
     */
    public String getStateNameButterfly(double rootHeight) {
        try {
            // ((cydno, timareta), numata);
            if (Trees.haveSameRootedTopology(coalescentTree, new STITree<>("((cydno, numata), timareta);"))) {
                return "CN";
            } else if (Trees.haveSameRootedTopology(coalescentTree, new STITree<>("((timareta, numata), cydno);"))) {
                return "TN";
            } else {
                assert Trees.haveSameRootedTopology(coalescentTree, new STITree<>("((cydno, timareta), numata);"));
                if (Trees.getInternalNodes(coalescentTree).get(0).getNodeHeight() >= rootHeight) {
                    return "CT2";
                } else if (Trees.getInternalNodes(coalescentTree).get(0).getNodeHeight() < rootHeight) {
                    return "CT1";
                }
                System.out.println("ERROR!!!!!!!");
                System.exit(1);
            }
        } catch (IOException | ParseException e) {
            e.printStackTrace();
        }
        return null;
    }
}
