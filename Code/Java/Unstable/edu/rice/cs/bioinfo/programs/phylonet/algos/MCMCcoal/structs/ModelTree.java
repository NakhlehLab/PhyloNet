package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.core.StateNode;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.IOException;
import java.util.Stack;


/**
 * Species tree to sample
 * This is one model parameter (others are population size hyper parameter and recombination rate)
 * Comparable to UltrametricNetwork in MCMC_SEQ
 * Created by Xinhao Liu on 11/1/19.
 */
public class ModelTree extends StateNode {
    private STITree<TreeNodeInfo> _tree;
    private RecombinationRate _recombRate;

    public ModelTree(String s, RecombinationRate rate, boolean heightIncluded) {
        this._recombRate = rate;
        try {
            this._tree = new STITree<>(s);
        } catch (IOException | ParseException e) {
            e.printStackTrace();
            System.exit(1);
        }
        buildSTITreeNodeHeight(_tree);
        // build pop sizes
        for (TNode node:_tree.postTraverse()) {
            STINode<TreeNodeInfo> stiNode = (STINode<TreeNodeInfo>) node;
            if (node.isLeaf()) {
                stiNode.setData(new TreeNodeInfo(Utils.N0));
            } else {
                stiNode.setData(new TreeNodeInfo(Utils.POP_SIZE_MEAN));
            }
        }
        // set label for each node for generating coalescent history
        int label = 1;
        for (TNode node:_tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight()); // might not be necessary here because parsed STITree already has parent distance set
            }
            if (!node.isLeaf()) {
                STINode<TreeNodeInfo> stiNode = (STINode<TreeNodeInfo>) node;
                stiNode.getData().setLabel(label++);
            }
        }
    }

    public ModelTree(String s, RecombinationRate rate) {
        try {
            this._tree = new STITree<>(s);
        } catch (IOException|ParseException e) {
            e.printStackTrace();
        }
        initTreeNodes();
        this._recombRate = rate;
        // continue?
    }

    /**
     * Initialize population size (in unit of individuals) and branch length (in unit of generations) for each branch
     * Initialize label of species tree (assume model tree topology does not change.)
     *
     * DO NOT USE THE NODE HEIGHTS INITIALIZED HERE!!! USE THE PREVIOUS CONSTRUCTOR TO PASS IN HEIGHTS
     */
    private void initTreeNodes() {
        Stack<Tuple<STINode<TreeNodeInfo>, Integer>> stack = new Stack<>();
        stack.add(new Tuple<>(_tree.getRoot(), Utils.DEFAULT_TREE_ROOT_HEIGHT));
        while (!stack.isEmpty()) {
            Tuple<STINode<TreeNodeInfo>, Integer> tup = stack.pop();
            int height = tup.Item2;
            tup.Item1.setNodeHeight(height);
            height *= Utils.TREE_INTI_SCALE;
            tup.Item1.setData(new TreeNodeInfo(Utils.POP_SIZE_MEAN));
            for (STINode<TreeNodeInfo> child:tup.Item1.getChildren()) {
                if (child.isLeaf()) {
                    child.setNodeHeight(Utils.DEFAULT_TREE_LEAF_HEIGHT);
                    child.setData(new TreeNodeInfo(Utils.POP_SIZE_MEAN));
                } else {
                    stack.add(new Tuple<>(child, height));
                }
            }
        }

        int label = 1;
        for (TNode node:_tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
            if (!node.isLeaf()) {
                STINode<TreeNodeInfo> stiNode = (STINode<TreeNodeInfo>) node;
                stiNode.getData().setLabel(label++);
            }
        }
    }

    /**
     * In setTreeBySample() in VariationalModel we would modify each node height and pop size.
     * Set each parent distance for the tree after node heights are modified so that the model tree is correct.
     */
    public void refresh() {
        // setparentdistance is actually not necessary. what else do we need?
        for (TNode node:_tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
    }

    /**
     * In setTreeBySample() in VariationalModel we would modify each pop size and some branch lengths.
     * Reset node heights and parent distances after some branch lengths are modified so that the model tree is correct.
     */
    public void refreshNewBranchLength() {
        // set node heights in _tree according to the branch lengths set in BranchLengthVariable()
        for (TNode n:_tree.postTraverse()) {
            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) n;
            if (node.isLeaf()) {
                node.setNodeHeight(0);
            } else {
                if (node.countLeafChildren() == 2) {
                    STINode<TreeNodeInfo> connectedLeaf = node.getChildren().iterator().next(); // TODO need to make sure this leaf is indeed the one passed to BranchLengthVariable
                    double height = connectedLeaf.getNodeHeight() + connectedLeaf.getParentDistance();
                    node.setNodeHeight(height);
                } else if (node.countLeafChildren() == 1) {
                    STINode<TreeNodeInfo> nonleafChild = null;
                    for (STINode<TreeNodeInfo> child : node.getChildren()) {
                        if (!child.isLeaf()) {
                            nonleafChild = child;
                            break;
                        }
                    }
                    double height = nonleafChild.getNodeHeight() + nonleafChild.getParentDistance();
                    node.setNodeHeight(height);
                } else if (node.countLeafChildren() == 0) {
                    STINode<TreeNodeInfo> leftChild = null;
                    for (STINode<TreeNodeInfo> child : node.getChildren()) {
                        leftChild = child;
                        break;
                    }
                    double height = leftChild.getNodeHeight() + leftChild.getParentDistance();
                    node.setNodeHeight(height);
                } else {
                    System.out.println("Shouldn't be here.");
                    System.exit(1);
                }
            }
        }
        // set all parent distances according to the newly set node heights
        for (TNode node:_tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
    }

    /**
     * STITree parsed from string does not have associated node height. It only has parent distance for each node.
     * This function builds the node height of each node of a newly parsed STITree.
     * Also used in HmmBuilder.
     */
    private void buildSTITreeNodeHeight(STITree<TreeNodeInfo> gt) {
        for (TNode node:gt.postTraverse()) {
            if (node.isLeaf()) {
                node.setNodeHeight(0);
            } else {
                double height = 0;
                for (TNode child:node.getChildren()) {
                    height = Math.max(height, child.getNodeHeight() + child.getParentDistance());
                }
                node.setNodeHeight(height);
            }
        }
    }

    @Override
    public double propose() {
        return 0;
    }

    @Override
    public void undo() {

    }

    @Override
    public void accept() {

    }

    @Override
    public void reject() {

    }

    @Override
    public double logDensity() {
        return 0;
    }


    @Override
    public boolean mayViolate() {
        return false;
    }

    @Override
    public boolean isValid() {
        return false;
    }

    public STITree<TreeNodeInfo> getTree() {
        return _tree;
    }

    public RecombinationRate getRecombRate() {
        return _recombRate;
    }

    //TEST
    public static void main(String[] args) {
        ModelTree tree = new ModelTree("((A,B),(C,D));", null);
        //ModelTree tree = new ModelTree("(((D,(F,E)),(C,(B,A))),(J,(K,(I,G))));", null);
        System.out.println(tree._tree.isRooted());
        System.out.println(tree._tree.toNewick());
        for (TNode plainNode:tree._tree.postTraverse()) {
            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) plainNode;
            System.out.println();
            System.out.println("Node name is " + node.getName());
            System.out.println("Node height is " + node.getNodeHeight());
            System.out.println("Node parent distance is " + node.getParentDistance());
            System.out.println("Node pop size is " + node.getData().getPopSize());
            System.out.println("Node label is " + node.getData().getLabel());
        }
//        System.out.println(tree.getTree().getLeafCount());
//        String[] leaves = tree.getTree().getLeaves();
//        for (String leaf:leaves) {
//            System.out.println(leaf);
//        }
        for (TNode leaf:tree.getTree().getRoot().getLeaves()) {
            System.out.println(leaf.getName());
        }
    }
}
