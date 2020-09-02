package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test;


import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HiddenState;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HmmBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HmmCore;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.RecombinationRate;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.TreeNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.Arrays;

/**
 * Build another Heliconious (butterfly) model for testing purposes.
 * This one has a higher lowest internal node.
 *
 * Created by Xinhao Liu on 06/02/20.
 */
public class NewButterflyModelBuilder {
    public static double trueRecombRate = 1E-8;

    public static ModelTree getButterflyModel() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("(((pachinus,cyd.gal.), melpomene), silvaniform);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(1000000);
        }
        for (TNode tnode:tree.postTraverse()) {
            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) tnode;
            if (!node.isLeaf() && !node.isRoot()) {
                if (node.getLeafCount() == 2) {
                    node.getData().setPopSize(2000000);
                    node.setNodeHeight(1520000);
                } else if (node.getLeafCount() == 3) {
                    node.getData().setPopSize(2000000);
                    node.setNodeHeight(5200000);
                }
            }
        }
        tree.getRoot().getData().setPopSize(2000000);
        tree.getRoot().setNodeHeight(8400000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getButterflyModelInit() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("(((pachinus,cyd.gal.), melpomene), silvaniform);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(1000000);
        }
        for (TNode tnode:tree.postTraverse()) {
            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) tnode;
            if (!node.isLeaf() && !node.isRoot()) {
                if (node.getLeafCount() == 2) {
                    node.getData().setPopSize(3000000);
                    node.setNodeHeight(3000000);
                } else if (node.getLeafCount() == 3) {
                    node.getData().setPopSize(3000000);
                    node.setNodeHeight(6000000);
                }
            }
        }
        tree.getRoot().getData().setPopSize(3000000);
        tree.getRoot().setNodeHeight(9000000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getButterflyModelMyValue() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("(((pachinus,cyd.gal.), melpomene), silvaniform);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(1000000);
        }
        for (TNode tnode:tree.postTraverse()) {
            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) tnode;
            if (!node.isLeaf() && !node.isRoot()) {
                if (node.getLeafCount() == 2) {
                    node.getData().setPopSize((int)2000000);
                    node.setNodeHeight(1215613.9161032941);
                } else if (node.getLeafCount() == 3) {
                    node.getData().setPopSize((int)2000000);
                    node.setNodeHeight(4822451.820160116);
                }
            }
        }
        tree.getRoot().getData().setPopSize((int)2000000);
        tree.getRoot().setNodeHeight(7761147.563291283);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getButterflyModelMyValue2() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("(((pachinus,cyd.gal.), melpomene), silvaniform);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(1000000);
        }
        for (TNode tnode:tree.postTraverse()) {
            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) tnode;
            if (!node.isLeaf() && !node.isRoot()) {
                if (node.getLeafCount() == 2) {
                    node.getData().setPopSize((int)2012551.378785898);
                    node.setNodeHeight(984663.9162390988);
                } else if (node.getLeafCount() == 3) {
                    node.getData().setPopSize((int)2185369.489401234);
                    node.setNodeHeight(4672820.261930535);
                }
            }
        }
        tree.getRoot().getData().setPopSize((int)2211527.458975877);
        tree.getRoot().setNodeHeight(7657788.510541027);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static void main(String[] args) {
        ModelTree model = getButterflyModel();
        HmmBuilder builder = new HmmBuilder(model.getTree(), model.getRecombRate());
        HmmCore hmm = builder.build();

        System.out.println(hmm.getStates().size());

        for (HiddenState state: hmm.getStates()) {
            System.out.println(state.getTree().toNewick());
        }
        System.out.println(Arrays.toString(hmm.getPi()));
        System.out.println(Arrays.deepToString(hmm.getA()));
    }
}
