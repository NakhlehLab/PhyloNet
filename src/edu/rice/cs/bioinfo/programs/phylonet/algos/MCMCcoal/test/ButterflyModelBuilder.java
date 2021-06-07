package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test;


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
 * Build Heliconious (butterfly) models for testing purposes.
 * Created by Xinhao Liu on 05/29/20.
 */
public class ButterflyModelBuilder {
    public static double trueRecombRate = 1E-8;

    public static ModelTree getButterflyModel() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("(((mel.agl.,mel.ama.), tim.ssp_nor.), silvaniform);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(2000000);
        }
        for (TNode tnode:tree.postTraverse()) {
            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) tnode;
            if (!node.isLeaf() && !node.isRoot()) {
                if (node.getLeafCount() == 2) {
                    node.getData().setPopSize(2000000);
                    node.setNodeHeight(400000);
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
        ModelTree model = new ModelTree("(((mel.agl.,mel.ama.), tim.ssp_nor.), silvaniform);", dummyRecombRate);
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
        ModelTree model = new ModelTree("(((mel.agl.,mel.ama.), tim.ssp_nor.), silvaniform);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(1000000);
        }
        for (TNode tnode:tree.postTraverse()) {
            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) tnode;
            if (!node.isLeaf() && !node.isRoot()) {
                if (node.getLeafCount() == 2) {
                    node.getData().setPopSize((int)1106757.2978256762);
                    node.setNodeHeight(974581.7822515282);
                } else if (node.getLeafCount() == 3) {
                    node.getData().setPopSize((int)1820047.324083828);
                    node.setNodeHeight(5660955.953852965);
                }
            }
        }
        tree.getRoot().getData().setPopSize((int)1673528.0050732098);
        tree.getRoot().setNodeHeight(8640261.550110878);
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
    }
}
