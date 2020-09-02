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
 * Build ABC models for testing bin number and simulation length.
 * Created by Xinhao Liu on 06/16/20.
 */
public class ABCModelBuilder {
    public static double trueRecombRate = 1.5E-7;

    public static ModelTree getABCModelMyValue() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((A, B), C);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize((int)38508.510267593214);
                node.setNodeHeight(89327.1359517113);
            }
        }
        tree.getRoot().getData().setPopSize((int)24839.588806066255);
        tree.getRoot().setNodeHeight(424145.79262059653);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getABC1Model() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((A, B), C);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(40000);
                node.setNodeHeight(100000);
            }
        }
        tree.getRoot().getData().setPopSize(40000);
        tree.getRoot().setNodeHeight(150000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getABC2Model() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((A, B), C);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(40000);
                node.setNodeHeight(100000);
            }
        }
        tree.getRoot().getData().setPopSize(40000);
        tree.getRoot().setNodeHeight(200000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getABC3Model() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((A, B), C);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(40000);
                node.setNodeHeight(100000);
            }
        }
        tree.getRoot().getData().setPopSize(40000);
        tree.getRoot().setNodeHeight(300000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getABC4Model() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((A, B), C);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(40000);
                node.setNodeHeight(100000);
            }
        }
        tree.getRoot().getData().setPopSize(40000);
        tree.getRoot().setNodeHeight(400000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }


    public static void main(String[] args) {
        ModelTree model = getABC1Model();
        HmmBuilder builder = new HmmBuilder(model.getTree(), model.getRecombRate());
        HmmCore hmm = builder.build();

        System.out.println(hmm.getStates().size());
//        System.out.println(Arrays.toString(hmm.getPi()));
//        System.out.println(Arrays.deepToString(hmm.getA()));
//        for (HiddenState state:hmm.getStates()) {
//            System.out.println("============");
//            System.out.println(state.getIndex());
//            System.out.println(state.getTree().toNewick());
//        }
    }


}
