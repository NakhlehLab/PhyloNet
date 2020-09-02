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
 * Build HCGO models for testing purposes.
 *
 * Created by Xinhao Liu on 06/05/20.
 */
public class HCGOModelBuilder {
    public static double trueRecombRate = 1.5E-7;

    public static ModelTree getHCGOModel() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("(((H, C), G), O);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (TNode tnode:tree.postTraverse()) {
            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) tnode;
            if (!node.isLeaf() && !node.isRoot()) {
                if (node.getLeafCount() == 2) {
                    node.getData().setPopSize(40000);
                    node.setNodeHeight(160000);
                } else if (node.getLeafCount() == 3) {
                    node.getData().setPopSize(40000);
                    node.setNodeHeight(220000);
                }
            }
        }
        tree.getRoot().getData().setPopSize(40000);
        tree.getRoot().setNodeHeight(720000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHCGOModelMyValue() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("(((H, C), G), O);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (TNode tnode:tree.postTraverse()) {
            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) tnode;
            if (!node.isLeaf() && !node.isRoot()) {
                if (node.getLeafCount() == 2) {
                    node.getData().setPopSize((int)65000);
                    node.setNodeHeight(135000);
                } else if (node.getLeafCount() == 3) {
                    node.getData().setPopSize((int)40000);
                    node.setNodeHeight(220000);
                }
            }
        }
        tree.getRoot().getData().setPopSize(40000);
        tree.getRoot().setNodeHeight(720000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHCGOModelMyValue2() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("(((H, C), G), O);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (TNode tnode:tree.postTraverse()) {
            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) tnode;
            if (!node.isLeaf() && !node.isRoot()) {
                if (node.getLeafCount() == 2) {
                    node.getData().setPopSize((int)61017.31114306007);
                    node.setNodeHeight(141094.52870456775);
                } else if (node.getLeafCount() == 3) {
                    node.getData().setPopSize((int)42609.40174792459);
                    node.setNodeHeight(204221.17285133447);
                }
            }
        }
        tree.getRoot().getData().setPopSize((int)67281.02212425758);
        tree.getRoot().setNodeHeight(673312.6340551431);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHCGOModelInit() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("(((H, C), G), O);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (TNode tnode:tree.postTraverse()) {
            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) tnode;
            if (!node.isLeaf() && !node.isRoot()) {
                if (node.getLeafCount() == 2) {
                    node.getData().setPopSize(50000);
                    node.setNodeHeight(180000);
                } else if (node.getLeafCount() == 3) {
                    node.getData().setPopSize(50000);
                    node.setNodeHeight(360000);
                }
            }
        }
        tree.getRoot().getData().setPopSize(50000);
        tree.getRoot().setNodeHeight(720000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static void main(String[] args) {
        ModelTree model = getHCGOModelMyValue2();
        HmmBuilder builder = new HmmBuilder(model.getTree(), model.getRecombRate());
        HmmCore hmm = builder.build();
//
//        System.out.println(builder.getCoalescentHistory("(4:20.598,(1:6.581,(2:5.690,3:5.690):0.890):14.018);"));

        System.out.println(hmm.getStates().size());
//        System.out.println(Arrays.toString(hmm.getPi()));
////        System.out.println(Arrays.deepToString(hmm.getA()));
        for (HiddenState state:hmm.getStates()) {
            System.out.println("============");
//            System.out.println(state.getIndex());
            System.out.println(state.getTree().toNewick());
        }
    }
}
