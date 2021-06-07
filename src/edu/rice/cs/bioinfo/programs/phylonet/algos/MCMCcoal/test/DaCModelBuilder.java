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
 * Build models for testing divide and conquer.
 * Created by Xinhao Liu on 06/12/20.
 */
public class DaCModelBuilder {
    public static double trueRecombRate = 1.5E-7;

    public static ModelTree getHCOModel() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((H,C), O);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(40000);
                node.setNodeHeight(160000);
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

    public static ModelTree getHCOModelInit() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((H,C), O);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(100000);
                node.setNodeHeight(400000);
            }
        }
        tree.getRoot().getData().setPopSize(100000);
        tree.getRoot().setNodeHeight(800000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHCOModelMyValue() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((H,C), O);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize((int) 37000.06571345537);
                node.setNodeHeight(160394.01177411538);
            }
        }
        tree.getRoot().getData().setPopSize((int) 63043.14247326787);
        tree.getRoot().setNodeHeight(663022.8438284714);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHGOModel() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((H,G), O);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(40000);
                node.setNodeHeight(220000);
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

    public static ModelTree getHGOModelInit() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((H,G), O);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(100000);
                node.setNodeHeight(400000);
            }
        }
        tree.getRoot().getData().setPopSize(100000);
        tree.getRoot().setNodeHeight(800000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHGOModelMyValue() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((H,G), O);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize((int)60000);
                node.setNodeHeight(200000);
            }
        }
        tree.getRoot().getData().setPopSize((int)60000); // 60000
        tree.getRoot().setNodeHeight(680000); // 680000
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHCGModel() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((H,C), G);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(40000);
                node.setNodeHeight(160000);
            }
        }
        tree.getRoot().getData().setPopSize(40000);
        tree.getRoot().setNodeHeight(220000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHCGModelInit() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((H,C),G);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(100000);
                node.setNodeHeight(200000);
            }
        }
        tree.getRoot().getData().setPopSize(100000);
        tree.getRoot().setNodeHeight(400000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHCGModelMyValue() {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((H,C),G);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(60000);
                node.setNodeHeight(140000);
            }
        }
        tree.getRoot().getData().setPopSize(40000);
        tree.getRoot().setNodeHeight(220000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static void main(String[] args) {
        ModelTree model = getHGOModelInit();
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
