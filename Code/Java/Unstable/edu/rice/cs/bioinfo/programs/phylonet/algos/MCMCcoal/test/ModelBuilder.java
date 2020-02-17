package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.RecombinationRate;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.TreeNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

/**
 * Build different models for testing purposes.
 * Created by Xinhao Liu on 02/08/20.
 */
public class ModelBuilder {
    public static ModelTree getHCGModel() {
        RecombinationRate dummyRecombRate = new RecombinationRate(1.5E-8);
        ModelTree model = new ModelTree("((H,C), G);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(30000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(40000);
                node.setNodeHeight(160000);
            }
        }
        tree.getRoot().getData().setPopSize(40000);
//        tree.getRoot().getData().setPopSize(400000);
        tree.getRoot().setNodeHeight(220000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHCGModelBad() {
        // Set HC and HCG population size to 60000, change HC divergence time to 100000 gens
        RecombinationRate dummyRecombRate = new RecombinationRate(1.5E-8);
        ModelTree model = new ModelTree("((H,C), G);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(30000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(60000);
                node.setNodeHeight(100000);
            }
        }
        tree.getRoot().getData().setPopSize(60000);
//        tree.getRoot().getData().setPopSize(400000);
        tree.getRoot().setNodeHeight(400000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHCGModelWrongRecomb(double recombRate) {
        RecombinationRate dummyRecombRate = new RecombinationRate(recombRate);
        ModelTree model = new ModelTree("((H,C), G);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(30000);
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

    public static ModelTree getHCGModelWrongTopo() {
        // Set topology to H(CG)
        RecombinationRate dummyRecombRate = new RecombinationRate(1.5E-8);
        ModelTree model = new ModelTree("(H,(C,G));", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(30000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(60000);
                node.setNodeHeight(100000);
            }
        }
        tree.getRoot().getData().setPopSize(60000);
//        tree.getRoot().getData().setPopSize(400000);
        tree.getRoot().setNodeHeight(400000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHCGModelSlightlyWrong() {
        // change HC population size to
        RecombinationRate dummyRecombRate = new RecombinationRate(1.5E-8);
        ModelTree model = new ModelTree("((H,C), G);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(30000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(500000);
                node.setNodeHeight(160000);
            }
        }
        tree.getRoot().getData().setPopSize(40000);
//        tree.getRoot().getData().setPopSize(400000);
        tree.getRoot().setNodeHeight(220000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHCGModelWrongHCSize(Integer HCSize) {
        // change HC population size to
        RecombinationRate dummyRecombRate = new RecombinationRate(1.5E-8);
        ModelTree model = new ModelTree("((H,C), G);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(30000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(HCSize);
                node.setNodeHeight(160000);
            }
        }
        tree.getRoot().getData().setPopSize(40000);
//        tree.getRoot().getData().setPopSize(400000);
        tree.getRoot().setNodeHeight(220000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }


}
