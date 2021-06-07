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
 * Build HCG models for testing purposes.
 * Created by Xinhao Liu on 02/08/20.
 */
public class HCGModelBuilder {
//    public static double trueRecombRate = 3.0E-7; //TODO: adjust this
//    public static double trueRecombRate = 1.5E-7;
    public static double trueRecombRate = 1.5E-8;
    //public static double trueRecombRate = 2E-9;

    public static ModelTree getHCGModel() {
        //RecombinationRate dummyRecombRate = new RecombinationRate(1.5E-8);
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

    public static ModelTree getHCGModelSetValue(double t_hc, double t_hcg, double n_hc, double n_hcg) {
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((H,C), G);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize((int) n_hc);
                node.setNodeHeight(t_hc);
            }
        }
        tree.getRoot().getData().setPopSize((int) n_hcg);
        tree.getRoot().setNodeHeight(t_hcg);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHCGModelBad() {
        // Set HC and HCG population size to 60000, change HC divergence time to 100000 gens
        //RecombinationRate dummyRecombRate = new RecombinationRate(1.5E-8);
        RecombinationRate dummyRecombRate = new RecombinationRate(3.0E-7); //TODO: adjust this
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
//        RecombinationRate dummyRecombRate = new RecombinationRate(1.5E-8);
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
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

    public static ModelTree getHCGModelIllegal() {
        //RecombinationRate dummyRecombRate = new RecombinationRate(1.5E-8);
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((H,C), G);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(30000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(-40000);
                node.setNodeHeight(240000);
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
        //RecombinationRate dummyRecombRate = new RecombinationRate(1.5E-8);
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((H,C), G);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize(50000);
                node.setNodeHeight(150000);
            }
        }
        tree.getRoot().getData().setPopSize(50000);
        tree.getRoot().setNodeHeight(300000);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static ModelTree getHCGModelMyValue() {
        //RecombinationRate dummyRecombRate = new RecombinationRate(1.5E-8);
        RecombinationRate dummyRecombRate = new RecombinationRate(trueRecombRate);
        ModelTree model = new ModelTree("((H,C), G);", dummyRecombRate);
        STITree<TreeNodeInfo> tree = model.getTree();
        for (String leafName:tree.getLeaves()) {
            STINode<TreeNodeInfo> leafNode = tree.getNode(leafName);
            leafNode.getData().setPopSize(10000);
        }
        for (STINode<TreeNodeInfo> node:tree.getNodes()) {
            if (!node.isLeaf() && !node.isRoot()) {
                node.getData().setPopSize((int)66096.92358712363);
                node.setNodeHeight(155712.35168420564);
            }
        }
        tree.getRoot().getData().setPopSize((int)39253.06129897994);
        tree.getRoot().setNodeHeight(214100.1565971283);
        for (TNode node:tree.postTraverse()) {
            if (!node.isRoot()) {
                node.setParentDistance(node.getParent().getNodeHeight() - node.getNodeHeight());
            }
        }
        return model;
    }

    public static void main(String[] args) {
        ModelTree model = getHCGModel();
        HmmBuilder builder = new HmmBuilder(model.getTree(), model.getRecombRate());
        HmmCore hmm = builder.build();

        System.out.println(hmm.getStates().size());
        System.out.println(Arrays.toString(hmm.getPi()));
        System.out.println(Arrays.deepToString(hmm.getA()));
        for (HiddenState state:hmm.getStates()) {
            System.out.println("============");
            System.out.println(state.getIndex());
            System.out.println(state.getTree().toNewick());
        }
    }



}
