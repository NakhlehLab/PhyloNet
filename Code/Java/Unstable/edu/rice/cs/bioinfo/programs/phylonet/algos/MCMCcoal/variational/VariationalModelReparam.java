package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.TreeNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.BranchLengthVariable;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.NodeHeightVariable;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.PopSizeVariable;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.VariationalVariable;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Variational approximation to the posterior of the reparametrized model.
 * Reparametrized model parameters are branch lengths and population sizes, instead of node heights and population sizes.
 *
 * Created by Xinhao Liu on 1/22/21.
 */
public class VariationalModelReparam {
    private ModelTree model;

    private List<VariationalVariable> branchLengthVariableList;

    private List<VariationalVariable> popSizeVariableList;

    // Two variational variables need to access the same tree node (pop size of branch and node height need to be stored in one tree node). If some tree node property
    // is modified in one variational variable, the tree node in the other variational variable is also changed.

    public VariationalModelReparam(ModelTree model) {
        branchLengthVariableList = new ArrayList<>();
        popSizeVariableList = new ArrayList<>();

        this.model = model;
        STITree<TreeNodeInfo> tree = model.getTree();

        for (TNode n:tree.postTraverse()) {
            if (n.isLeaf()) {
                continue;
            }

            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) n;
            // leaf nodes does not incur branch lengths
            // if internal node has two leaf children, incur a branch length of the branch connecting the node to an arbitrary leaf
            // if internal node has one leaf child, incur a branch length of the branch connecting the node to the non-leaf child
            // if internal node has non leaf child, incur a branch length of the branch connecting the node to the left child
            if (node.countLeafChildren() == 2) {
                branchLengthVariableList.add(new BranchLengthVariable(node.getNodeHeight(), Utils.BRANCH_LENGTH_INIT_STDDEV, node, node.getChildren().iterator().next()));
            } else if (node.countLeafChildren() == 1) {
                STINode<TreeNodeInfo> nonleafChild = null;
                for (STINode<TreeNodeInfo> child : node.getChildren()) {
                    if (!child.isLeaf()) {
                        nonleafChild = child;
                        break;
                    }
                }
                branchLengthVariableList.add(new BranchLengthVariable(nonleafChild.getParentDistance(), Utils.BRANCH_LENGTH_INIT_STDDEV, node, nonleafChild));
            } else if (node.countLeafChildren() == 0) {
                // TODO: how to find which one is the left child? for now choose the first one in children list
                STINode<TreeNodeInfo> leftChild = null;
                for (STINode<TreeNodeInfo> child : node.getChildren()) {
                    leftChild = child;
                    break;
                }
                branchLengthVariableList.add(new BranchLengthVariable(leftChild.getParentDistance(), Utils.BRANCH_LENGTH_INIT_STDDEV, node, leftChild));
            } else {
                System.out.println("Shouldn't be here.");
                System.exit(1);
            }

            popSizeVariableList.add(new PopSizeVariable(node.getData().getPopSize(), Utils.POP_SIZE_INIT_STDDEV, node));
        }
    }

    /**
     * Sample from the variational posterior.
     */
    public Map<VariationalVariable, Double> sample() {
        // sample each variable in turn
        Map<VariationalVariable, Double> sample = new HashMap<>();
        for (VariationalVariable var:branchLengthVariableList) {
            sample.put(var, var.sample());
        }
        for (VariationalVariable var:popSizeVariableList) {
            sample.put(var, var.sample());
        }
        return sample;
    }

    public void setTreeBySample(Map<VariationalVariable, Double> sample) {
//        // BELOW IS DEBUG CODE
//        System.out.println("samples for branch lengths:");
//        for (VariationalVariable var:branchLengthVariableList) {
//            System.out.println(sample.get(var));
//        }
//        System.out.println("samples for pop sizes:");
//        for (VariationalVariable var:popSizeVariableList) {
//            System.out.println(sample.get(var));
//        }
//        // ABOVE IS DEBUG CODE
        for (Map.Entry<VariationalVariable, Double> entry: sample.entrySet()) {
            VariationalVariable var = entry.getKey();
            double value = entry.getValue();
            var.setVariableValue(value);
        }
        model.refreshNewBranchLength();
//        // BELOW IS DEBUG CODE
//        STITree<TreeNodeInfo> tree = model.getTree();
//        System.out.println(tree.toNewick());
//        for (TNode node:tree.postTraverse()) {
//            STINode<TreeNodeInfo> stiNode = (STINode<TreeNodeInfo>) node;
//            System.out.println(stiNode.getNodeHeight());
//            System.out.println(stiNode.getData().getPopSize());
//            System.out.println("---");
//        }
//        // ABOVE IS DEBUG CODE
    }

    public ModelTree getModel() {
        return model;
    }

    public List<VariationalVariable> getBranchLengthVariableList() {
        return branchLengthVariableList;
    }

    public List<VariationalVariable> getPopSizeVariableList() {
        return popSizeVariableList;
    }
}
