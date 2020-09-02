package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.TreeNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.NodeHeightVariable;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.PopSizeVariable;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.RecombRateVariable;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.VariationalVariable;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Variational approximation to the posterior of the entire model
 *
 * Created by Xinhao Liu on 3/14/20.
 */
public class VariationalModel {
    private ModelTree model;

    private List<VariationalVariable> nodeHeightVariableList;

    private List<VariationalVariable> popSizeVariableList;

    private VariationalVariable recombRateVariable;

    private VariationalVariable popSizeHyperParameter; // not needed anymore

    // Two variational variables need to access the same tree node (pop size of branch and node height need to be stored in one tree node). If some tree node property
    // is modified in one variational variable, the tree node in the other variational variable is also changed.

    public VariationalModel(ModelTree model) {
        nodeHeightVariableList = new ArrayList<>();
        popSizeVariableList = new ArrayList<>();

        this.model = model;
        STITree<TreeNodeInfo> tree = model.getTree();

        for (TNode n:tree.postTraverse()) {
            if (n.isLeaf()) {
                continue;
            }

            STINode<TreeNodeInfo> node = (STINode<TreeNodeInfo>) n;
            nodeHeightVariableList.add(new NodeHeightVariable(node.getNodeHeight(), Utils.NODE_HEIGHT_INIT_STDDEV, node));
            popSizeVariableList.add(new PopSizeVariable(node.getData().getPopSize(), Utils.POP_SIZE_INIT_STDDEV, node));
//            // for rerunning experiments
//            if (!node.isRoot()) {
//                if (node.getLeafCount() == 2) {
//                    nodeHeightVariableList.add(new NodeHeightVariable(node.getNodeHeight(), 10011.247868849025, node, 0.033828001710504006, 0.025428084436187555));
//                    popSizeVariableList.add(new PopSizeVariable(node.getData().getPopSize(), 3016.599224635356, node, 0.19164748035022744, 0.3607854600907661));
//                } else if (node.getLeafCount() == 3) {
//                    nodeHeightVariableList.add(new NodeHeightVariable(node.getNodeHeight(), 10003.181075402601, node, 0.17554798359388676, 0.054143807247309134));
//                    popSizeVariableList.add(new PopSizeVariable(node.getData().getPopSize(), 3003.9218789065776, node, 0.45660731841409624, 0.2763947980863776));
//                }
//            } else {
//                nodeHeightVariableList.add(new NodeHeightVariable(node.getNodeHeight(), 21222.721879009703, node, 0.003344919972527076, 0.0051103078990821635));
//                popSizeVariableList.add(new PopSizeVariable(node.getData().getPopSize(), 8104.702014262722, node, 0.025133714761294428, 0.06397856909265269));
//            }
        }

//        recombRateVariable = new RecombRateVariable(model.getRecombRate().getRecombRate() / Utils.RECOMB_RATE_SCALE, Utils.RECOMB_RATE_INIT_STDDEV, model.getRecombRate());
    }

    /**
     * Sample from the variational posterior.
     */
    public Map<VariationalVariable, Double> sample() {
        // sample each variable in turn
        Map<VariationalVariable, Double> sample = new HashMap<>();
        for (VariationalVariable var:nodeHeightVariableList) {
            sample.put(var, var.sample());
        }
        for (VariationalVariable var:popSizeVariableList) {
            sample.put(var, var.sample());
        }
        //sample.put(recombRateVariable, recombRateVariable.sample());
//        // BELOW IS DEBUG CODE
//        System.out.println("samples for node heights:");
//        for (VariationalVariable var:nodeHeightVariableList) {
//            System.out.println(sample.get(var));
//        }
//        System.out.println("samples for pop sizes:");
//        for (VariationalVariable var:popSizeVariableList) {
//            System.out.println(sample.get(var));
//        }
//        System.out.println("sample for recomb rate:");
//        System.out.println(sample.get(recombRateVariable));
//        System.out.println("VariationalModel: finished generating one sample");
        return sample;
    }

    public void setTreeBySample(Map<VariationalVariable, Double> sample) {
        for (Map.Entry<VariationalVariable, Double> entry: sample.entrySet()) {
            VariationalVariable var = entry.getKey();
            double value = entry.getValue();
            var.setVariableValue(value);
        }
        model.refresh();
    }

    public ModelTree getModel() {
        return model;
    }

    public List<VariationalVariable> getNodeHeightVariableList() {
        return nodeHeightVariableList;
    }

    public List<VariationalVariable> getPopSizeVariableList() {
        return popSizeVariableList;
    }

    public VariationalVariable getRecombRateVariable() {
        return recombRateVariable;
    }
}
