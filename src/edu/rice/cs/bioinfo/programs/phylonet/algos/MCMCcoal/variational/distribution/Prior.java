package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.TreeNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * Prior of the entire model
 *
 * Created by Xinhao Liu on 3/14/20.
 */
public class Prior {
    private NormalDistribution recombRatePrior;
    private GammaDistribution popSizePrior;

    public Prior() {
        this(new NormalDistribution(3, 1));
    }

    public Prior(NormalDistribution recombRatePrior) {
        this.recombRatePrior = recombRatePrior;
        this.popSizePrior = new GammaDistribution(Utils.GAMMA_SHAPE, Utils.POP_SIZE_MEAN / Utils.GAMMA_SHAPE);
    }

    public double logPrior(ModelTree model) {
        double prior = getNodeHeightPrior(model) + getPopSizePrior(model) + getRecombRatePrior(model);

        if (Utils.DISABLE_ALL_PRIOR) {
            prior = 0.0;
        }

        return prior;
    }

    private double getNodeHeightPrior(ModelTree model) {
        /* Uniform distribution */
        return 0.0;
    }

    private double getPopSizePrior(ModelTree model) {
        //return 0.0;
        double logP = 0.0;
        for (TNode node:model.getTree().postTraverse()) {
            if (!node.isLeaf()) {
                STINode<TreeNodeInfo> stiNode = (STINode<TreeNodeInfo>) node;
                int popSize = stiNode.getData().getPopSize();
                logP += popSizePrior.logDensity(popSize);
                if (popSize <= 0) {
                    return 0.0;
                }
            }
        }
        return logP;
    }

    private double getRecombRatePrior(ModelTree model) {
        //return recombRatePrior.logDensity(model.getRecombRate().getRecombRate() / Utils.RECOMB_RATE_SCALE);
        return 0.0;
    }
}
