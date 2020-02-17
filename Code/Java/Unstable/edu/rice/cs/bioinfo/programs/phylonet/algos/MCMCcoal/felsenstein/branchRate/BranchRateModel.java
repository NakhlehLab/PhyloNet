package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.branchRate;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.core.StateNode;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HiddenState;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;

/**
 * Created by wendingqiao on 5/4/16.
 * Defines a mean rate for each branch in the tree
 */
public interface BranchRateModel {

    public double getRateForBranch(HiddenState tree, TNode node);

    public abstract class Base extends StateNode implements BranchRateModel {
        protected double _meanRate = 1.0;
    }
}