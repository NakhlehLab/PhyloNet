package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.branchRate;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.core.StateNode;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;

/**
 * Created by wendingqiao on 5/4/16.
 * Defines a mean rate for each branch in the tree
 */
public interface BranchRateModel {

    public double getRateForBranch(UltrametricTree tree, TNode node);

    public abstract class Base extends StateNode implements BranchRateModel {
        protected double _meanRate = 1.0;
    }
}