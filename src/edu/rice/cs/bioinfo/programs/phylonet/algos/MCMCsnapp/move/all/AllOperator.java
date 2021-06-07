package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.all;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;

import java.util.List;

/**
 * Created by wendingqiao on 4/29/16.
 * Abstract operator that changes trees and network
 */
public abstract class AllOperator extends Operator {

    protected List<UltrametricTree> _allTrees;
    protected UltrametricNetwork _speciesNet;
    protected boolean _violate; // if the operation would violate the temporal constraint

    public AllOperator(List<UltrametricTree> trees, UltrametricNetwork net)  {
        this._allTrees = trees;
        this._speciesNet = net;
        this._violate = false;
    }

    public Utils.MOVE_TYPE getCategory() {
        return Utils.MOVE_TYPE.ALL;
    }

    public boolean mayViolate() {
        return _violate;
    }

    public void optimize (double logAlpha) {}

}

