package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.all;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;

import java.util.List;

/**
 * Created by wendingqiao on 4/29/16.
 * Abstract operator that changes trees and network
 */
public abstract class AllOperator extends Operator {

    protected UltrametricNetwork _speciesNet;
    protected boolean _violate; // if the operation would violate the temporal constraint

    public AllOperator(UltrametricNetwork net)  {
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

