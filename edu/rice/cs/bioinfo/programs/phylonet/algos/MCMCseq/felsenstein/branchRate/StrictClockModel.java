package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.branchRate;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;

/**
 * Created by wendingqiao on 5/4/16.
 */
public class StrictClockModel extends BranchRateModel.Base {

    private double _upperBound = Double.MAX_VALUE;
    private double _lowerBound = 0.0;

    public StrictClockModel() {
        _lowerBound = _meanRate;
        _upperBound = _meanRate;
    }

    public StrictClockModel(double mu) {
        _meanRate = mu;
        _upperBound = mu;
        _lowerBound = mu;
    }

    public StrictClockModel(double mu, double muUpperBound, double muLowerBound) {
        _meanRate = mu;
        _upperBound = Math.max(mu, muUpperBound);
        _lowerBound = Math.max(0, muLowerBound);
    }

    @Override
    public double getRateForBranch(final UltrametricTree tree, final TNode node) {
        return _meanRate;
    }

    @Override
    public double propose() {
        return 0;
    }

    @Override
    public void undo() {

    }

    @Override
    public void accept() {

    }

    @Override
    public void reject() {

    }

    @Override
    public double logDensity() {
        return 0;
    }

    @Override
    public boolean mayViolate() {
        return false;
    }

    @Override
    public boolean isValid() {
        return true;
    }
}