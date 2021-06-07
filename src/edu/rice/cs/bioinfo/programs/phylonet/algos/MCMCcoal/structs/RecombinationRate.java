package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.core.StateNode;

import java.math.BigDecimal;

/**
 * Recombination rate of the species tree
 * This is one model parameter (others are population size hyper parameter and species tree)
 * Created by Xinhao Liu on 11/8/19
 */
public class RecombinationRate extends StateNode {

    private double _recombRate;

    public RecombinationRate(double rate) {
        this._recombRate = rate;
    }

    public double getRecombRate() {
        return this._recombRate;
    }

    public void setRecombRate(double newRate) {
        this._recombRate = newRate;
    }

    public String toString() {
        return BigDecimal.valueOf(_recombRate).toPlainString();
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
        return false;
    }
}
