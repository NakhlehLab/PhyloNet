package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.core;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Randomizer;

import java.util.Arrays;

/**
 * Created by Xinhao Liu on 11/1/19.
 * Modified from edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core.StateNode.java
 */
public abstract class StateNode {
    protected Operator _operator = null;
    protected boolean _dirty = false;

    // propose and undo are for moved node
    public abstract double propose();
    public abstract void undo();

    // accept and reject are for dirty node
    public abstract void accept();
    public abstract void reject();

    public abstract double logDensity();

    public abstract boolean mayViolate();

    public abstract boolean isValid();

    public Operator getOperation() {
        return _operator;
    }

    public void setDirty(boolean dirty) {
        this._dirty = dirty;
    }

    public boolean isDirty() {
        return _dirty;
    }

    public static Operator getOp(Operator[] operators, double[] weights) {
        if(operators.length < weights.length) {
            throw new RuntimeException("Operator-Weight pair doesn't match " + Arrays.toString(weights));
        }
        double rand = Randomizer.getRandomDouble();
        for(int i = 0; i < weights.length; i++) {
            if(rand < weights[i]) {
                return operators[i];
            }
        }
        throw new RuntimeException("Error propose Operation " + rand); // should never reach here
    }
}
