package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;

/**
 * Created by wendingqiao on 2/15/16.
 */
public abstract class Operator {

    public int _acCounter = 0;
    public int _rejCounter = 0;
    public int _acCorrectionCounter = 0;
    public int _rejCorrectionCounter = 0;
    public Utils.Transform _transform = Utils.Transform.None;

    /**
     * make a new move
     * @return hastings ratio
     */
    public abstract double propose ();

    /**
     * undo the move
     */
    public abstract void undo ();

    /**
     * called after every invocation of this operator to see whether a parameter
     * can be optimised for better acceptance hence faster mixing
     *
     * @param logAlpha difference in posterior between previous state & proposed state + hasting ratio
     */
    public abstract void optimize (double logAlpha);

    /**
     * @return  The category of operation
     */
    public abstract Utils.MOVE_TYPE getCategory();

    /**
     * @return  If the operator may violate the temporal constraints
     */
    public abstract boolean mayViolate();

    /**
     * @return  The name of the operator
     */
    public abstract String getName();
}
