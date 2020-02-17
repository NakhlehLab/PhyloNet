package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.move;

/**
 * Created by Xinhao Liu on 11/1/19.
 * Modified from edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.Operator
 */
public abstract class Operator {
    public int _acCounter = 0;
    public int _rejCounter = 0;
    public int _acCorrectionCounter = 0;
    public int _rejCorrectionCounter = 0;
    //TODO: what is this?
    //public Utils.Transform _transform = Utils.Transform.None;

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
    //TODO: what is this?
    public abstract void optimize (double logAlpha);

    /**
     * @return  The category of operation
     */
    //public abstract Utils.MOVE_TYPE getCategory();

    /**
     * @return  If the operator may violate the temporal constraints
     */
    //TODO: what is this?
    public abstract boolean mayViolate();

    /**
     * @return  The name of the operator
     */
    public abstract String getName();
}
