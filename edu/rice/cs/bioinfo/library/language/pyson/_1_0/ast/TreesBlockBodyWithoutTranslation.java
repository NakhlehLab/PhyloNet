package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/31/11
 * Time: 5:46 PM
 * To change this template use File | Settings | File Templates.
 */
public class TreesBlockBodyWithoutTranslation extends RNewichAssignmentsBlockBodyBase<TreeAssignment> implements TreesBlockBody {

    public TreesBlockBodyWithoutTranslation(Iterable<TreeAssignment> assignments)
    {
        super(assignments);
    }

    public <R, T, E extends Exception> R execute(BlockAlgo<R, T, E> algo, T input) throws E {
        return algo.forTreesBlock(this, input);
    }
}
