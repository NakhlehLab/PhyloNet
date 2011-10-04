package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/31/11
 * Time: 5:46 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class RNewichAssignmentsBlockBodyBase<T extends RNewickAssigmentContainer> implements RNewichAssignmentsBlockBody<T>
{
    public final Iterable<T> Assignments;

    public Iterable<T> getAssignments()
    {
       return Assignments;
    }

    RNewichAssignmentsBlockBodyBase(Iterable<T> assigments)
    {
        Assignments =  assigments;
    }
}
