package edu.rice.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/31/11
 * Time: 5:46 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class TreesBlockBodyBase implements TreesBlockBody
{
    public final Iterable<TreeAssignment> Assigments;

    TreesBlockBodyBase(Iterable<TreeAssignment> assigments)
    {
        Assigments =  assigments;
    }
}
