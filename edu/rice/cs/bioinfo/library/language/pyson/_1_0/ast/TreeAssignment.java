package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/31/11
 * Time: 6:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class TreeAssignment implements PySONNode
{
    public final boolean IsRooted;

    public final RichNewickAssignment Assignment;

    public TreeAssignment(boolean isRooted, RichNewickAssignment assigment)
    {
        IsRooted = isRooted;
        Assignment = assigment;
    }
}
