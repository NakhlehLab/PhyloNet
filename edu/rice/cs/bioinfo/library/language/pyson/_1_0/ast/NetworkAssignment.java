package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/20/11
 * Time: 11:17 AM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkAssignment implements PySONNode, RNewickAssigmentContainer
{
    public final RichNewickAssignment Assignment;

    public NetworkAssignment(RichNewickAssignment assigment)
    {
        Assignment = assigment;
    }

    public RichNewickAssignment getAssignment() {
        return Assignment;
    }
}