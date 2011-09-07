package edu.rice.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 2:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class PhyloNetCommand implements PySONNode {

    public Iterable<PhyloNetCommandPart> Parts;

    public PhyloNetCommand(Iterable<PhyloNetCommandPart> parts)
    {
        Parts = parts;
    }

}
