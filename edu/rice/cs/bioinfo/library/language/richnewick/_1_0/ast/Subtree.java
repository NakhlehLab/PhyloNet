package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/14/11
 * Time: 6:47 PM
 * To change this template use File | Settings | File Templates.
 */
public class Subtree implements AbstractSyntaxNode
{
    public final DescendantList Descendants;

    public final NetworkInfo NetworkInfo;

    public Subtree(DescendantList descendants, NetworkInfo nodeInfo)
    {
        Descendants = descendants;
        NetworkInfo = nodeInfo;
    }

}
