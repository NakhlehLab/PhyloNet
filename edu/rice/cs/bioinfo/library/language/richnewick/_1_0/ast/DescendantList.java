package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/14/11
 * Time: 6:44 PM
 * To change this template use File | Settings | File Templates.
 */
public class DescendantList implements AbstractSyntaxNode
{
    public static final DescendantList EMPTY_DESCENDANT_LIST = new DescendantList(new ArrayList<Subtree>());

    public final Iterable<Subtree> Subtrees;

    public DescendantList(Iterable<Subtree> subtrees)
    {
        Subtrees = subtrees;
    }
}
