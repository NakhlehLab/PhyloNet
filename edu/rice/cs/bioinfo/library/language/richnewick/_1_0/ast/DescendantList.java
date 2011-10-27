package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

import com.sun.xml.internal.ws.api.model.Parameter;

import java.util.ArrayList;
import java.util.LinkedList;

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

    public DescendantList(NetworkInfo... leafs)
    {
        LinkedList<Subtree> trees = new LinkedList<Subtree>();

        for(NetworkInfo n : leafs)
        {
            trees.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, n));
        }

        Subtrees = trees;
    }
}
