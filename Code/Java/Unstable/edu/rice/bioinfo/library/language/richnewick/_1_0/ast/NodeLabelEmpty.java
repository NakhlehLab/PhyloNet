package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/14/11
 * Time: 6:01 PM
 * To change this template use File | Settings | File Templates.
 */
public final class NodeLabelEmpty implements NodeLabel
{
    public static final NodeLabelEmpty Singleton = new NodeLabelEmpty();

    private  NodeLabelEmpty()
    {
    }

    public <R, T, E extends Exception> R execute(NodeLabelAlgo<R, T, E> algo, T input) throws E {
        return algo.forNodeLabelEmpty(this, input);
    }
}
