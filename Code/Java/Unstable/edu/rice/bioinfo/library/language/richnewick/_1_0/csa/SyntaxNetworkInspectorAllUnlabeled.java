package edu.rice.bioinfo.library.language.richnewick._1_0.csa;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/25/11
 * Time: 2:15 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class SyntaxNetworkInspectorAllUnlabeled<T> implements SyntaxNetworkInspector<T>
{
    public String getNodeLabelText(T node)
    {
        return null;
    }

    public int getNodeLabelTextLineNumber(T node)
    {
        return -1;
    }

    public int getNodeLabelTextColumnNumber(T node)
    {
        return -1;
    }
}
