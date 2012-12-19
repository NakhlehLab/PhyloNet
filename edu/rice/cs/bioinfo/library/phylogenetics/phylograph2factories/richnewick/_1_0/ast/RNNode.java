package edu.rice.cs.bioinfo.library.phylogenetics.phylograph2factories.richnewick._1_0.ast;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/1/12
 * Time: 2:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class RNNode
{
    public final String Label;

    public final HybridNodeType HybridType;

    public RNNode(String label, HybridNodeType hybridType)
    {
        Label = label;
        HybridType = hybridType;
    }

     public RNNode(String label)
    {
        this(label, null);
    }

    @Override
    public boolean equals(Object obj)
    {
        try
        {
            return equals((RNNode)obj);
        }
        catch(ClassCastException e)
        {
            return false;
        }
    }

    public boolean equals(RNNode other)
    {
        if(other.Label == null)
        {
            if(this.Label == null)
                return other == this;
            else
                return false;
        }
        else
        {
            if(this.Label == null)
                return false;
            else
                return other.Label.equals(this.Label);
        }

    }

    @Override
    public String toString()
    {
        return Label == null ? super.toString() : Label.toString();
    }

    @Override
    public int hashCode()
    {
        return Label == null ? super.hashCode() : Label.hashCode();
    }
}