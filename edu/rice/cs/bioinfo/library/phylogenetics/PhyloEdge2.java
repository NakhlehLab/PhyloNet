package edu.rice.cs.bioinfo.library.phylogenetics;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/28/12
 * Time: 11:43 AM
 * To change this template use File | Settings | File Templates.
 */
public class PhyloEdge2<T,D> {

    private D _branchLength;

    public D getBranchLength()
    {
        return _branchLength;
    }

    public void setBranchLength(D newBranchLength)
    {
        _branchLength = newBranchLength;
    }

    private D _support;

    public D getSupport()
    {
        return _support;
    }

    public void setSupport(D newSupport)
    {
        _support = newSupport;
    }

    private D _probability;

    public D getProbabilty()
    {
        return _probability;
    }

    public void setProbability(D newProbability)
    {
        _probability = newProbability;
    }

    public final T Source;

    public final T Destination;

    public PhyloEdge2(T source, T destination)
    {
        Source = source;
        Destination = destination;
    }

    public boolean equals(Object candidate)
    {
        try
        {
            return equals((PhyloEdge2)candidate);
        }
        catch (ClassCastException e)
        {
            return false;
        }
    }

    public boolean equals(PhyloEdge2 candidate)
    {
        return candidate.Source.equals(this.Source) &&
               candidate.Destination.equals(this.Destination) &&
               candidate.getBranchLength() == this.getBranchLength() &&
               candidate.getProbabilty() == this.getProbabilty() &&
               candidate.getSupport() == this.getSupport();
    }

    @Override
    public String toString()
    {
        return "(" + Source + "," + Destination + ")";
    }

}
