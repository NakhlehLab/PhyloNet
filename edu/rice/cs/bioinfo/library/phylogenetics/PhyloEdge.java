package edu.rice.cs.bioinfo.library.phylogenetics;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 3:09 PM
 * To change this template use File | Settings | File Templates.
 */
public class PhyloEdge<T>
{
    private double _branchLength = 1.0;

    public double getBranchLength()
    {
        return _branchLength;
    }

    public void setBranchLength(double newBranchLength)
    {
        _branchLength = newBranchLength;
    }

    private double _support = 100;

    public double getSupport()
    {
        return _support;
    }

    public void setSupport(double newSupport)
    {
        if(newSupport <= 100 && newSupport >= 0.0)
        {
            _support = newSupport;
        }
        else
        {
            throw new IllegalArgumentException("Support value must be between zero and a hundred inclusive.  Found: " + newSupport);
        }
    }

    private double _probability = 1.0;

    public double getProbabilty()
    {
        return _probability;
    }

    public PhyloEdge<T> setProbability(double newProbability)
    {
        if(newProbability <= 1.0 && newProbability >= 0.0)
        {
            _probability = newProbability;
        }
        else
        {
            throw new IllegalArgumentException("Probability value must be between zero and one inclusive.  Found: " + newProbability);
        }

        return this;
    }

    public final T Source;

    public final T Destination;

    public PhyloEdge(T source, T destination)
    {
        Source = source;
        Destination = destination;
    }

    public PhyloEdge(T source, T destination, double branchLength)
    {
        this(source,destination);
        setBranchLength(branchLength);
    }

    public boolean equals(Object candidate)
    {
        try
        {
            return equals((PhyloEdge)candidate);
        }
        catch (ClassCastException e)
        {
            return false;
        }
    }

    public boolean equals(PhyloEdge candidate)
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
