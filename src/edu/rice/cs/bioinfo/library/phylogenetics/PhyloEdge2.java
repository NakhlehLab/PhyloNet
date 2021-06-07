package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/28/12
 * Time: 11:43 AM
 * To change this template use File | Settings | File Templates.
 */
public class PhyloEdge2<T,D> {

    public static <T,D> Func1<PhyloEdge2<T,D>, D> GetBranchLength()
    {
        return new Func1<PhyloEdge2<T, D>, D>() {
            public D execute(PhyloEdge2<T, D> edge) {
                return  edge.getBranchLength();
            }
        };
    }

    private D _branchLength;

    public D getBranchLength()
    {
        return _branchLength;
    }

    public PhyloEdge2<T,D> setBranchLength(D newBranchLength)
    {
        _branchLength = newBranchLength;
        return this;
    }

    private D _support;

    public D getSupport()
    {
        return _support;
    }

    public PhyloEdge2<T,D> setSupport(D newSupport)
    {
        _support = newSupport;
        return this;
    }

    private D _probability;

    public D getProbability()
    {
        return _probability;
    }

    public PhyloEdge2<T,D> setProbability(D newProbability)
    {
        _probability = newProbability;
        return this;
    }

    public final T Source;

    public final T Destination;

    public final Tuple<T,T> NodesOfEdge;

    public PhyloEdge2(T source, T destination)
    {
        Source = source;
        Destination = destination;
        NodesOfEdge = new Tuple<T, T>(source,destination);
    }

    public PhyloEdge2(T source, T destination, D branchLength)
    {
        this(source,destination);
        setBranchLength(branchLength);
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
               candidate.Destination.equals(this.Destination);
               //candidate.getBranchLength() == this.getBranchLength() &&
               //candidate.getProbability() == this.getProbability() &&
               //candidate.getSupport() == this.getSupport();
    }


    @Override
    public String toString()
    {
        return "(" + Source + "," + Destination + ")";
    }

}
