package edu.rice.cs.bioinfo.library.phylogenetics;

/**
* Created by IntelliJ IDEA.
* User: Matt
* Date: 5/13/12
* Time: 3:09 PM
* To change this template use File | Settings | File Templates.
*/
public class PhyloEdge<T> extends PhyloEdge2<T,Double>
{
    public PhyloEdge(T source, T destination)
    {
        super(source, destination);
    }

    public PhyloEdge(T source, T destination, double branchLength)
    {
        super(source,destination,branchLength);
    }

    public PhyloEdge<T> setBranchLength(double newBranchLength)
    {
        super.setBranchLength(newBranchLength);
        return this;
    }
    
    public int hashCode(){
        return this.Source.hashCode() + this.Destination.hashCode();
    }


    public PhyloEdge<T> setSupport(double newSupport)
    {
        super.setSupport(newSupport);
        return this;
    }

    public PhyloEdge<T> setProbability(double newProbability)
    {
        super.setProbability(newProbability);
        return this;
    }
}
