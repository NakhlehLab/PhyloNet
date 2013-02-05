package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

public interface SnitkinEdge<P,D>
{
    public P getSource();

    public P getDestination();

    public int getGeneticDistance();

    public int getEpidemiologicalDistance();

    public D getDistance();
}
