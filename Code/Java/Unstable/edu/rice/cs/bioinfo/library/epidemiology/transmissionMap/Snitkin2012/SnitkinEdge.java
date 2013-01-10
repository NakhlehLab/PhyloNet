package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import edu.rice.cs.bioinfo.library.programming.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/19/12
 * Time: 12:58 PM
 * To change this template use File | Settings | File Templates.
 */
public interface SnitkinEdge<P,D>
{
    public P getSource();

    public P getDestination();

    public D getGeneticDistance();

    public D getEpidemiologicalDistance();

    public D getDistance();
}
