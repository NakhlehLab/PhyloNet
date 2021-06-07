package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/19/12
 * Time: 1:59 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class SnitkinEdgeBase<P,D> implements SnitkinEdge<P,D> {

    private final P _source;
    private final P _destination;
    private final int _geneticDistance;
    private final int _epidemiologicalDistance;

    public SnitkinEdgeBase(P source, P destination, int geneticDistance, int epidemiologicalDistance) {
        _source = source;
        _destination = destination;
        _geneticDistance = geneticDistance;
        _epidemiologicalDistance = epidemiologicalDistance;
    }


    public P getSource() {
        return _source;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public P getDestination() {
        return _destination;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getGeneticDistance() {
        return _geneticDistance;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getEpidemiologicalDistance() {
        return _epidemiologicalDistance;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String toString()
    {
        return getSource().toString() + "->" + getDestination().toString();
    }
}
