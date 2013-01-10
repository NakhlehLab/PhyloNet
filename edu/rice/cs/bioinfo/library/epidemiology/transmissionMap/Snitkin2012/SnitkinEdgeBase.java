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
    private final D _geneticDistance;
    private final D _epidemiologicalDistance;

    public SnitkinEdgeBase(P source, P destination, D geneticDistance, D epidemiologicalDistance) {
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

    public D getGeneticDistance() {
        return _geneticDistance;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public D getEpidemiologicalDistance() {
        return _epidemiologicalDistance;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
