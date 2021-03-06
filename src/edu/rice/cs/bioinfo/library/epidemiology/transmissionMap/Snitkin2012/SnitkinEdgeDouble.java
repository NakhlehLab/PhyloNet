package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/19/12
 * Time: 1:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class SnitkinEdgeDouble<P> extends SnitkinEdgeBase<P,Double> {

    private final Double _distance;

    public SnitkinEdgeDouble(P source, P destination, int geneticDistance, int epidemiologicalDistance, int eMax) {
        super(source, destination, geneticDistance, epidemiologicalDistance);


        int epiDistance = this.getEpidemiologicalDistance();
        _distance = geneticDistance + (( 999.0 * ( ((double)epiDistance) / ((double)eMax))) / 100000.0);
    }

    public Double getDistance() {
        return _distance;
    }
}
