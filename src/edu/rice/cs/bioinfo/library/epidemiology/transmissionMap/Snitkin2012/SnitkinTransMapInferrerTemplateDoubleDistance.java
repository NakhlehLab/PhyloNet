package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import org.joda.time.LocalDate;

import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/20/12
 * Time: 11:42 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class SnitkinTransMapInferrerTemplateDoubleDistance<E, S> extends SnitkinTransMapInferrerTemplate<E, S,Double>
{

    public SnitkinTransMapInferrerTemplateDoubleDistance(Map<E, Map<LocalDate, Object>> patientTraces,
                                                         Map<E, LocalDate> patientToFirstPositiveDate, Map<S, E> sequencingToPatient) {
        super(patientTraces, patientToFirstPositiveDate, sequencingToPatient);
    }

    @Override
    protected Double add(Double a, Double b) {
        return a + b;
    }

    @Override
    protected Double subtract(Double a, Double b) {
        return a - b;
    }

    @Override
    protected Double makeDistance(int intDistance) {
        return new Double(intDistance);
    }

    @Override
    protected SnitkinEdge<E, Double> makeEdge(E source, E destination, int geneticDistance, int epidemiologicalDistance) {
        return new SnitkinEdgeDouble<E>(source, destination, geneticDistance, epidemiologicalDistance, this.getEMax());
    }
}
