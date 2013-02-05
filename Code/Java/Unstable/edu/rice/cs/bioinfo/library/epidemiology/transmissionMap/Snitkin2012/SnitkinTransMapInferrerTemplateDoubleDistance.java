package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import org.joda.time.LocalDate;

import java.math.BigDecimal;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/20/12
 * Time: 11:42 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class SnitkinTransMapInferrerTemplateDoubleDistance<P, S> extends SnitkinTransMapInferrerTemplate<P, S,Double>
{

    public SnitkinTransMapInferrerTemplateDoubleDistance(Map<P, Map<LocalDate, Object>> patientTraces,
                                                         Map<P, LocalDate> patientToFirstPositiveDate, Map<S, P> sequencingToPatient) {
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
    protected SnitkinEdge<P, Double> makeEdge(P source, P destination, int geneticDistance, int epidemiologicalDistance) {
        return new SnitkinEdgeDouble<P>(source, destination, geneticDistance, epidemiologicalDistance, this.getEMax());
    }
}
