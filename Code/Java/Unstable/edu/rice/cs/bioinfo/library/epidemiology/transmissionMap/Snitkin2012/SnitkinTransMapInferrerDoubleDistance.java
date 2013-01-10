package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import edu.rice.cs.bioinfo.library.programming.Func2;
import org.joda.time.LocalDate;

import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/20/12
 * Time: 11:42 AM
 * To change this template use File | Settings | File Templates.
 */
public class SnitkinTransMapInferrerDoubleDistance<P, S> extends SnitkinTransMapInferrer<P, S,Double>
{

    public SnitkinTransMapInferrerDoubleDistance(P rootPatient, List<S> genomes, Func2<S, S, Double> getGeneticDistance, Map<P, Map<LocalDate, Object>> patientTraces,
                                                 Map<P, LocalDate> patientToFirstPositiveDate, Map<S, P> genomeToPatient) {
        super(rootPatient, genomes, getGeneticDistance, patientTraces, patientToFirstPositiveDate, genomeToPatient);
    }

    @Override
    protected Double subtract(Double a, Double b) {
        return a - b;
    }

    @Override
    protected Double makeZero() {
        return 0.0;
    }

    @Override
    protected Double makeMax() {
        return Double.MAX_VALUE;
    }

    @Override
    protected Double makeDistance(int intDistance) {
        return new Double(intDistance);
    }

    @Override
    protected SnitkinEdge<P, Double> makeEdge(P source, P destination, Double geneticDistance, Double epidemiologicalDistance) {
        return new SnitkinEdgeDouble<P>(source, destination, geneticDistance, epidemiologicalDistance);
    }
}
