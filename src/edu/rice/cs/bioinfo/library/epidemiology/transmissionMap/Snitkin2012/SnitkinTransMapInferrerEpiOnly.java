package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Func4;
import org.joda.time.LocalDate;

import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/22/13
 * Time: 7:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class SnitkinTransMapInferrerEpiOnly<P,S,D extends Comparable<D>> extends SnitkinTransMapInferrer<P,S,D>
{
    public SnitkinTransMapInferrerEpiOnly(Map<P, Map<LocalDate, Object>> patientTraces, Map<P, LocalDate> patientToFirstPositiveDate, Map<S, P> genomeToPatient,
                                          Func2<D, D, D> add, Func2<D, D, D> subtract, final Func1<Integer, D> makeDistance, final Func<D> makeMaxDistance)
    {
        super(patientTraces, patientToFirstPositiveDate, genomeToPatient,

                new Func2<S, S, Integer>() {
                    public Integer execute(S input1, S input2) {
                        return 0;
                    }}

                , add, subtract, makeDistance, makeMaxDistance,

                new Func4<P, P, Integer, Integer, SnitkinEdge<P, D>>() {
                    public SnitkinEdge<P, D> execute(P donor, P recipient, Integer geneticDistance, final Integer epidemiologicalDistance)
                    {
                        return new SnitkinEdgeBase<P, D>(donor, recipient, geneticDistance, epidemiologicalDistance) {

                            public D getDistance() {
                                return makeDistance.execute(epidemiologicalDistance);
                            }
                        };
                    }
                });
    }
}
