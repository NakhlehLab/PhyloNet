package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Func4;
import org.joda.time.LocalDate;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/25/13
 * Time: 4:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class SnitkinTransMapInferrerGeneticOnly<P,S,D extends Comparable<D>> extends SnitkinTransMapInferrer<P,S,D>
{
    private static <P> Map<P, LocalDate> makePatientToFirstPositiveDate(Set<P> patients)
    {
        Map<P, LocalDate> patientToFirstPositiveDate = new HashMap<P, LocalDate>();

        for(P patient : patients)
        {
            patientToFirstPositiveDate.put(patient, null);
        }

        return patientToFirstPositiveDate;
    }

    public SnitkinTransMapInferrerGeneticOnly(Map<S, P> genomeToPatient, Func2<S, S, Integer> getGeneticDistance, Func2<D, D, D> add, Func2<D, D, D> subtract,
                                              final Func1<Integer, D> makeDistance, final Func<D> getMaxDistance)
    {
        super(new HashMap<P, Map<LocalDate, Object>>(), makePatientToFirstPositiveDate(new HashSet<P>(genomeToPatient.values())), genomeToPatient, getGeneticDistance, add, subtract, makeDistance, getMaxDistance,

               new Func4<P, P, Integer, Integer, SnitkinEdge<P, D>>()
               {

                   public SnitkinEdge<P, D> execute(P donor, P recipient, final Integer geneticDistance, Integer epidemiologicalDistance) {

                       return new SnitkinEdgeBase<P, D>(donor, recipient, geneticDistance, 0) {
                           public D getDistance() {
                               return makeDistance.execute(geneticDistance);
                           }
                       };
                   }
               });

    }

    @Override
    protected int getEpidemiologicalDistance(P donor, P recipient)
    {
        return 0;
    }

}
