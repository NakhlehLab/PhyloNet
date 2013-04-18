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
public class SnitkinTransMapInferrerGeneticOnly<E,S,D extends Comparable<D>> extends SnitkinTransMapInferrer<E,S,D>
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

    public SnitkinTransMapInferrerGeneticOnly(Map<S, E> genomeToPatient, Func2<S, S, Integer> getGeneticDistance, Func2<D, D, D> add, Func2<D, D, D> subtract,
                                              final Func1<Integer, D> makeDistance, final Func<D> getMaxDistance)
    {
        super(new HashMap<E, Map<LocalDate, Object>>(), makePatientToFirstPositiveDate(new HashSet<E>(genomeToPatient.values())), genomeToPatient, getGeneticDistance, add, subtract, makeDistance, getMaxDistance,

               new Func4<E, E, Integer, Integer, SnitkinEdge<E, D>>()
               {

                   public SnitkinEdge<E, D> execute(E donor, E recipient, final Integer geneticDistance, Integer epidemiologicalDistance) {

                       return new SnitkinEdgeBase<E, D>(donor, recipient, geneticDistance, 0) {
                           public D getDistance() {
                               return makeDistance.execute(geneticDistance);
                           }
                       };
                   }
               });

    }

    @Override
    protected int getEpidemiologicalDistance(E donor, E recipient)
    {
        return 0;
    }

}
