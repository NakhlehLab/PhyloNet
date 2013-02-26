package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import org.joda.time.LocalDate;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/20/12
 * Time: 11:55 AM
 * To change this template use File | Settings | File Templates.
 */
public class SnitkinTransMapInferrerTemplateDoubleDistanceTest extends SnitkinTransMapInferrerTemplateTestBase<Double> {


    @Override
    protected Double score(Collection<SnitkinEdge<Integer, Double>> edges) {

        double accum = 0;
        for(SnitkinEdge<Integer,Double> edge : edges)
        {
            accum += edge.getDistance();
        }

        return accum;
    }

    @Override
    protected SnitkinTransMapInferrerTemplate<Integer, Sequencing, Double> makeInferrer(HashMap<Sequencing, Integer> genomeToPatient, Map<Integer, Map<LocalDate, Object>> patientTraceData, Map<Integer, LocalDate> patientToFirstPositive) {
         return new SnitkinTransMapInferrerTemplateDoubleDistance<Integer,Sequencing>(patientTraceData, patientToFirstPositive, genomeToPatient)
         {
             @Override
             protected Double getMaxDistance() {
                 return 1000000.0;
             }

             @Override
             protected int getGeneticDistance(Sequencing sequence1, Sequencing sequence2) {
                 return sequence1.getGeneticDistance(sequence2);
             }


         };

    }
}
