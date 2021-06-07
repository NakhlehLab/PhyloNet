package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Predicate1;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import junit.framework.Assert;
import org.junit.Test;

import java.util.Collection;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/22/13
 * Time: 7:36 PM
 * To change this template use File | Settings | File Templates.
 */
public class SnitkinTransMapInferrerEpiOnlyTest extends NIHOutbreakDataTestBase
{
    @Test
    public void testInferMap()
    {

        SnitkinTransMapInferrerTemplate<Integer,NIHOutbreakDataTestBase.Sequencing,Double> transMapInferrer =
                new SnitkinTransMapInferrerEpiOnly(this.patientTraceData, this.patientIdToFirstPositive, this.sequenceToPatient,
                        new Func2<Double,Double,Double>()
                        {

                            public Double execute(Double input1, Double input2) {
                                return input1 + input2;
                            }
                        },
                        new Func2<Double,Double,Double>()
                        {

                            public Double execute(Double input1, Double input2) {
                                return input1 - input2;
                            }
                        },
                        new Func1<Integer,Double>()
                        {
                            public Double execute(Integer input) {
                                return new Double(input);
                            }
                        },
                        new Func<Double>()
                        {

                            public Double execute() {
                                return 100000.0;
                            }
                        });
        Set<SnitkinEdge<Integer,Double>> patientGraph = transMapInferrer.makeCompletePatientGraph(1);
        Set<SnitkinEdge<Integer,Double>> epiMapEdges = transMapInferrer.inferMaps(1).getSolutions().iterator().next();
        Collection<SnitkinEdge<Integer,Double>> paperEdges = IterableHelp.filter(patientGraph, new Predicate1<SnitkinEdge<Integer, Double>>() {
            public boolean execute(SnitkinEdge<Integer, Double> someEdge) {

                Integer source = someEdge.getSource();
                Integer dest = someEdge.getDestination();

                return (source.equals(5) && dest.equals(2)) ||
                        (source.equals(1) && dest.equals(3)) ||
                        (source.equals(2) && dest.equals(4)) ||
                        (source.equals(3) && dest.equals(5)) ||
                        (source.equals(10) && dest.equals(6)) ||
                        (source.equals(6) && dest.equals(7)) ||
                        (source.equals(11) && dest.equals(8)) ||
                        (source.equals(12) && dest.equals(9)) ||
                        (source.equals(9) && dest.equals(10)) ||
                        (source.equals(9) && dest.equals(11)) ||
                        (source.equals(3) && dest.equals(12)) ||
                        (source.equals(14) && dest.equals(13)) ||
                        (source.equals(11) && dest.equals(14)) ||
                        (source.equals(9) && dest.equals(15)) ||
                        (source.equals(15) && dest.equals(16)) ||
                        (source.equals(14) && dest.equals(17)) ||
                        (source.equals(17) && dest.equals(18));
            }
        });

        Double score =  score(epiMapEdges);
        Double paperScore = score(paperEdges);
        Assert.assertEquals(paperScore, score);
        Assert.assertEquals(17, epiMapEdges.size());
        findEdge(5,2,epiMapEdges);
        findEdge(1,3,epiMapEdges);
        findEdge(2,4,epiMapEdges);
        findEdge(3,5,epiMapEdges);
        findEdge(10,6,epiMapEdges);
        findEdge(6,7,epiMapEdges);
        findEdge(11,8,epiMapEdges);
        findEdge(12,9,epiMapEdges);
        findEdge(9,10,epiMapEdges);
        findEdge(9,11,epiMapEdges);
        findEdge(3,12,epiMapEdges);
        findEdge(14,13,epiMapEdges);
        findEdge(11,14,epiMapEdges);
        findEdge(9,15,epiMapEdges);
        findEdge(15,16,epiMapEdges);
        findEdge(14,17,epiMapEdges);
        findEdge(17,18,epiMapEdges);

    }

    protected Double score(Collection<SnitkinEdge<Integer, Double>> edges) {

        double accum = 0;
        for(SnitkinEdge<Integer,Double> edge : edges)
        {
            accum += edge.getDistance();
        }

        return accum;
    }

    private SnitkinEdge<Integer, Double> findEdge(final int parent, final int child, Set<SnitkinEdge<Integer, Double>> edges) {

        return IterableHelp.filter(edges, new Predicate1<SnitkinEdge<Integer, Double>>() {
            public boolean execute(SnitkinEdge<Integer, Double> edge) {
                return edge.getSource().intValue() == parent && edge.getDestination().intValue() == child;
            }
        }).iterator().next();

    }
}
