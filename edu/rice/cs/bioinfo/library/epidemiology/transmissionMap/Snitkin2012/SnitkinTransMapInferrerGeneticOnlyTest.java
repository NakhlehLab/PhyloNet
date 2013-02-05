package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import junit.framework.Assert;
import org.junit.Test;

import java.util.Collection;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/25/13
 * Time: 5:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class SnitkinTransMapInferrerGeneticOnlyTest extends NIHOutbreakDataTestBase
{


    @Test
    public void testInferMap()
    {

        SnitkinTransMapInferrerTemplate<Integer,NIHOutbreakDataTestBase.Sequencing,Integer> transMapInferrer =
                new SnitkinTransMapInferrerGeneticOnly(this.sequenceToPatient,

                        new Func2<Sequencing,Sequencing,Integer>()
                        {
                            public Integer execute(Sequencing s1, Sequencing s2) {
                                return s1.getGeneticDistance(s2);
                            }
                        },

                        new Func2<Integer,Integer,Integer>()
                        {

                            public Integer execute(Integer input1, Integer input2) {
                                return input1 + input2;
                            }
                        },
                        new Func2<Integer,Integer,Integer>()
                        {

                            public Integer execute(Integer input1, Integer input2) {
                                return input1 - input2;
                            }
                        },
                        new Func1Identity<Integer>(),
                        new Func<Integer>()
                        {
                            public Integer execute() {
                                return 100000000;
                            }
                        });
        Set<SnitkinEdge<Integer,Integer>> patientGraph = transMapInferrer.makeCompletePatientGraph(1);
        Set<SnitkinEdge<Integer,Integer>> transMapEdges = transMapInferrer.inferMaps(1).iterator().next();
        Collection<SnitkinEdge<Integer,Integer>> paperEdges = IterableHelp.filter(patientGraph, new Predicate1<SnitkinEdge<Integer, Integer>>() {
            public boolean execute(SnitkinEdge<Integer, Integer> someEdge) {

                Integer source = someEdge.getSource();
                Integer dest = someEdge.getDestination();

                return (source.equals(3) && dest.equals(2)) ||
                        (source.equals(1) && dest.equals(3)) ||
                        (source.equals(1) && dest.equals(4)) ||
                        (source.equals(3) && dest.equals(5)) ||
                        (source.equals(4) && dest.equals(6)) ||
                        (source.equals(6) && dest.equals(7)) ||
                        (source.equals(1) && dest.equals(8)) ||
                        (source.equals(4) && dest.equals(9)) ||
                        (source.equals(6) && dest.equals(10)) ||
                        (source.equals(6) && dest.equals(11)) ||
                        (source.equals(13) && dest.equals(12)) ||
                        (source.equals(6) && dest.equals(13)) ||
                        (source.equals(6) && dest.equals(14)) ||
                        (source.equals(6) && dest.equals(15)) ||
                        (source.equals(6) && dest.equals(16)) ||
                        (source.equals(6) && dest.equals(17)) ||
                        (source.equals(13) && dest.equals(18));
            }
        });

        Integer myScore =  score(transMapEdges);
        Integer paperScore = score(paperEdges);
        Assert.assertEquals(myScore, paperScore);
    }

    protected Integer score(Collection<SnitkinEdge<Integer, Integer>> edges) {

        int accum = 0;
        for(SnitkinEdge<Integer,Integer> edge : edges)
        {
            accum += edge.getDistance();
        }

        return accum;
    }
}
