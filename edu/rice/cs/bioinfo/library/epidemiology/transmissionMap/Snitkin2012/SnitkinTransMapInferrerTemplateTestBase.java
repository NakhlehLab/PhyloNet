package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import static ch.lambdaj.Lambda.*;
import static org.hamcrest.CoreMatchers.*;

import edu.rice.cs.bioinfo.library.language.dot_2013_1.printing.DigraphDotPrinter;
import edu.rice.cs.bioinfo.library.programming.Predicate1;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import junit.framework.Assert;
import junit.framework.AssertionFailedError;
import org.joda.time.LocalDate;
import org.junit.Test;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/19/12
 * Time: 4:52 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class SnitkinTransMapInferrerTemplateTestBase<D extends Comparable<D>> extends NIHOutbreakDataTestBase
{
    @Test
    public void testInferMapsFromEdgesDynamic()
    {

    }

    @Test
    public void testInferMapsNIHOutbreak()
    {

        SnitkinTransMapInferrerTemplate<Integer,Sequencing,D> transMapInferrer = makeInferrer(this.sequenceToPatient, patientTraceData, patientIdToFirstPositive);


        Set<SnitkinEdge<Integer,D>> patientGraph = transMapInferrer.makeCompletePatientGraph(1);
        Set<SnitkinEdge<Integer,D>> transMapEdges = transMapInferrer.inferMaps(1).iterator().next();
        Collection<SnitkinEdge<Integer,D>> paperEdges = IterableHelp.filter(patientGraph, new Predicate1<SnitkinEdge<Integer, D>>() {
            public boolean execute(SnitkinEdge<Integer, D> someEdge) {

                Integer source = someEdge.getSource();
                Integer dest = someEdge.getDestination();

                return (source.equals(5) && dest.equals(2)) ||
                        (source.equals(1) && dest.equals(3)) ||
                        (source.equals(1) && dest.equals(4)) ||
                        (source.equals(3) && dest.equals(5)) ||
                        (source.equals(10) && dest.equals(6)) ||
                        (source.equals(6) && dest.equals(7)) ||
                        (source.equals(1) && dest.equals(8)) ||
                        (source.equals(4) && dest.equals(9)) ||
                        (source.equals(4) && dest.equals(10)) ||
                        (source.equals(10) && dest.equals(11)) ||
                        (source.equals(7) && dest.equals(12)) ||
                        (source.equals(12) && dest.equals(13)) ||
                        (source.equals(17) && dest.equals(14)) ||
                        (source.equals(14) && dest.equals(15)) ||
                        (source.equals(15) && dest.equals(16)) ||
                        (source.equals(7) && dest.equals(17)) ||
                        (source.equals(13) && dest.equals(18));
            }
        });

        verifyEdgesEpiDistance(transMapEdges, transMapInferrer.getEMax());
        verifyEdgesEpiDistance(paperEdges, transMapInferrer.getEMax());

        String mySolutionDot = new DigraphDotPrinter<Integer,SnitkinEdge<Integer,D>>()
        {
            @Override
            protected Integer getSource(SnitkinEdge<Integer, D> edge) {
                return edge.getSource();
            }

            @Override
            protected Integer getDestination(SnitkinEdge<Integer, D> edge) {
                return edge.getDestination();
            }

            @Override
            protected String getNodeLabel(Integer node)
            {
                return node.toString();
            }

            @Override
            protected String getEdgeLabel(SnitkinEdge<Integer, D> edge)
            {
                return "GD = " + edge.getGeneticDistance() + ", ED = " + edge.getEpidemiologicalDistance();
            }
        }.toDot(new HashSet<SnitkinEdge<Integer, D>>( paperEdges));

        int myGeneticScore = scoreGenetics(transMapEdges);
        int paperGeneticScore = scoreGenetics(paperEdges);

        int myEpiScore = scoreEpiDistance(transMapEdges);
        int paperEpiScore = scoreEpiDistance(paperEdges);

        D myScore =  score(transMapEdges);
        D paperScore = score(paperEdges);
        Assert.assertEquals(myScore, paperScore);


    }

    private void verifyEdgesEpiDistance(Collection<SnitkinEdge<Integer,D>> edges, int eMax)
    {
        Map<Integer,Map<Integer,Integer>> epiDistanceMap = this.makeEpiDistanceMap(eMax);

        for(SnitkinEdge<Integer,D> edge : edges)
        {
            int expected = epiDistanceMap.get(edge.getSource()).get(edge.getDestination()).intValue();
            int found = edge.getEpidemiologicalDistance();
            //  try
            //  {
            Assert.assertEquals(expected, found);
            //  }
            //   catch(AssertionFailedError e)
            //   {
            //       throw e;
            //   }
        }
    }

    private int scoreEpiDistance(Collection<SnitkinEdge<Integer, D>> edges) {
        int score = 0;
        for(SnitkinEdge<Integer,D> edge : edges)
        {
            score += edge.getEpidemiologicalDistance();
        }

        return score;
    }

    private int scoreGenetics(Collection<SnitkinEdge<Integer, D>> edges)
    {
        int score = 0;
        for(SnitkinEdge<Integer,D> edge : edges)
        {
            score += edge.getGeneticDistance();
        }

        return score;
    }



    protected abstract D score(Collection<SnitkinEdge<Integer,D>> edges);


    protected abstract SnitkinTransMapInferrerTemplate<Integer,Sequencing,D> makeInferrer(HashMap<Sequencing, Integer> genomeToPatient, Map<Integer, Map<LocalDate, Object>> patientTraceData, Map<Integer, LocalDate> patientToFirstPositive);




}
