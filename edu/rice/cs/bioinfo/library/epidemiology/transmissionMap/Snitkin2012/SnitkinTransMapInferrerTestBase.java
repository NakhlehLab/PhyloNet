package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import static ch.lambdaj.Lambda.*;
import static org.hamcrest.CoreMatchers.*;

import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.uci.ics.jung.graph.DirectedGraph;
import junit.framework.Assert;
import org.joda.time.LocalDate;
import org.junit.Test;

import java.math.BigDecimal;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/19/12
 * Time: 4:52 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class SnitkinTransMapInferrerTestBase<D extends Comparable<D>>
{
    class Genome
    {
        public final String Sequence;

        public Genome(String sequence)
        {
            Sequence = sequence;
        }
    }

    @Test
    public void testInferMapsNIHOutbreak()
    {
        List<Genome> genomes = Arrays.asList(
                new Genome("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 1
                new Genome("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 1
                new Genome("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 1
                new Genome("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 1
                new Genome("AAAAAAAAAAAAATTTTAAAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 5
                new Genome("AAAAAAAAAAAAATTTTAAAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 3
                new Genome("AAAAAAAAAAAAATTTTAAAAATAAAAAAAAAAAAAAAAAAA"),   // patient 2
                new Genome("AAAAAAAAAAAAAATTTAAAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 1
                new Genome("TTTAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAA"),   // patient 1
                new Genome("TTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 1
                new Genome("TTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTT"),   // patient 8
                new Genome("TTTAAATTTTTTTAAAAAAAAAAAAAAAAAAATTAAAAAAAA"),   // patient 4
                new Genome("TTTAAATTTTTTTAAAAAAATAAATTTAAAAAAAAAAAAAAA"),   // patient 9
                new Genome("TTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 14
                new Genome("TTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 6
                new Genome("TTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 15
                new Genome("TTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 17
                new Genome("TTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 10
                new Genome("TTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 7
                // skip ventilator (if you are following along in the dataset figure in the NIH paper)
                new Genome("TTTTTTTTTTTTTAAAAAAAAAAAAAATTAAAAAAAAAAAAA"),   // patient 16
                new Genome("TTTTTTTTTTTTTAAAAAAAAAAAAAAAATTTAAAAAAAAAA"),   // patient 11
                new Genome("TTTTTTTTTTTTTAAAATTAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 13
                new Genome("TTTTTTTTTTTTTAAAATTAAAAAAAAAAAAAAAAAAAAAAA"),   // patient 12
                new Genome("TTTTTTTTTTTTTAAAATTTTAAAAAAAAAAAAAAAAAAAAA"));  // patient 18

        HashMap<Genome,Integer> genomeToPatient = new HashMap<Genome, Integer>();
        int genomesIndex = 0;
        genomeToPatient.put(genomes.get(genomesIndex++), 1);
        genomeToPatient.put(genomes.get(genomesIndex++), 1);
        genomeToPatient.put(genomes.get(genomesIndex++), 1);
        genomeToPatient.put(genomes.get(genomesIndex++), 1);
        genomeToPatient.put(genomes.get(genomesIndex++), 5);
        genomeToPatient.put(genomes.get(genomesIndex++), 3);
        genomeToPatient.put(genomes.get(genomesIndex++), 2);
        genomeToPatient.put(genomes.get(genomesIndex++), 1);
        genomeToPatient.put(genomes.get(genomesIndex++), 1);
        genomeToPatient.put(genomes.get(genomesIndex++), 1);
        genomeToPatient.put(genomes.get(genomesIndex++), 8);
        genomeToPatient.put(genomes.get(genomesIndex++), 4);
        genomeToPatient.put(genomes.get(genomesIndex++), 9);
        genomeToPatient.put(genomes.get(genomesIndex++), 14);
        genomeToPatient.put(genomes.get(genomesIndex++), 6);
        genomeToPatient.put(genomes.get(genomesIndex++), 15);
        genomeToPatient.put(genomes.get(genomesIndex++), 17);
        genomeToPatient.put(genomes.get(genomesIndex++), 10);
        genomeToPatient.put(genomes.get(genomesIndex++), 7);
        genomeToPatient.put(genomes.get(genomesIndex++), 16);
        genomeToPatient.put(genomes.get(genomesIndex++), 11);
        genomeToPatient.put(genomes.get(genomesIndex++), 13);
        genomeToPatient.put(genomes.get(genomesIndex++), 12);
        genomeToPatient.put(genomes.get(genomesIndex++), 18);

        Func2<Genome,Genome,Double> genDistance = new Func2<Genome,Genome,Double>()
        {
            public Double execute(Genome input1, Genome input2) {
                return getGeneticDistance(input1.Sequence, input2.Sequence);
            }
        };

        Object wardD = new Object();
        Object wardB = new Object();
        Object wardC = new Object();
        Object icu = new Object();
        Object wardF = new Object();

        Map<Integer, Map<LocalDate,Object>> patientTraceData = new HashMap<Integer, Map<LocalDate, Object>>();
        patientTraceData.put(1, new HashMap<LocalDate, Object>());
        patientTraceData.put(2, new HashMap<LocalDate, Object>());
        patientTraceData.put(3, new HashMap<LocalDate, Object>());
        patientTraceData.put(5, new HashMap<LocalDate, Object>());
        patientTraceData.put(6, new HashMap<LocalDate, Object>());
        patientTraceData.put(7, new HashMap<LocalDate, Object>());
        patientTraceData.put(8, new HashMap<LocalDate, Object>());
        patientTraceData.put(9, new HashMap<LocalDate, Object>());
        patientTraceData.put(10, new HashMap<LocalDate, Object>());
        patientTraceData.put(11, new HashMap<LocalDate, Object>());
        patientTraceData.put(12, new HashMap<LocalDate, Object>());
        patientTraceData.put(13, new HashMap<LocalDate, Object>());
        patientTraceData.put(14, new HashMap<LocalDate, Object>());
        patientTraceData.put(15, new HashMap<LocalDate, Object>());
        patientTraceData.put(16, new HashMap<LocalDate, Object>());
        patientTraceData.put(17, new HashMap<LocalDate, Object>());
        patientTraceData.put(18, new HashMap<LocalDate, Object>());

        Map<LocalDate,Object> patient4Trace = new HashMap<LocalDate, Object>();
        patientTraceData.put(4, patient4Trace);
        addRange(wardD, new LocalDate(2011,6,23), new LocalDate(2011,8,11), patient4Trace);
        addRange(icu,   new LocalDate(2011,8,12), new LocalDate(2011,8,13), patient4Trace);
        addRange(wardD, new LocalDate(2011,8,14), new LocalDate(2011,8,29), patient4Trace);
        addRange(icu,   new LocalDate(2011,8,30), new LocalDate(2011,9,2),  patient4Trace);


        Map<Integer,LocalDate> patientToFirstPositive = new HashMap<Integer, LocalDate>();
        patientToFirstPositive.put(1, new LocalDate(2011,6,14));
        patientToFirstPositive.put(2, new LocalDate(2011,8,5));
        patientToFirstPositive.put(3, new LocalDate(2011,8,15));
        patientToFirstPositive.put(4, new LocalDate(2011,8,23));
   //     patientToFirstPositive.put(5, new LocalDate(2011,8,29));
        patientToFirstPositive.put(6, new LocalDate(2011,9,15));
     //   patientToFirstPositive.put(7, new LocalDate(2011,9,19));
        patientToFirstPositive.put(8, new LocalDate(2011,9,22));
        patientToFirstPositive.put(9, new LocalDate(2011,9,26));
     //   patientToFirstPositive.put(10, new LocalDate(2011,10,6));
        patientToFirstPositive.put(11, new LocalDate(2011,10,17));
     //   patientToFirstPositive.put(12, new LocalDate(2011,10,17));
        patientToFirstPositive.put(13, new LocalDate(2011,11,3));
     /*   patientToFirstPositive.put(14, new LocalDate(2011,11,10));
        patientToFirstPositive.put(15, new LocalDate(2011,11,17));  */
        patientToFirstPositive.put(16, new LocalDate(2011,11,18));
   /*     patientToFirstPositive.put(17, new LocalDate(2011,11,27));
        patientToFirstPositive.put(18, new LocalDate(2011,12,14)); */

        SnitkinTransMapInferrer<Integer,Genome,D> transMapInferrer = makeInferrer(1, genomes, genomeToPatient, genDistance, patientTraceData, patientToFirstPositive);
        DirectedGraph<Integer,SnitkinEdge<Integer,D>> transMap = transMapInferrer.inferMaps().iterator().next();

        Assert.assertEquals(1, select(transMap.getEdges(), having(on(SnitkinEdge.class).getSource(), equalTo(1)).and(   // edge 1 -> 8
                having(on(SnitkinEdge.class).getDestination(), equalTo(8)))).size());

        Assert.assertEquals(1, select(transMap.getEdges(), having(on(SnitkinEdge.class).getSource(),      equalTo(1)).and(   // edge 1 -> 3
                                                           having(on(SnitkinEdge.class).getDestination(), equalTo(3)))).size());

    }

    protected abstract SnitkinTransMapInferrer<Integer,Genome,D> makeInferrer(Integer rootPatient, List<Genome> genomes, HashMap<Genome, Integer> genomeToPatient,
                                                                              Func2<Genome, Genome, Double> genDistance, Map<Integer, Map<LocalDate, Object>> patientTraceData, Map<Integer, LocalDate> patientToFirstPositive);

    private void addRange(Object location, LocalDate first, LocalDate last, Map<LocalDate, Object> trace) {

        for(LocalDate i = first; i.isBefore(last.plusDays(1)); i = i.plusDays(1))
        {
            trace.put(i, location);
        }

    }

    private double getGeneticDistance(String genome1, String genome2)
    {
        if(genome1.length() != genome2.length())
        {
            throw new IllegalArgumentException("Genomes must have same length");
        }

        double variantsAccum = 0;
        for(int i = 0; i<genome1.length(); i++)
        {
            if(genome1.charAt(i) != genome2.charAt(i))
                variantsAccum++;
        }

        return variantsAccum;
    }
}
